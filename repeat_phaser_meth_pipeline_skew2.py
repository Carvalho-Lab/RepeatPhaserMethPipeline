#!/usr/bin/env python3
# Author: Chris Seward
# Script Name: RepeatPhaserMethPipeline.py
# Version: 3.0
# Prototyping assistance provided by Google Gemini 2.5 Pro AI
# Date: June 2, 2025 (Updated)
# Description: This script processes multiple BAM files to:
#              1. Phase long reads based on repeat element copy numbers.
#              2. Tag reads in an output BAM file with haplotype information (HP:i:0, 1, or 2).
#              3. Run modbamtools plot and calcMeth for methylation analysis.
#              4. Perform differential methylation testing (Welch's t-test) between determined haplotypes.
#              5. Calculate methylation percent change skewness.
#              6. Output phased BAM and analysis files into structured directories per input BAM.
#              7. Generate a consolidated summary report for all processed BAMs with detailed phasing stats.

# --- Standard library imports ---
import argparse
import datetime
import os
import re
import subprocess
import sys

# --- Third-party library imports ---
# These must be installed via pip: pip install pysam numpy scikit-learn scipy pandas
import numpy as np
import pandas as pd
import pysam
from scipy.stats import mannwhitneyu, chi2, ttest_ind_from_stats
from sklearn.cluster import KMeans

# --- Function Definitions ---

def parse_genomic_region(region_str):
    """
    Parses a genomic region string (e.g., "chr1:100-200") into a standardized format.

    Args:
        region_str (str): Genomic region string in "chr:start-end" format (1-based inclusive).

    Returns:
        tuple: (chromosome, 0-based_start, 1-based_inclusive_end) or None if parsing fails.
               The end coordinate is suitable for both BED format and as the exclusive end for pysam's 0-based fetching.
    """
    try:
        chrom_part, range_part = region_str.split(':', 1)
        start_str, end_str = range_part.split('-')
        start_1based = int(start_str)
        end_1based = int(end_str) 
        if start_1based <= 0 or end_1based < start_1based:
            raise ValueError("Invalid region coordinates.")
        return chrom_part, start_1based - 1, end_1based 
    except ValueError as e:
        sys.stderr.write(f"Error parsing region string '{region_str}': {e}\n")
        return None

def sanitize_region_for_filename(region_str):
    """ 
    Sanitizes a region string (e.g., "chrX:100-200") to a format suitable for filenames 
    (e.g., "chrX_100_200").
    """
    return re.sub(r'[:-]', '_', region_str)

def get_aligned_read_segment_in_region(read, region_chrom, region_start_0, region_end_0_exclusive):
    """
    Extracts the portion of the read's sequence that aligns within a specified reference region.
    This function correctly handles insertions and deletions in the read's CIGAR string.

    Args:
        read (pysam.AlignedSegment): The read object.
        region_chrom (str): Chromosome of the region.
        region_start_0 (int): 0-based start coordinate of the region.
        region_end_0_exclusive (int): 0-based exclusive end coordinate of the region.

    Returns:
        tuple: (aligned_segment_string, count_of_aligned_reference_bases_in_region)
    """
    if read.reference_name != region_chrom:
        return "", 0
    aligned_segment_chars = []
    aligned_ref_bp_count = 0
    # get_aligned_pairs maps each base in the read to its corresponding reference position.
    for q_idx, r_idx in read.get_aligned_pairs(matches_only=False, with_seq=False):
        if r_idx is not None and region_start_0 <= r_idx < region_end_0_exclusive:
            aligned_ref_bp_count +=1 
            if q_idx is not None: # If q_idx is None, it's a deletion in the read.
                aligned_segment_chars.append(read.query_sequence[q_idx])
    return "".join(aligned_segment_chars), aligned_ref_bp_count

def count_repeats_in_read_segment(read, rep_chrom, rep_start_0, rep_end_0_exclusive, 
                                  motif, min_overlap_bp):
    """
    Counts non-overlapping occurrences of a motif within the part of the read aligned to the RE region.

    Args:
        read (pysam.AlignedSegment): The read object.
        rep_chrom (str): Chromosome of the repetitive element.
        rep_start_0 (int): 0-based start of the repetitive element region.
        rep_end_0_exclusive (int): 0-based exclusive end of the repetitive element region.
        motif (str): The motif sequence to count (e.g., "CAG").
        min_overlap_bp (int): Minimum number of base pairs the read must align within the RE region.

    Returns:
        int or np.nan: The count of the motif, or np.nan if criteria are not met.
    """
    if read.is_unmapped or read.reference_name != rep_chrom:
        return np.nan
    read_ref_start, read_ref_end = read.reference_start, read.reference_end
    if not (max(read_ref_start, rep_start_0) < min(read_ref_end, rep_end_0_exclusive)):
        return np.nan # No overlap
    
    aligned_segment, aligned_ref_bases = get_aligned_read_segment_in_region(
        read, rep_chrom, rep_start_0, rep_end_0_exclusive
    )
    if aligned_ref_bases < min_overlap_bp or not aligned_segment:
        return np.nan
    return aligned_segment.upper().count(motif.upper())

def run_modbamtools_plot(phased_bam_path_abs, plot_region_str, analysis_dir, plot_basename_prefix):
    """
    Calls the external `modbamtools plot` command to generate methylation plots.
    The command is run from within the specified analysis_dir to keep outputs organized.
    """
    sys.stdout.write(f"\n  Running modbamtools plot for region {plot_region_str}...\n")
    cmd = ["modbamtools", "plot", phased_bam_path_abs, "-r", plot_region_str, "-o", ".", "-p", plot_basename_prefix, "--hap"]
    sys.stdout.write(f"    Executing in '{analysis_dir}': {' '.join(cmd)}\n")
    try:
        process = subprocess.run(cmd, cwd=analysis_dir, check=True, capture_output=True, text=True)
        sys.stdout.write("    modbamtools plot completed successfully.\n")
        if process.stdout.strip(): sys.stdout.write("    Plot stdout:\n" + process.stdout + "\n")
        if process.stderr and process.stderr.strip(): sys.stderr.write("    Plot stderr:\n" + process.stderr + "\n") 
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"    Error running modbamtools plot:\n{e}\n    Stdout:\n{e.stdout}\n    Stderr:\n{e.stderr}\n")
    except FileNotFoundError:
        sys.stderr.write("    Error: modbamtools command not found. Is it installed and in your PATH?\n")

def run_modbamtools_calc_meth(phased_bam_path_abs, bed_file_path_abs, analysis_dir, calc_meth_output_filename):
    """
    Calls the external `modbamtools calcMeth` command to calculate methylation statistics.
    The command is run from within the specified analysis_dir.
    Returns the full path to the expected output .txt file, or None if an error occurs.
    """
    sys.stdout.write(f"\n  Running modbamtools calcMeth using BED file {bed_file_path_abs}...\n")
    cmd = ["modbamtools", "calcMeth", "-b", bed_file_path_abs, "-hp", "-o", calc_meth_output_filename, phased_bam_path_abs]
    expected_meth_file_path = os.path.join(analysis_dir, calc_meth_output_filename)
    sys.stdout.write(f"    Executing in '{analysis_dir}': {' '.join(cmd)}\n")
    process_stderr_output = "" 
    try:
        process = subprocess.run(cmd, cwd=analysis_dir, check=True, capture_output=True, text=True)
        process_stderr_output = process.stderr 
        sys.stdout.write("    modbamtools calcMeth completed successfully (exit code 0).\n")
        if process.stdout.strip(): sys.stdout.write("    calcMeth stdout:\n" + process.stdout + "\n")
        if process_stderr_output and process_stderr_output.strip(): sys.stderr.write("    calcMeth stderr:\n" + process_stderr_output + "\n")
        if os.path.exists(expected_meth_file_path): return expected_meth_file_path
        else:
            sys.stderr.write(f"    Error: Expected calcMeth output file '{expected_meth_file_path}' not found (exit code 0).\n")
            if not (process_stderr_output and process_stderr_output.strip()):
                 sys.stderr.write("      (modbamtools calcMeth produced no stderr but failed to create output file)\n")
            try: dir_contents = os.listdir(analysis_dir); sys.stderr.write(f"      Contents of '{analysis_dir}': {dir_contents}\n")
            except Exception as e_ls: sys.stderr.write(f"      Could not list contents of '{analysis_dir}': {e_ls}\n")
            return None
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"    Error running modbamtools calcMeth (non-zero exit code):\n{e}\n    Stdout:\n{e.stdout}\n    Stderr:\n{e.stderr}\n"); return None
    except FileNotFoundError: sys.stderr.write("    Error: modbamtools command not found.\n"); return None

def analyze_calc_meth_output(meth_file_path, min_coverage_diff_meth=5, low_meth_threshold=0.5, epsilon=0.001):
    """
    Parses `modbamtools calcMeth` output, performs Welch's t-test for differential methylation,
    and calculates methylation percent change skewness.
    Returns a list of dictionaries, each representing a row for the summary.
    """
    if not meth_file_path or not os.path.exists(meth_file_path):
        sys.stderr.write(f"calcMeth output file not found: {meth_file_path}\n"); return []
    sys.stdout.write(f"\n  Analyzing differential methylation and skewness from {meth_file_path}...\n")
    
    results_for_summary = []
    expected_column_names = [
        'chr', 'start', 'end', 'name', 
        'metric_all_reads', 'std_all_reads', 'coverage_all_reads',
        'metric_hap1', 'std_hap1', 'coverage_hap1',
        'metric_hap2', 'std_hap2', 'coverage_hap2']
    try:
        df = pd.read_csv(meth_file_path, sep='\t', comment='#', header=None, names=expected_column_names, na_filter=False)
        if df.empty: sys.stdout.write("    calcMeth output file is empty. No analysis performed.\n"); return []
            
        col_widths = { "Region": 30, "H1_Avg_Meth": 15, "H1_StdDev":10, "H1_Coverage": 12, 
                       "H2_Avg_Meth": 15, "H2_StdDev":10, "H2_Coverage": 12, 
                       "T_test_P_Value": 15, "Meth_Skew_(%)": 20}
        header_parts = [
            f"{'Region':<{col_widths['Region']}}",
            f"{'H1_Avg_Meth':>{col_widths['H1_Avg_Meth']}}", f"{'H1_StdDev':>{col_widths['H1_StdDev']}}", f"{'H1_Cov':>{col_widths['H1_Coverage']}}",
            f"{'H2_Avg_Meth':>{col_widths['H2_Avg_Meth']}}", f"{'H2_StdDev':>{col_widths['H2_StdDev']}}", f"{'H2_Cov':>{col_widths['H2_Coverage']}}",
            f"{'T-test_P':>{col_widths['T_test_P_Value']}}", f"{'Meth_Skew_%':>{col_widths['Meth_Skew_(%)']}}"
        ]
        header_str = " ".join(header_parts)
        sys.stdout.write("    " + header_str + "\n"); sys.stdout.write("    " + "-" * len(header_str) + "\n")

        for _, row in df.iterrows():
            try:
                chrom, start, end = row['chr'], int(row['start']), int(row['end'])     
                mean1 = float(row['metric_hap1']) if pd.notna(row['metric_hap1']) else np.nan
                std1 = float(row['std_hap1']) if pd.notna(row['std_hap1']) else np.nan
                nobs1 = int(row['coverage_hap1']) if pd.notna(row['coverage_hap1']) else 0
                mean2 = float(row['metric_hap2']) if pd.notna(row['metric_hap2']) else np.nan
                std2 = float(row['std_hap2']) if pd.notna(row['std_hap2']) else np.nan
                nobs2 = int(row['coverage_hap2']) if pd.notna(row['coverage_hap2']) else 0
            except (KeyError, ValueError) as e_parse: 
                sys.stderr.write(f"    Warning: Could not parse row values: {row}. Error: {e_parse}\n    Skipping.\n"); continue
            
            region_str_out = f"{chrom}:{start+1}-{end}" 
            p_value_ttest = np.nan; p_value_str = "NA"
            methylation_skew_value = np.nan; methylation_skew_display_str = "NA" 

            can_calculate_stats = not (pd.isna(mean1) or pd.isna(mean2) or \
                                   nobs1 < min_coverage_diff_meth or nobs2 < min_coverage_diff_meth)
            
            if can_calculate_stats:
                # Welch's t-test for differential methylation
                if not (pd.isna(std1) or pd.isna(std2) or std1 < 0 or std2 < 0): 
                    try:
                        if nobs1 >= 2 and nobs2 >= 2 : 
                             _, p_value_ttest = ttest_ind_from_stats(mean1, std1, nobs1, mean2, std2, nobs2, equal_var=False)
                             p_value_str = f"{p_value_ttest:.4g}"
                        else: p_value_str = "NA (n<2)"
                    except ZeroDivisionError: p_value_str = "NA (var=0)" 
                    except Exception: p_value_str = f"NA (T-test error)"
                else: p_value_str = "NA (std issue)"
                
                # Methylation Percent Change Skewness
                m1_calc, m2_calc = mean1, mean2 
                if m1_calc < low_meth_threshold and m2_calc < low_meth_threshold:
                    methylation_skew_display_str = "LowMeth"
                    methylation_skew_value = 0.0 
                elif (m1_calc < low_meth_threshold and m2_calc >= low_meth_threshold) or \
                     (m2_calc < low_meth_threshold and m1_calc >= low_meth_threshold):
                    methylation_skew_display_str = "Extreme"
                    if min(m1_calc, m2_calc) < epsilon: 
                        methylation_skew_value = 99999.0 
                    else:
                        methylation_skew_value = ((max(m1_calc, m2_calc) / (min(m1_calc, m2_calc) + epsilon)) - 1) * 100.0
                else: 
                    min_meth = min(m1_calc, m2_calc)
                    max_meth = max(m1_calc, m2_calc)
                    if min_meth + epsilon == 0: 
                        methylation_skew_display_str = "DivByZero"
                        methylation_skew_value = np.nan
                    else:
                        methylation_skew_value = ((max_meth / (min_meth + epsilon)) - 1) * 100.0
                        methylation_skew_display_str = f"{methylation_skew_value:.1f}"
            else:
                p_value_str = "NA (low cov/NaN)"
                methylation_skew_display_str = "NA (low cov/NaN)"
            
            mean1_str = f"{mean1:.2f}" if not pd.isna(mean1) else "NA"
            std1_str = f"{std1:.2f}" if not pd.isna(std1) else "NA"
            mean2_str = f"{mean2:.2f}" if not pd.isna(mean2) else "NA"
            std2_str = f"{std2:.2f}" if not pd.isna(std2) else "NA"

            sys.stdout.write(f"{region_str_out:<{col_widths['Region']}} "
                             f"{mean1_str:>{col_widths['H1_Avg_Meth']}} {std1_str:>{col_widths['H1_StdDev']}} {str(nobs1):>{col_widths['H1_Coverage']}} "
                             f"{mean2_str:>{col_widths['H2_Avg_Meth']}} {std2_str:>{col_widths['H2_StdDev']}} {str(nobs2):>{col_widths['H2_Coverage']}} "
                             f"{p_value_str:>{col_widths['T_test_P_Value']}} {methylation_skew_display_str:>{col_widths['Meth_Skew_(%)']}}\n")
            results_for_summary.append({
                "Region": region_str_out, 
                "H1_Avg_Meth": mean1, "H1_StdDev": std1, "H1_Coverage": nobs1,
                "H2_Avg_Meth": mean2, "H2_StdDev": std2, "H2_Coverage": nobs2, 
                "Diff_Meth_P_Value": p_value_ttest,
                "Methylation_Percent_Change_Skew": methylation_skew_value 
            })
    except pd.errors.EmptyDataError: sys.stderr.write(f"  Warning: pandas.read_csv resulted in empty DataFrame for '{meth_file_path}'.\n")
    except Exception as e: sys.stderr.write(f"  Error processing calcMeth output '{meth_file_path}': {e}\n"); import traceback; traceback.print_exc()
    return results_for_summary

def generate_bed_for_calcmeth(plot_region_str, output_dir):
    """
    Generates a BED file for `modbamtools calcMeth`.
    """
    parsed_region_info = parse_genomic_region(plot_region_str)
    if not parsed_region_info: return None
    chrom, start_0, end_1based_inclusive = parsed_region_info
    bed_file_name = f"region_for_calcMeth_{sanitize_region_for_filename(plot_region_str)}.bed"
    bed_file_path = os.path.join(output_dir, bed_file_name)
    try:
        with open(bed_file_path, 'w') as bf:
            bf.write(f"{chrom}\t{start_0}\t{end_1based_inclusive}\tPlotRegion\n") 
        sys.stdout.write(f"  Generated BED file for calcMeth: {bed_file_path}\n")
        return bed_file_path
    except IOError as e: sys.stderr.write(f"  Error writing BED file {bed_file_path}: {e}\n"); return None

def process_single_bam(input_bam_path, args, parsed_rep_elements_template, num_re_elements, 
                       output_region_chrom_orig, output_region_start_0_pysam, output_region_end_pysam_exclusive):
    """
    Encapsulates the full processing pipeline for a single input BAM file.
    """
    sys.stdout.write(f"\n======================================================================\n")
    sys.stdout.write(f"Processing Input BAM: {input_bam_path}\n")
    sys.stdout.write(f"======================================================================\n")

    input_bam_basename = os.path.splitext(os.path.basename(input_bam_path))[0]
    output_region_sanitized_for_dir = sanitize_region_for_filename(args.output_region) 
    analysis_output_dir = f"{input_bam_basename}_meth_analysis_{output_region_sanitized_for_dir}"
    try:
        os.makedirs(analysis_output_dir, exist_ok=True)
        sys.stdout.write(f"Using output directory for this BAM: {os.path.abspath(analysis_output_dir)}\n")
    except OSError as e:
        sys.stderr.write(f"Error creating output directory {analysis_output_dir} for {input_bam_path}: {e}\n")
        return None 

    auto_output_bam_name = f"{input_bam_basename}_phased_in_{output_region_sanitized_for_dir}.bam"
    full_output_bam_path = os.path.join(analysis_output_dir, auto_output_bam_name)
    sys.stdout.write(f"Phased output BAM will be saved to: {full_output_bam_path}\n")

    all_reads_in_output_region_data = [] 
    infile = None; outfile = None; phased_bam_created_successfully = False; reads_written_count_final = 0
    phasing_summary_data = {
        "input_bam": input_bam_path, "reads_in_output_region": 0, "reads_for_phasing_definition": 0,
        "H1_reads_phased": 0, "H2_reads_phased": 0, "phasing_combined_p_value": np.nan, 
        "phasing_stats_per_re": [], 
        "differential_methylation_results": [] }

    try:
        infile = pysam.AlignmentFile(input_bam_path, "rb"); bam_references = infile.references
        def reconcile_chrom_name(orig_chrom_name, bam_refs, region_label=""):
            if orig_chrom_name in bam_refs: return orig_chrom_name
            prefix_chr = "chr" + orig_chrom_name;
            if prefix_chr in bam_refs and not orig_chrom_name.startswith("chr"): return prefix_chr
            strip_chr = orig_chrom_name[3:] if orig_chrom_name.startswith("chr") else None
            if strip_chr and strip_chr in bam_refs: return strip_chr
            for ref_in_bam in bam_refs:
                if ref_in_bam.lower() == orig_chrom_name.lower():
                    sys.stdout.write(f"INFO: Matched {region_label} chr '{orig_chrom_name}' to BAM ref '{ref_in_bam}'.\n")
                    return ref_in_bam
            return None

        output_region_chrom = reconcile_chrom_name(output_region_chrom_orig, bam_references, "Output region")
        if not output_region_chrom:
            sys.stderr.write(f"Error: Output region chr '{output_region_chrom_orig}' not found in BAM '{input_bam_path}'.\n"); return phasing_summary_data
        
        current_parsed_re_elements = []
        for re_info_template in parsed_rep_elements_template: 
            re_info = re_info_template.copy()
            re_info["chrom"] = reconcile_chrom_name(re_info["orig_chrom"], bam_references, f"RE '{re_info['orig_region_str']}'")
            if not re_info["chrom"]:
                sys.stderr.write(f"Error: RE chr '{re_info['orig_chrom']}' not found for BAM '{input_bam_path}'.\n"); return phasing_summary_data
            current_parsed_re_elements.append(re_info)
        
        outfile = pysam.AlignmentFile(full_output_bam_path, "wb", header=infile.header) 
        sys.stdout.write(f"\nProcessing reads in primary output region: {output_region_chrom}:{output_region_start_0_pysam+1}-{output_region_end_pysam_exclusive}...\n")
        for read in infile.fetch(output_region_chrom, output_region_start_0_pysam, output_region_end_pysam_exclusive): 
            if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
            counts_vector = [np.nan] * num_re_elements
            for re_idx, re_info_dyn in enumerate(current_parsed_re_elements): 
                counts_vector[re_idx] = count_repeats_in_read_segment(
                    read, re_info_dyn['chrom'], re_info_dyn['start_0'], re_info_dyn['end_0_exclusive'], 
                    re_info_dyn['motif'], args.min_overlap_bp_in_rep)
            all_reads_in_output_region_data.append({'read_obj': read, 'qname': read.query_name, 'counts_vector': counts_vector})
        
        phasing_summary_data["reads_in_output_region"] = len(all_reads_in_output_region_data)
        sys.stdout.write(f"  Collected {phasing_summary_data['reads_in_output_region']} primary reads from the output region.\n")
        
        if not all_reads_in_output_region_data: sys.stdout.write("  No reads found. Output BAM will be empty.\n")
        
        reads_for_clustering = [item for item in all_reads_in_output_region_data if not all(np.isnan(c) for c in item['counts_vector'])]
        phasing_summary_data["reads_for_phasing_definition"] = len(reads_for_clustering)
        sys.stdout.write(f"  {phasing_summary_data['reads_for_phasing_definition']} reads have >=1 valid RE count and will be used for defining clusters.\n")

        qname_to_hap_int = {}; clustering_ok = False; h1_label_idx = -1 
        if len(reads_for_clustering) < args.min_reads_for_clustering:
            sys.stdout.write(f"  Not enough reads ({len(reads_for_clustering)}) for clustering (min {args.min_reads_for_clustering}). Reads will be tagged HP:i:0.\n")
        else:
            feature_matrix = np.array([item['counts_vector'] for item in reads_for_clustering])
            col_means_for_imputation = np.nanmean(feature_matrix, axis=0)
            for i in range(num_re_elements): 
                if np.isnan(col_means_for_imputation[i]): 
                    col_means_for_imputation[i] = 0 
                    sys.stderr.write(f"  Warning: RE {i+1} ('{current_parsed_re_elements[i]['motif']}') had no valid counts. Imputing with 0.\n") 
            nan_indices = np.where(np.isnan(feature_matrix)); feature_matrix[nan_indices] = np.take(col_means_for_imputation, nan_indices[1])
            has_variation = any(len(np.unique(feature_matrix[:, i])) > 1 for i in range(num_re_elements) if feature_matrix.shape[1] > i)
            if not has_variation and num_re_elements > 0:
                 sys.stdout.write("  Not enough variation in repeat counts. Reads will be tagged HP:i:0.\n")
            else:
                sys.stdout.write("  Attempting K-Means clustering (K=2)...\n")
                try:
                    kmeans = KMeans(n_clusters=2, random_state=42, n_init='auto').fit(feature_matrix)
                    cluster_labels_for_subset = kmeans.labels_
                    means_of_first_re_for_labeling = feature_matrix[:, 0] 
                    cluster_0_mean_re0 = np.mean(means_of_first_re_for_labeling[cluster_labels_for_subset == 0]) if len(means_of_first_re_for_labeling[cluster_labels_for_subset == 0]) > 0 else np.nan
                    cluster_1_mean_re0 = np.mean(means_of_first_re_for_labeling[cluster_labels_for_subset == 1]) if len(means_of_first_re_for_labeling[cluster_labels_for_subset == 1]) > 0 else np.nan
                    if np.isnan(cluster_0_mean_re0) and np.isnan(cluster_1_mean_re0): sys.stderr.write("    Warning: Both cluster means for RE0 NaN.\n")
                    elif np.isnan(cluster_0_mean_re0): h1_label_idx = 1 ; sys.stderr.write("    Warning: Cluster 0 mean for RE0 NaN. H1 -> Cluster 1.\n")
                    elif np.isnan(cluster_1_mean_re0): h1_label_idx = 0 ; sys.stderr.write("    Warning: Cluster 1 mean for RE0 NaN. H1 -> Cluster 0.\n")
                    elif cluster_0_mean_re0 <= cluster_1_mean_re0: h1_label_idx = 0
                    else: h1_label_idx = 1
                    if h1_label_idx != -1: 
                        clustering_ok = True; h2_label_idx = 1 - h1_label_idx
                        sys.stdout.write("\n--- Phasing Results for this BAM ---\n"); individual_p_values = []
                        for idx, item in enumerate(reads_for_clustering):
                            if cluster_labels_for_subset[idx] == h1_label_idx: qname_to_hap_int[item['qname']] = 1 
                            elif cluster_labels_for_subset[idx] == h2_label_idx: qname_to_hap_int[item['qname']] = 2
                        phasing_summary_data["H1_reads_phased"] = list(qname_to_hap_int.values()).count(1)
                        phasing_summary_data["H2_reads_phased"] = list(qname_to_hap_int.values()).count(2)
                        re_stats_list_for_summary = []
                        for re_idx in range(num_re_elements):
                            counts_h1_orig = [item['counts_vector'][re_idx] for item in reads_for_clustering if qname_to_hap_int.get(item['qname']) == 1]
                            counts_h2_orig = [item['counts_vector'][re_idx] for item in reads_for_clustering if qname_to_hap_int.get(item['qname']) == 2]
                            valid_h1 = np.array([c for c in counts_h1_orig if not np.isnan(c)]); valid_h2 = np.array([c for c in counts_h2_orig if not np.isnan(c)])
                            p_value = np.nan
                            if len(valid_h1)>=1 and len(valid_h2)>=1:
                                try: 
                                    if not (len(np.unique(valid_h1))==1 and len(np.unique(valid_h2))==1 and np.unique(valid_h1)[0]==np.unique(valid_h2)[0]):
                                        _, p_value = mannwhitneyu(valid_h1, valid_h2, alternative='two-sided', nan_policy='propagate')
                                        if not np.isnan(p_value) and 0<p_value<1: individual_p_values.append(p_value)
                                except ValueError: pass 
                            re_stats_list_for_summary.append({
                                "motif": current_parsed_re_elements[re_idx]['motif'],
                                "h1_stats_str": f"{np.mean(valid_h1):.1f} ({int(np.min(valid_h1))}-{int(np.max(valid_h1))})" if len(valid_h1)>0 else "NA",
                                "h2_stats_str": f"{np.mean(valid_h2):.1f} ({int(np.min(valid_h2))}-{int(np.max(valid_h2))})" if len(valid_h2)>0 else "NA",
                                "mwu_p_value": p_value
                            })
                        phasing_summary_data['phasing_stats_per_re'] = re_stats_list_for_summary
                        if len(individual_p_values) > 0:
                            valid_p_fish = [p for p in individual_p_values if 0 < p < 1] 
                            if len(valid_p_fish) == 1: phasing_summary_data["phasing_combined_p_value"] = valid_p_fish[0]
                            elif len(valid_p_fish) > 1:
                                fish_s = -2 * np.sum(np.log(valid_p_fish)); df_f = 2 * len(valid_p_fish)
                                phasing_summary_data["phasing_combined_p_value"] = chi2.sf(fish_s, df_f)
                    else: clustering_ok = False; sys.stdout.write("  Phasing failed (H1/H2 labels).\n")
                except Exception as e: sys.stderr.write(f"  Clustering error: {e}\n"); import traceback; traceback.print_exc(); clustering_ok = False
        
        sys.stdout.write("\nWriting output BAM file...\n"); reads_written_count_final = 0; final_counts = {1:0, 2:0, 0:0} 
        for item in all_reads_in_output_region_data: 
            read = item['read_obj']; qname = item['qname']
            hap_int = qname_to_hap_int.get(qname, 0) 
            try: read.set_tag(args.hp_tag, hap_int, value_type='i'); outfile.write(read); reads_written_count_final += 1; final_counts[hap_int] +=1
            except Exception as e: sys.stderr.write(f"  Warning: Tag/write failed for {qname}: {e}\n")
        sys.stdout.write(f"Finished phasing step. {reads_written_count_final} reads written to {full_output_bam_path}.\n") 
        for hap_val, count_val in final_counts.items(): sys.stdout.write(f"  HP:i:{hap_val} reads written: {count_val}\n")
        if reads_written_count_final > 0: phased_bam_created_successfully = True 
    except FileNotFoundError: sys.stderr.write(f"Error: Input BAM file '{input_bam_path}' not found.\n"); return phasing_summary_data 
    except ImportError as e: sys.stderr.write(f"Error: Missing package: {e}.\n"); return phasing_summary_data
    except Exception as e: sys.stderr.write(f"Error during phasing for {input_bam_path}: {e}\n"); import traceback; traceback.print_exc(); return phasing_summary_data
    finally:
        if infile: infile.close()
        if outfile: outfile.close() 

    if phased_bam_created_successfully: 
        abs_phased_bam_path = os.path.abspath(full_output_bam_path) 
        plot_basename_prefix = f"{input_bam_basename}_plot_{sanitize_region_for_filename(args.plot_region)}" 
        calc_meth_filename = f"{input_bam_basename}_calcMeth_{sanitize_region_for_filename(args.plot_region)}.txt" 
        auto_generated_bed_path = generate_bed_for_calcmeth(args.plot_region, analysis_output_dir)
        if auto_generated_bed_path:
            auto_generated_bed_path_abs = os.path.abspath(auto_generated_bed_path)
            try:
                sys.stdout.write(f"\nIndexing output BAM for modbamtools: {full_output_bam_path}...\n") 
                pysam.index(full_output_bam_path); sys.stdout.write("  Indexing complete.\n")
            except Exception as e: sys.stderr.write(f"  Error indexing output BAM: {e}.\n")
            run_modbamtools_plot(abs_phased_bam_path, args.plot_region, analysis_output_dir, plot_basename_prefix)
            meth_file_path = run_modbamtools_calc_meth(abs_phased_bam_path, auto_generated_bed_path_abs, analysis_output_dir, calc_meth_filename) 
            if meth_file_path and os.path.exists(meth_file_path): 
                phasing_summary_data["differential_methylation_results"] = analyze_calc_meth_output(meth_file_path, args.min_coverage_diff_meth)
        else: sys.stderr.write("Skipping modbamtools calcMeth for {input_bam_path} as BED generation failed.\n")
    else: sys.stdout.write(f"\nSkipping methylation analysis for {input_bam_path}: Phased BAM not created/empty.\n")
    return phasing_summary_data

def main():
    parser = argparse.ArgumentParser(
        description="Version 3.0: Phase reads using repeat copy numbers, tag BAM, run modbamtools, analyze differential methylation, and summarize."
    )
    parser.add_argument("-i", "--input_bams", required=True, nargs='+', help="Path(s) to input BAM file(s) (must be indexed).")
    parser.add_argument("--output_region", required=True, 
                        help="Primary genomic region for output BAM content (e.g., 'chr1:1000-2000').")
    parser.add_argument("--rep_element_regions", required=True, nargs='+',
                        help="Space-separated list of repetitive element genomic regions.")
    parser.add_argument("--motifs", required=True, nargs='+', 
                        help="Space-separated list of repetitive motif sequences.")
    parser.add_argument("--min_overlap_bp_in_rep", type=int, default=10, 
                        help="Min bp overlap in RE region for valid count (default: 10).")
    parser.add_argument("--min_reads_for_clustering", type=int, default=5, 
                        help="Min reads with valid RE count for clustering (default: 5).")
    parser.add_argument("--hp_tag", default="HP", help="BAM tag for haplotype (default: HP).")
    parser.add_argument("--plot_region", type=str, 
                        help="Region for modbamtools plot and calcMeth BED. Defaults to --output_region.")
    parser.add_argument("--min_coverage_diff_meth", type=int, default=5,
                        help="Min coverage for H1/H2 in BED region for diff meth test (default: 5).")
    args = parser.parse_args()

    if not args.plot_region: 
        args.plot_region = args.output_region 
        sys.stdout.write(f"INFO: --plot_region not set, defaulting to --output_region for methylation analysis: {args.plot_region}\n")
    if len(args.rep_element_regions) != len(args.motifs):
        sys.stderr.write("Error: Number of --rep_element_regions must match number of --motifs.\n"); sys.exit(1)
    if not all(args.motifs):
        sys.stderr.write("Error: All motifs must be non-empty strings.\n"); sys.exit(1)

    output_region_parsed_tuple = parse_genomic_region(args.output_region) 
    if not output_region_parsed_tuple: sys.exit(1)
    output_region_chrom_orig, output_region_start_0_pysam, output_region_end_pysam_exclusive = output_region_parsed_tuple[0], output_region_parsed_tuple[1], output_region_parsed_tuple[2]

    num_re_elements = len(args.rep_element_regions)
    parsed_rep_elements_template = [] 
    for i in range(num_re_elements):
        parsed_re_tuple = parse_genomic_region(args.rep_element_regions[i])
        if not parsed_re_tuple: sys.exit(1)
        parsed_rep_elements_template.append({
            "orig_region_str": args.rep_element_regions[i], "orig_chrom": parsed_re_tuple[0],
            "chrom": parsed_re_tuple[0], "start_0": parsed_re_tuple[1],
            "end_0_exclusive": parsed_re_tuple[2], "motif": args.motifs[i] 
        })
    
    all_bam_summaries = []
    script_run_dir = os.getcwd() 

    for input_bam_file_path in args.input_bams:
        summary_for_this_bam = process_single_bam(
            input_bam_file_path, args, parsed_rep_elements_template, num_re_elements,
            output_region_chrom_orig, output_region_start_0_pysam, output_region_end_pysam_exclusive
        )
        if summary_for_this_bam:
            all_bam_summaries.append(summary_for_this_bam)
        else:
            sys.stderr.write(f"Warning: Processing returned no summary for {input_bam_file_path}.\n")

    if all_bam_summaries:
        generate_master_summary_report(all_bam_summaries, script_run_dir)
    else:
        sys.stdout.write("\nNo BAM files were successfully processed to generate a master summary.\n")

def generate_master_summary_report(all_summaries, base_output_dir):
    """
    Generates and prints a master summary report for all processed BAM files.
    Writes the summary to a timestamped TSV file.
    """
    sys.stdout.write("\n\n====================================================================================================================================================\n") 
    sys.stdout.write("                                                                    Master Summary Report\n")
    sys.stdout.write("====================================================================================================================================================\n") 

    report_lines = []
    static_header = ["Input_BAM", "Reads_in_Region", "Reads_for_Phasing", "H1_Phased_Reads", "H2_Phased_Reads", "Phasing_P_Value"]
    dynamic_header = []
    if all_summaries and all_summaries[0].get('phasing_stats_per_re'):
        for i, re_stat in enumerate(all_summaries[0]['phasing_stats_per_re']):
            re_label = f"RE{i+1}_{re_stat['motif']}"
            dynamic_header.extend([f"{re_label}_H1_Counts", f"{re_label}_H2_Counts", f"{re_label}_MWU_P"])
    meth_header = [ "BED_Region", "H1_Avg_Meth(%)", "H1_StdDev", "H1_Coverage", "H2_Avg_Meth(%)", "H2_StdDev", "H2_Coverage", 
                    "Diff_Meth_P(T-test)", "Meth_Pct_Change_Skew(%)"] 
    full_header = static_header + dynamic_header + meth_header
    report_lines.append("\t".join(full_header))
    
    col_widths = {"Input_BAM": 30, "Reads_in_Region":10, "Reads_for_Phasing":10, "H1_Phased_Reads":10, "H2_Phased_Reads":10, "Phasing_P_Value": 12, 
                  "BED_Region": 25, "H1_Avg_Meth(%)": 15, "H1_StdDev":10, "H1_Coverage": 12, 
                  "H2_Avg_Meth(%)": 15, "H2_StdDev":10, "H2_Coverage": 12, 
                  "Diff_Meth_P(T-test)": 18, "Meth_Pct_Change_Skew(%)":25} 
    
    console_header_parts = [f"{'Input_BAM':<{col_widths['Input_BAM']}}",
                            f"{'Total_Reads':>{col_widths['Reads_in_Region']}}", 
                            f"{'Phasing_Pool':>{col_widths['Reads_for_Phasing']}}",
                            f"{'H1_Count':>{col_widths['H1_Phased_Reads']}}", 
                            f"{'H2_Count':>{col_widths['H2_Phased_Reads']}}", 
                            f"{'Phasing_Pval':>{col_widths['Phasing_P_Value']}}"]
    
    if dynamic_header:
        for i in range(len(all_summaries[0]['phasing_stats_per_re'])):
            re_label = f"RE{i+1}"
            col_widths[f"RE{i+1}_H1"] = 18
            col_widths[f"RE{i+1}_H2"] = 18
            col_widths[f"RE{i+1}_P"] = 12
            console_header_parts.extend([
                f"{re_label}_H1_Counts".center(col_widths[f'RE{i+1}_H1']),
                f"{re_label}_H2_Counts".center(col_widths[f'RE{i+1}_H2']),
                f"{re_label}_Pval".center(col_widths[f'RE{i+1}_P'])
            ])

    console_header_parts.extend([f"{'BED_Region':<{col_widths['BED_Region']}}",
                                 f"{'H1_Meth%':>{col_widths['H1_Avg_Meth(%)']}}",
                                 f"{'H1_StdD':>{col_widths['H1_StdDev']}}",
                                 f"{'H1_Cov':>{col_widths['H1_Coverage']}}",
                                 f"{'H2_Meth%':>{col_widths['H2_Avg_Meth(%)']}}",
                                 f"{'H2_StdD':>{col_widths['H2_StdDev']}}",
                                 f"{'H2_Cov':>{col_widths['H2_Coverage']}}",
                                 f"{'Diff_Meth_Pval':>{col_widths['Diff_Meth_P(T-test)']}}",
                                 f"{'Meth_Pct_Change%':>{col_widths['Meth_Pct_Change_Skew(%)']}}"]) 

    console_header_str = " ".join(console_header_parts)
    sys.stdout.write(console_header_str + "\n")
    sys.stdout.write("-" * len(console_header_str) + "\n")

    for summary in all_summaries:
        if not summary: continue 
        bam_name_short = os.path.basename(summary["input_bam"])
        reads_out_reg = str(summary.get("reads_in_output_region", "NA"))
        reads_phasing = str(summary.get("reads_for_phasing_definition", "NA"))
        h1_phased_count = str(summary.get("H1_reads_phased", "NA"))
        h2_phased_count = str(summary.get("H2_reads_phased", "NA"))
        phasing_p = f"{summary['phasing_combined_p_value']:.3g}" if pd.notna(summary['phasing_combined_p_value']) else "NA"
        
        re_stats = summary.get('phasing_stats_per_re', [])
        
        diff_meth_results = summary.get("differential_methylation_results", [])
        if diff_meth_results:
            for i, meth_res in enumerate(diff_meth_results):
                row_parts_console = []
                row_parts_file = []
                
                if i == 0:
                    row_parts_console.extend([f"{bam_name_short:<{col_widths['Input_BAM']}}",
                                              f"{reads_out_reg:>{col_widths['Reads_in_Region']}}",
                                              f"{reads_phasing:>{col_widths['Reads_for_Phasing']}}",
                                              f"{h1_phased_count:>{col_widths['H1_Phased_Reads']}}",
                                              f"{h2_phased_count:>{col_widths['H2_Phased_Reads']}}",
                                              f"{phasing_p:>{col_widths['Phasing_P_Value']}}"])
                    row_parts_file.extend([bam_name_short, reads_out_reg, reads_phasing, h1_phased_count, h2_phased_count, phasing_p])
                    for re_idx, re_stat in enumerate(re_stats):
                        mwu_p_val_str = f"{re_stat.get('mwu_p_value', np.nan):.3g}" if pd.notna(re_stat.get('mwu_p_value')) else "NA"
                        row_parts_console.extend([
                            f"{re_stat.get('h1_stats_str', 'NA'):>{col_widths[f'RE{re_idx+1}_H1']}}",
                            f"{re_stat.get('h2_stats_str', 'NA'):>{col_widths[f'RE{re_idx+1}_H2']}}",
                            f"{mwu_p_val_str:>{col_widths[f'RE{re_idx+1}_P']}}"
                        ])
                        row_parts_file.extend([re_stat.get('h1_stats_str', 'NA'), re_stat.get('h2_stats_str', 'NA'), mwu_p_val_str])
                else: 
                    row_parts_console.append(f"{'':<{col_widths['Input_BAM']}}") 
                    row_parts_console.extend([''] * (5 + 3*len(re_stats))) 
                    row_parts_file.extend([bam_name_short, reads_out_reg, reads_phasing, h1_phased_count, h2_phased_count, phasing_p])
                    row_parts_file.extend([''] * (3*len(re_stats))) 
                

                region = meth_res["Region"]
                h1_meth = f"{meth_res['H1_Avg_Meth']:.2f}" if pd.notna(meth_res['H1_Avg_Meth']) else "NA"
                h1_std = f"{meth_res['H1_StdDev']:.2f}" if pd.notna(meth_res['H1_StdDev']) else "NA"
                h1_cov = str(meth_res['H1_Coverage'])
                h2_meth = f"{meth_res['H2_Avg_Meth']:.2f}" if pd.notna(meth_res['H2_Avg_Meth']) else "NA"
                h2_std = f"{meth_res['H2_StdDev']:.2f}" if pd.notna(meth_res['H2_StdDev']) else "NA"
                h2_cov = str(meth_res['H2_Coverage'])
                diff_p = f"{meth_res['Diff_Meth_P_Value']:.3g}" if pd.notna(meth_res['Diff_Meth_P_Value']) else "NA"
                
                meth_skew_val = meth_res['Methylation_Percent_Change_Skew']
                meth_skew_display_str = "NA"
                if pd.notna(meth_skew_val):
                    if meth_skew_val == 99999.0: meth_skew_display_str = "Extreme"
                    elif meth_skew_val == 0.0: meth_skew_display_str = "LowMeth/Equal"
                    else: meth_skew_display_str = f"{meth_skew_val:.1f}"
                
                row_parts_console.extend([f"{region:<{col_widths['BED_Region']}}",
                                          f"{h1_meth:>{col_widths['H1_Avg_Meth(%)']}}", f"{h1_std:>{col_widths['H1_StdDev']}}", f"{h1_cov:>{col_widths['H1_Coverage']}}",
                                          f"{h2_meth:>{col_widths['H2_Avg_Meth(%)']}}", f"{h2_std:>{col_widths['H2_StdDev']}}", f"{h2_cov:>{col_widths['H2_Coverage']}}",
                                          f"{diff_p:>{col_widths['Diff_Meth_P(T-test)']}}", f"{meth_skew_display_str:>{col_widths['Meth_Pct_Change_Skew(%)']}}"])
                row_parts_file.extend([region, h1_meth, h1_std, h1_cov, h2_meth, h2_std, h2_cov, diff_p, meth_skew_display_str])
                
                report_lines.append("\t".join(row_parts_file))
                sys.stdout.write(" ".join(row_parts_console) + "\n")
        else: 
            row_parts_console = [f"{bam_name_short:<{col_widths['Input_BAM']}}",
                                 f"{reads_out_reg:>{col_widths['Reads_in_Region']}}", f"{reads_phasing:>{col_widths['Reads_for_Phasing']}}",
                                 f"{h1_phased_count:>{col_widths['H1_Phased_Reads']}}", f"{h2_phased_count:>{col_widths['H2_Phased_Reads']}}",
                                 f"{phasing_p:>{col_widths['Phasing_P_Value']}}"]
            row_parts_file = [bam_name_short, reads_out_reg, reads_phasing, h1_phased_count, h2_phased_count, phasing_p]
            
            for re_idx, re_stat in enumerate(re_stats):
                mwu_p_val_str = f"{re_stat.get('mwu_p_value', np.nan):.3g}" if pd.notna(re_stat.get('mwu_p_value')) else "NA"
                row_parts_console.extend([
                    f"{re_stat.get('h1_stats_str', 'NA'):>{col_widths[f'RE{re_idx+1}_H1']}}",
                    f"{re_stat.get('h2_stats_str', 'NA'):>{col_widths[f'RE{re_idx+1}_H2']}}",
                    f"{mwu_p_val_str:>{col_widths[f'RE{re_idx+1}_P']}}"
                ])
                row_parts_file.extend([re_stat.get('h1_stats_str', 'NA'), re_stat.get('h2_stats_str', 'NA'), mwu_p_val_str])

            num_meth_cols = 9 
            row_parts_console.extend([f"{'NA':<{col_widths['BED_Region']}}"] + [f"{'NA':>{list(col_widths.values())[i]}}" for i in range(-num_meth_cols,0)])
            row_parts_file.extend(['NA'] * num_meth_cols)
            
            report_lines.append("\t".join(row_parts_file))
            sys.stdout.write(" ".join(row_parts_console) + "\n")

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_filename = f"master_phasing_meth_summary_{timestamp}.tsv"
    summary_filepath = os.path.join(base_output_dir, summary_filename) 
    try:
        with open(summary_filepath, 'w') as f_summary:
            for line in report_lines:
                f_summary.write(line + "\n")
        sys.stdout.write(f"\nMaster summary report written to: {summary_filepath}\n")
    except IOError as e:
        sys.stderr.write(f"Error writing master summary report to {summary_filepath}: {e}\n")

if __name__ == "__main__":
    main()
