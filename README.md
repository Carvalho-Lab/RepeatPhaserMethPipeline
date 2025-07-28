# **RepeatPhaserMethPipeline**

A Python pipeline to phase long-read sequencing data using tandem repeat polymorphisms, tag reads by haplotype, and perform downstream methylation analysis to assess allelic imbalance and X-inactivation.

## **Description**

RepeatPhaserMethPipeline.py is a command-line tool designed for haplotype-specific methylation analysis, particularly for studying X-chromosome inactivation. It leverages informative tandem repeat polymorphisms present on long reads (e.g., from PacBio HiFi) to phase them into their respective haplotypes.

The script performs the following steps for one or more input BAM files:

1. **Collects Reads:** Identifies all primary long reads that overlap a user-defined genomic region.
2. **Quantifies Repeats:** For each read, it counts the copy number of one or more user-specified tandem repeat motifs within corresponding genomic loci.
3. **Phases Haplotypes:** Uses K-Means clustering (K=2) on the repeat count data to automatically partition reads into two haplotype groups (H1 and H2).
4. **Validates Phasing:** Calculates p-values using the Mann-Whitney U test for each repeat marker and combines them using Fisher's method to provide an overall confidence score for the phasing.
5. **Generates Phased BAM:** Creates a new, indexed BAM file for each input sample, where every read from the target region is tagged with its haplotype (HP\:i:1 for H1, HP\:i:2 for H2, or HP\:i:0 for unphased).  This allows for easy viewing in IGV or other programs.
6. **Automates Methylation Analysis:** Automatically runs modbamtools plot and modbamtools calcMeth on the newly phased BAM file.
7. **Calculates Differential Methylation & Skewness:** Performs a Welch's t-test to determine if methylation is significantly different between the two haplotypes and calculates a methylation skewness metric.
8. **Generates Summary Report:** Produces a consolidated TSV report summarizing the phasing and methylation statistics for all processed BAM files.

## **Dependencies**

This script requires Python 3 and several third-party libraries.

1. **Required Python Libraries:**

* pysam
* numpy
* scikit-learn
* scipy
* pandas

You can install these using pip\:pip install pysam numpy scikit-learn scipy pandas



1. **External Tools:**

* **modbamtools**: This tool must be installed and available in your system's PATH. Please follow the installation instructions on the [<u>modbamtools GitHub page</u>](https://www.google.com/search?q=https://github.com/r937/modbamtools).
* **samtools**: Required by pysam for indexing BAM files. Ensure it is installed and in your PATH.

## **Usage**

```bash
python3 repeat_phaser_meth_pipeline.py -i <input1.bam> [<input2.bam> ...] \
                                       --output_region <chr:start-end> \
                                       --rep_element_regions <chr:start-end> [<chr:start-end> ...] \
                                       --motifs <motif1> [<motif2> ...] \
                                       [OPTIONS]
```


### **Required Arguments:**

* `-i, --input_bams`: One or more paths to the input BAM files (must be indexed).
* `--output_region`: The primary genomic region of interest (e.g., chrX:10000-20000). All primary reads overlapping this region will be processed and included in the output BAM.
* `--rep_element_regions`: A space-separated list of genomic regions containing the tandem repeats to be used for phasing (e.g., chrX:12000-12500 chrX:18000-18300).
* `--motifs`: A space-separated list of repeat motifs corresponding to each --rep\_element\_regions (e.g., CAG GGC). The number of motifs must match the number of regions.

### **Optional Arguments:**

* `--plot_region`: The specific sub-region to use for generating plots and methylation statistics with modbamtools (e.g., chrX:15000-16000). If not provided, this defaults to the --output\_region.
* `--min_reads_for_clustering`: Minimum number of reads with at least one valid repeat count required to attempt clustering. **Default: 5**.
* `--min_overlap_bp_in_rep`: Minimum number of base pairs a read must align within a repetitive element region for its repeat count to be considered valid. **Default: 10**.
* `--min_coverage_diff_meth`: Minimum read coverage required for both H1 and H2 in a BED region to perform the differential methylation t-test. **Default: 5**.
* `--hp_tag`: The BAM tag to use for storing haplotype information. **Default: HP**.

## **Example Command**

This command processes two BAM files, phasing them based on a CAG repeat and a GGC repeat on chromosome X. It then generates plots and methylation statistics for a specific sub-region.

```bash
python3 repeat_phaser_meth_pipeline.py \
    --input_bams sampleA.bam sampleB.bam \
    --output_region chrX:67543951-67546579 \
    --rep_element_regions chrX:67545314-67545420 chrX:67546511-67546568 \
    --motifs GCA GGC \
    --plot_region chrX:67544000-67546000`
```



## **Output**

The script generates the following outputs for each input BAM file:

1. **Analysis Directory:** A new directory is created for each input BAM, named like \[input\_bam\_basename]\_meth\_analysis\_\[output\_region].
2. **Phased BAM File:** Inside the analysis directory, a new BAM file (\[input\_bam\_basename]\_phased\_in\_\[output\_region].bam) is created. It contains only reads from the --output\_region, tagged with HP\:i:1, HP\:i:2, or HP\:i:0.
3. **modbamtools Plot:** An HTML plot file generated by modbamtools plot is saved in the analysis directory.
4. **modbamtools Methylation Stats:** A .txt file with methylation statistics from modbamtools calcMeth.
5. **Auto-generated BED File:** A BED file used as input for modbamtools calcMeth.

### **Master Summary Report**

After all BAM files are processed, a single, timestamped TSV file (e.g., master\_phasing\_meth\_summary\_YYYYMMDD\_HHMMSS.tsv) is created in the directory where the script was run. This file contains a consolidated summary of the key phasing and methylation statistics for all samples, making it easy to compare results.
