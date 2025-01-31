# Tutorial: Metagenomic Analysis of Bacterial Genomes Using Illumina Short-Read Data

### Learning Objectives
This tutorial provides a step-by-step command-line workflow for analyzing bacterial genomes from Illumina short-read metagenomic data. By following these steps, you’ll be able to:
1. Conduct quality control on Illumina reads.
2. Assemble short reads into contigs using **metaSPAdes**.
3. Bin contigs into genomes (Metagenome-Assembled Genomes or MAGs).
4. Taxonomically classify and annotate these genomes.
5. Assess antimicrobial resistance (AMR) potential and visualize taxonomic data.

### Prerequisites
To complete this tutorial, you need:
- **Basic Linux command-line** knowledge.
- **Installed tools**: FastQC, Fastp, seqtk, Bbmap, metaSPAdes, MetaBAT2, Kraken2, CheckM, Prokka, Abricate (for AMR prediction), and Pavian or Krona for visualization.

### Dataset

Raw files were grabbed from Guo, X., Tang, N., Lei, H., Fang, Q., Liu, L., Zhou, Q., & Song, C. (2021). Metagenomic Analysis of Antibiotic Resistance Genes in Untreated Wastewater From Three Different Hospitals. Frontiers in microbiology, 12, 709051. https://doi.org/10.3389/fmicb.2021.709051. Selected ones were spiked with CPE plasmids.

  ```
  # Raw
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR142/072/SRR14297772/SRR14297772_1.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR142/072/SRR14297772/SRR14297772_2.fastq.gz
  # Spiked: c/o Aarthi and Arun
  # Filenames: SRR14297772_cpe107_1.fastq.gz and SRR14297772_cpe107_2.fastq.gz

  ```


## Step 1: Quality Control of Illumina Reads

 
Quality control ensures your Illumina reads are suitable for assembly. **FastQC** will identify quality issues, and **Fastp** will trim low-quality bases and adapters.

 

1. **Run FastQC to Assess Read Quality**:
   ``` 
   fastqc raw_reads/*.fastq.gz -o fastqc_output/
   ```

   - **What It Does**: FastQC generates quality metrics and HTML reports, including GC content, read length distribution, and base quality scores.
   - **Output**: HTML reports in the `fastqc_output` directory. Open these to check for any quality issues.

2. **Clean Reads Using Fastp and Hostile**:
   ``` 
   fastp -i raw_reads/sample_R1.fastq -I raw_reads/sample_R2.fastq    -o trimmed_sample_R1.fastq -O trimmed_sample_R2.fastq    -h fastp_report.html -j fastp_report.json --length_required 50
   ```
   
   ```
   # We will use the spiked pair
   fastp -i SRR14297772_cpe107_1.fastq.gz -I SRR14297772_cpe107_2.fastq.gz    -o trimmed_SRR14297772_cpe107_1.fastq.gz -O trimmed_SRR14297772_cpe107_2.fastq.gz -h fastp_report.html -j fastp_report.json --length_required 50 
   ```

   - **What It Does**: Fastp removes adapters, trims low-quality bases, and discards reads shorter than 50 bp. Hostile gets rid of host (mainly) human reads.
   - **Key Options for Fastp**:
     - `-h` and `-j`: Generate HTML and JSON reports with trimming metrics.
     - `--length_required`: Discards reads shorter than 50 bp.
   - **Hostile**: https://github.com/bede/hostile
      ```
      # Create and activate a conda env 
      conda create -y -n hostile -c conda-forge -c bioconda hostile
      conda activate hostile
      conda activate --stack amr # to access packages from the amr env
      ```
      ```
      # Run Hostile on paired short reads
      hostile clean --fastq1 SRR14297772_cpe107_1.fastq.gz --fastq2 SRR14297772_cpe107_2.fastq.gz -o - > SRR14297772_cpe107.interleaved.fastq
      ```
      ```
      # Bin interleaved fastq files into clean.fastq1 and clean.fastq2 using seqtk
      seqtk seq -1 SRR14297772_cpe107.interleaved.fastq > clean.SRR14297772_cpe107_1.fastq
      seqtk seq -2 SRR14297772_cpe107.interleaved.fastq > clean.SRR14297772_cpe107_2.fastq
      ```
      ```
      #Compress all fastq files
      gzip SRR14297772_cpe107.interleaved.fastq
      gzip clean.SRR14297772_cpe107_1.fastq
      gzip clean.SRR14297772_cpe107_2.fastq
      ```

3. **Verify Trimming Results**:
   Rerun FastQC and run Multiqc on clean fastq files to confirm improvements. 


## Step 2: Downsample Reads (Optional)
Reduce the number of reads in the dataset while preserving the diversity of the sample.

Use **seqtk** sample for random subsampling of reads to a desired percentage or absolute number.

```
seqtk sample -s100 input.fastq 0.1 > downsampled.fastq

#-s100 specifies a random seed for reproducibility.
#0.1 specifies 10% of the total reads (adjust according to needs).
```

Use **bbnorm.sh** in BBMap for read normalization, which reduces redundancy while keeping unique reads.
```
bbnorm.sh in=input.fastq out=downsampled.fastq target=20 min=2

# target=20 controls coverage normalization depth (can be adjusted based on data).
```

## Step 3:  Metagenome Assembly with metaSPAdes

 
**metaSPAdes** is optimized for metagenomic data and assembles reads into contigs, reconstructing genome fragments from complex microbial communities.


1. **Run metaSPAdes**:
   ``` 
   metaspades.py -1 clean.SRR14297772_cpe107_1.fastq.gz -2 clean.SRR14297772_cpe107_1.fastq.gz    -o SRR14297772_cpe107_metaspades_output/ --only-assembler
   ```

   - **What It Does**: metaSPAdes assembles contigs by building a de Bruijn graph adapted for metagenomic data.
   - **Output**: Assembled contigs are saved in the `metaspades_output/` directory.

2. **Review Assembly Results**:
   Inspect `contigs.fasta` in `metaspades_output/` to check contig lengths and quality.


## Step 4: Binning the Contigs

 
Binning groups contigs into bins representing putative genomes. **MetaBAT2** performs binning based on sequence composition and read coverage.


1. **Run MetaBAT2**:
   ``` 
   metabat2 -i metaspades_output/contigs.fasta -o bins_folder/ -m 1500
   ```

   - **What It Does**: MetaBAT2 clusters contigs into bins that represent draft genomes.
   - **Key Option**:
     - `-m 1500`: Sets the minimum contig length to 1500 bp for binning.

2. **Examine Binning Results**:
   The binned genomes are saved in `bins_folder/`, with each bin corresponding to a draft genome.



## Step 5: Quality Assessment of Bins (Optional)

 
Use **CheckM** to evaluate the quality of binned genomes, assessing completeness and contamination based on conserved marker genes.


1. **Run CheckM**:
   ``` 
   checkm lineage_wf bins_folder/ checkm_output/
   ```

   - **What It Does**: CheckM evaluates each bin for genome completeness and contamination.
   - **Interpret Results**: Bins with >90% completeness and <5% contamination are considered high-quality.


## Step 6: Taxonomic Classification of Contigs (MAGs) and/or Reads

 
Classify contigs to identify their taxonomic origin with **Kraken2**, which compares contigs to a taxonomic database.
 

1. **Run Kraken2 to identify microbial diversity**: 
   ``` 
   kraken2 --db kraken2_db --threads 8 --output kraken_output.txt --report kraken_report.txt metaspades_output/contigs.fasta
   ```

   - **What It Does**: Kraken2 assigns taxonomic classifications by matching sequences against a reference database.
   - **Output**: Results are saved in `kraken_output.txt` with a summary report in `kraken_report.txt`.
   - Alternatively, run Kraken 2 on the clean reads.
     ```
     kraken2 --db kraken2_db SRR14297772_cpe107.interleaved.fastq
     ```

## Step 7: Genome Annotation

 
Annotate each bin to identify genes and other genomic features using **Prokka**.

 

1. **Run Prokka for Genome Annotation**:
   ``` 
   prokka --outdir annotation_output --prefix bin_1_annotation bins_folder/bin_1.fa
   ```

   - **What It Does**: Prokka annotates genes and functional elements in each bin.
   - **Output**: Annotations are saved in `annotation_output/`.



## Step 8: AMR Prediction

 
Identify antimicrobial resistance (AMR) genes using **ABRicate**, which screens genomes against known AMR gene databases.

 

1. **Run ABRicate for AMR Prediction**:
   ``` 
   abricate --db resfinder bins_folder/bin_1.fa > abricate_output.txt
   ```

   - **What It Does**: ABRicate searches for known AMR genes by comparing genome sequences against the ResFinder database.
   - **Output**: Results in `abricate_output.txt` list detected AMR genes, their identities, and resistance classes.



## Step 9: Visualization of Taxonomy with Pavian or Krona

 
Visualize taxonomic classifications interactively using **Pavian** or **Krona**:

- **Pavian** provides an interactive web-based interface.
- **Krona** produces circular, hierarchical plots for exploring multi-level taxonomic data.

### Option 1: Visualization with Pavian

1. **Set Up Pavian for Visualization**:
   ``` 
   pavian server
   # or upload the kraken reports on the web server: https://shiny.hiplot.cn/pavian/
   ```

   - **What It Does**: Pavian launches a local server for visualizing taxonomic classifications in your web browser.
   - **Upload**: Load `kraken_report.txt` into Pavian for an interactive visualization.

### Option 2: Visualization with Krona

1. **Convert Kraken2 Output for Krona**:
   ``` 
   cut -f2,3 kraken_output.txt > krona_input.txt
   ```

   - **What It Does**: Extracts only the taxonomic ID and classification from Kraken2 output, creating a format compatible with Krona.

2. **Generate Krona Plot**:
   ``` 
   ktImportTaxonomy krona_input.txt -o krona_output.html
   ```

   - **What It Does**: Krona generates an HTML file (`krona_output.html`) with a multi-level circular plot for exploring taxonomic data.
   - **View the Plot**: Open `krona_output.html` in any web browser.

3. **Exploring the Krona Plot**:
   - Click on different sections to zoom into taxonomic levels.
   - Hover over sections to view specific taxonomic details.




### Exercises
1. Test this workflow on the sample Illumina dataset we provided.
2. Experiment with different assembly parameters in metaSPAdes.
3. Compare taxonomic classifications generated by Kraken2 with other tools, such as Centrifuge.
4. Visualize different datasets in Pavian or Krona to identify any microbial community trends or outliers.



