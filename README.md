<!DOCTYPE html>
<html>
<head>


</head>

<body>

<h1>VeloPro</h1>
<p>A pipeline for association analysis of translation velocity and protein structure features in diverse organisms</p>

[![视频截图](https://i9.ytimg.com/vi/DA8MI_CFrB8/mqdefault.jpg?v=6538a8d2&sqp=COC156kG&rs=AOn4CLDszNoNVREm2Ra7ks6NiYQtExdOyg)](https://www.youtube.com/watch?v=DA8MI_CFrB8)



<h1>Introduction</h1>
<p>Translation velocity plays an important role in modulating co-translational protein folding and protein functional integrity. Here we provide an easy-to-use and modifiable pipeline that enables the association analysis between translation velocity and many protein structure features in diverse organisms. This pipeline contains translation velocity quantification procedures and protein structure feature determination procedures with step-by-step explanations.</p>

<h1>Installation</h1>
<p><strong>Required Software:</strong></p>
<p>This pipeline was tested on a Linux-based system with Python 3.7.2 and R 3.5.2. The following packages are required, most of them can be installed via conda or bioconda.</p>

<ul>
  <li>fastq-dump
    <ul>
      <li>Install following the instructions <a href="https://www.metagenomics.wiki/tools/short-read/ncbi-sra-file-format/sra-tools-install">here</a></li>
    </ul>
  </li>
  <li>gffread
    <ul>
      <li>Install via bioconda: <code>conda install -c bioconda gffread</code></li>
    </ul>
  </li>
  <li>fastp
    <ul>
      <li>Install via bioconda: <code>conda install -c bioconda fastp</code></li>
    </ul>
  </li>
  <li>fastqc
    <ul>
      <li>Install via conda: <code>conda install fastqc</code></li>
    </ul>
  </li>
  <li>cutadapt
    <ul>
      <li>Install via conda: <code>conda install cutadapt</code></li>
    </ul>
  </li>
  <li>bowtie
    <ul>
      <li>Install via conda: <code>conda install bowtie</code></li>
    </ul>
  </li>
    <li>EMBOSS
    <ul>
      <li>Install via bioconda: <code>conda install -c bioconda emboss</code></li>
    </ul>
  </li>
  <li>samtools
    <ul>
      <li>Install via conda: <code>conda install samtools</code></li>
    </ul>
  </li>
  <li>bedtools
    <ul>
      <li>Install via conda: <code>conda install bedtools</code></li>
    </ul>
  </li>
</ul>

<h1>Example Procedure</h1>
<p>We use publicly available Ribo-seq data of Saccharomyces cerevisiae in <a href="https://doi.org/10.15252/msb.20145524">Pop et al. 2014</a> to illustrate the use of this pipeline. Other organisms can be analyzed in a similar way given the corresponding Ribo-seq data.</p>

<h1>Quantification of Translation Velocity from Ribo-seq data</h1>
<p>The R scripts used in this step were those modified from <a href="https://github.com/AlexanderBartholomaeus/MiMB_ribosome_profiling"> MiMB_ribosome_profiling pipeline</a> in <a href="https://doi.org/10.1007/978-1-0716-1150-0_12">Bartholomaus et al. 2021</a>.</p>

<h2>1. Ribo-seq Data and Reference Genome Preparation</h2>

<h3>1.1 Download Raw FASTQ Files from NCBI Sequence Read Archive (SRA)</h3>

<p>To obtain the ribosome-protected footprints (RPF) of Saccharomyces cerevisiae S288c, you can download the raw read file using the following link and then convert the raw data to FASTQ format by the following command:</p>

<p>Access the data at this <a href="https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1688545/SRR1688545">link</a>.</p>

<pre><code>fastq-dump --gzip SRR1688545</code></pre>


<h3>1.2 Download Reference Genome File and Annotation File</h3>

<p>Download the reference genome and annotation files for Saccharomyces cerevisiae from the <a href="https://www.ncbi.nlm.nih.gov/genome/?term=Saccharomyces+cerevisiae">NCBI Genome website</a>.</p>

<p>Generate a reference transcriptome file using the following command:</p>

<pre><code>gffread yeast.genomic.gff -g yeast.genomic.fa -w yeast.transcriptome.fa</code></pre>

<p>These steps will prepare the necessary data for your Ribo-seq analysis.</p>

<h2>2. Ribo-seq Quality Check</h2>

<p>For quality check of Ribo-seq data, you can use the following command:</p>

<pre><code>fastqc -o ./quality_check/ -t 6 ./SRR1688545.fastq.gz</code></pre>

<h2>3. Adapter Removal and Transcriptome Mapping</h2>

<p>To perform adapter trimming and transcriptome mapping, execute the following R script (Note: this step will take more than 10 minutes):</p>

<pre><code>Rscript mapping.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/ /home/bbian/smORFer_test/MiMB_ribosome_profiling/out /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/genome_data/yeast.transcriptome.fa</code></pre>

<p> Another option for adapters removal by using fastp:</p>

<pre><code> fastp --adapter_sequence AGATCGGAAGAGCACACGTCT </code></pre>

<h2>4. Statistics and Filtering Ribosomal Reads</h2>

<p>Generate an overview of the occurrence of each read length, and filter mapped reads by lengths using the following R script:</p>

<pre><code>Rscript read_length_distribution.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out</code></pre>

<h2>5. Calibration of RPFs</h2>

<h3>5.1 Determine the genes with high expression levels in yeast dataset </h3>

<pre><code>Rscript highest_expressed_genes.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/yeast_Riboseq.sort.bam /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/genome_data/yeast.genome.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/ 0.1</code></pre>

<h3>5.2 Prepare annotation files to include 50 nucleotides upstream and downstream from start and stop codon, respectively</h3>
  
<p> Prepare the BED file:</p>

<pre><code>Rscript prepare_coverage.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes.bed </code></pre>

<h3>5.3 Split Reads by Length and Assign 5' End</h3>

<p>As suggested by the original pipeline <a href="https://github.com/AlexanderBartholomaeus/MiMB_ribosome_profiling"> MiMB_ribosome_profiling pipeline</a> in <a href="https://doi.org/10.1007/978-1-0716-1150-0_12">Bartholomaus et al. 2021</a>, the offset values for calibration may vary depending on different read lengths. Thus, we first split the reads based on their lengths. To split reads by length and perform 5' assignment for read lengths of 27 to 31 nucleotides:</p>

<pre><code>Rscript split_by_length.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/ 5prime 27-31</code></pre>

<h3>5.4 Generate Calibration Plots</h3>

<p>Generate plots to calibrate the RPFs:</p>

<pre><code>Rscript calibration_count_plot.R /home/bbian/smORFer_test/MiMB_ribosome_profiling/out/calibration/split_by_length/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/</code></pre>
<img src="https://github.com/ytksailab-org/VeloPro/blob/main/yeast_to_calibrate_plot.png" alt="GitHub Logo">


<h3>5.5 Calibrate the Reads According to Manually Determined Offsets</h3>

<p>Calibrate the reads according to manually determined offsets using the calibration_5prime_config.csv file:</p>

<pre><code>Rscript calibrate_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibration_5prime_config.csv /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out</code></pre>

<h3>5.6 Merge Different Read Lengths After Calibration</h3>

<p>Merge different read lengths after calibration into one calibrated BAM format file:</p>

<pre><code>Rscript merge_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/</code></pre>

<h3>5.7 Generate a Coverage Plot Using Precisely Calibrated Reads</h3>

<pre><code>Rscript coverage_start_stop.original.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibrated_coverage</code></pre>
<img src="https://github.com/ytksailab-org/VeloPro/blob/main/yeast_calibrated_plot.png" alt="GitHub Logo">

<h2>6. Calculate Scaled Footprints from the Merged Calibrated BAM File</h2>

<h3>6.1 Convert the Merged Calibrated BAM File to BED File</h3>

<p>Convert the merged calibrated BAM file to a BED file containing read counts at each position:</p>

<pre><code>perl countReads3Edge_calibrated.readscount.pl yeast.Ribo.seq_calibrated.bam > inforam_trans.bed positive_strand_info_trans.wig negative_strand_info_trans.wig</code></pre>

<h3>6.2 Convert the BED File to Nucleotide-Wise Footprint Wave </h3>

<p>Most reads should be mapped to the sense strand as the Ribo-seq library is strand-specific. Therefore, in transcriptome mapping, reads mapped to the minus strand should be considered as errors. Consequently, we only select reads marked with a '+' for analysis. Convert the BED file to the nucleotide-wise wave in each transcript:</p>

<pre><code>cat inforam_trans.bed | perl -anle '$F[5] eq "+" and print' | bedtools groupby -g 1 -c 2,5 -o collapse,collapse> Yeast_Riboseq_calbrited.nucle.wave</code></pre>

<h3>6.3 Convert the Nucleotide-Wise Footprint Wave to Codon-Wise footprint Wave</h3>

<p>Convert the nucleotide-wise wave to the codon-wise wave by summing the footprints in each codon:</p>

<pre><code>python nuc2codon.transcriptome.mapping.py yeast.transcriptome.fa yeast_Riboseq_calbrited.nucle.wave yeast_Riboseq_calbrited.codon.wave</code></pre>

<h3>6.4 Select one isoform for each gene based on the Similarity to Sequences in AlphaFold database</h3>

<p>Select the most similar isoform for each gene by comparing to protein sequences in AlphaFold database:</p>

<pre><code>python multivar_adjust_format_add_codon.py yeast_bowtie_Riboseq.codon.wave yeast.geneid2name.table.txt yeast.genename2uniprot.tab.txt yeast | tee yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon</code></pre>

<h3>6.5 Guarantee the Length of the Selected Isoform is the Same as that in AlphaFold database</h3>

<p>Remove the genes with different lengths of each line in codon wise wave:</p>

<pre><code>python filter.different.length.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length</code></pre>

<h3>6.6 Filter the Genes with Low-Sequencing Coverage in Codon-wise Wave Data</h3>

<p>Exclude the genes with sequencing coverage less than 60% in codon-wise wave data:</p>

<pre><code>python filter_low_reads_multivarant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60</code></pre>

<h3>6.7 Scale the Ribo-seq Footprints</h3>

<p>The codon-wise wave are scaled by the mean coverage within each gene:</p>

<pre><code>python scale_Riboseq_multivariant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60 yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h1>Determination of Protein Structure Features from AlphaFold database</h1>

<h3>1. Download all the proteins in the PDB Format from AlphaFold Database for Yeast</h3>
<p>Access the data at this <a href="https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar">link</a>.</p>

<h3>1.1 Unzip the file:</h3>

<pre><code>tar –xvf UP000002311_559292_YEAST_v4.tar</code></pre>

<h3>1.2 Unzip All the PDB Files</h3>

<p>Unzip all the PDB files:</p>

<pre><code>for i in *.pdb.gz ; do
  n=`basename $i .gz`
  gzip -dc $i > $n
  echo $n
done > output.txt</code></pre>

<h3>2. Extract the Positive Charge and Proline Residue for All Proteins in Yeast</h3>

<p>Extract the positive charge and proline residue using the EMBOSS pepinfo tool:</p>

<pre><code>pepinfo -sequence $filename -graph ps -outfile $filename.pepinfo</code></pre>

<h3>3. Calculate the IDR (is defined as sequence segments that do not form specific 3D structures) Score for All Proteins in Yeast</h3>

<pre><code>for pdb in *.pdb ; do
  fa=`basename $pdb .pdb`.fa
  python pdb2fa.py $pdb $fa
done</code></pre>

<pre><code>for f in /Your/work/path/Alpha-fold/yeast/*.fa ; do
        n=`basename $f`
        python iupred2a.py $f long >/Your/work/path/IDRs_calculate/iupred2a/yeast/$n.txt
done</code></pre>

<h3>4. Calculate the Local Relative Contact Order and Local Absolute Contact Order for All Proteins</h3>
<p>The perl script used in this step was those modified from <a href="https://depts.washington.edu/bakerpg/contact_order/"> Determine a Protein's Contact Order</a> in <a href="https://doi.org/10.1006/jmbi.1998.1645">Plaxco et al. 1998</a>.</p>

<pre><code>perl contactOrder_local_fast_use.pl $filename 1>$filename.contact.order.fast.use 2>$filename.contact.order.fast.use.progress.txt</code></pre>

<h1>Association Analysis Between Translation Velocity and Protein Structure Features</h1>

<h3>1. Metagene Analysis for Positive Charge and Proline Residue</h3>

<h4>1.1 Merge Charge/Proline Residue and Footprints Together</h4>
<pre><code>python merge_charge_reads_each_gene.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>1.2 Metagene Analysis for Positive Charge with 95% Confidence Interval</h4>
<pre><code>python metagene_charge_plot.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>1.3 Metagene Analysis for Proline Residue with 95% Confidence Interval</h4>
<pre><code>python metagene_proline_plot.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>
<img src="https://github.com/ytksailab-org/VeloPro/blob/main/yeast_proline_metagene_plot.png" alt="GitHub Logo">

<h3>2. Comparison of Translation Velocity Between Different Secondary Structure Elements in Diverse Organisms</h3>

<h4>2.1 Extract Secondary Structure for All the Proteins in Yeast (Note: this step may take more than 10 minutes)</h4>
<pre><code>python collect_allresidue_ASA_ss_multvariant.new2.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all yeast.extract_SS_relativeASA.txt</code></pre>

<h4>2.2 Make the Boxplot</h4>
<pre><code>Rscript Boxplot_for_SS.R </code></pre>
In this script, we hard-code the input and output filenames. Please modify the script as necessary for the use in your environment. An example box plot looks like below:
<img src="https://github.com/ytksailab-org/VeloPro/blob/main/yeast.SS.Rplot.png" alt="GitHub Logo">

<h3>3. Partial Correlation Analysis</h3>
<h4>3.1 Merge Codon Usage Frequency and Footprints Together for All the Genes</h4>
<pre><code>python codon_usage_analysis_withoutAMBIGUOUS_CODON.py yeast.codon.usage.table.txt yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.2 Merge rASA (A measure of the solvent molecule contact of amino acid residues in a protein) and Footprints Together for All the Genes</h4>
<pre><code>python ASA_single_gene_correlation_multivariant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.3 Merge IDRs (Intrinsically Disordered Regions)and Footprints Together for All the Genes</h4>
<pre><code>python merge_idrs_reads_nig.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.4 Local Contact Order (A measure of the locality of the inter-amino acid contacts in the protein's native state tertiary structure) Normalization </h4>
<pre><code>python normalized_CO_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all out 2> log.txt</code></pre>

<h4>3.5 Merge Normalized Local Contact Order and Footprints for All the Genes</h4>
<pre><code>python normalized_CO_counts_merge.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.6 Merge Local Relative Contact Order and Footprints for All the Genes</h4>
<pre><code>python relative_contact_order_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.7 Merge Local Absolute Contact Order and Footprints for All the Genes</h4>
<pre><code>python absolute_contact_order_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.8 Merge All the Data Together</h4>
<pre><code>python partial_correlation_merge_batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.9 Partial Correlation for Codon Usage Frequency</h4>
<pre><code>python spearman_partial_correlation_nig_codonusage.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.10 Partial Correlation for ASA (Accessible Surface Area)</h4>
<pre><code>python spearman_partial_correlation_nig_rASA.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.11 Partial Correlation for IDRs (Intrinsically Disordered Regions)</h4>
<pre><code>python spearman_partial_correlation_nig_idr.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.12 Partial Correlation for Normalized Local Contact Order</h4>
<pre><code>python spearman_partial_correlation_nig_normalizedCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.13 Partial Correlation for Local Relative Contact Order</h4>
<pre><code>python spearman_partial_correlation_nig_relativeCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.14 Partial Correlation for Local Absolute Contact Order</h4>
<pre><code>python spearman_partial_correlation_nig_absoluteCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.15 One Sample t-test for All the Features</h4>
<pre><code>python ttest_for_all.py</code></pre>










































   









