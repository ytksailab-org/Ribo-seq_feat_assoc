<!DOCTYPE html>
<html>
<head>


</head>

<body>

<h1>VeloPro</h1>

## Here we developed a pipeline: VeloPro to decipher the association pattern between translation velocity and protein structure features in diverse organisms including bacteria, fungi, protozoa, nematode, plants, insect, and mammals

<h1>Introduction</h1>
<p>Translation velocity plays an important role in modulating co-translational protein folding and protein functional integrity. Here we provide an easy-to-use and modifiable pipeline that enables conducting the association analysis between translation velocity and many protein structure features in diverse organisms. This pipeline contains translation velocity quantification procedures and protein structure determination procedures in detail with step-by-step explanations.</p>

<h1>Software and Installation</h1>
<p><strong>Required Software:</strong></p>
<p>This pipeline was tested on NIG supercomputer at RIOS National Institute of Genetics. For most of the tools, we installed via conda or bioconda, and the following software is required:</p>

<ul>
  <li>fastq-dump
    <ul>
      <li>fastq-dump is installed by following the instructions from <a href="https://www.metagenomics.wiki/tools/short-read/ncbi-sra-file-format/sra-tools-install">this link</a></li>
    </ul>
  </li>
  <li>gffread
    <ul>
      <li>Install via bioconda: <code>conda install -c bioconda gffread</code></li>
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

<h1>The Example Procedure</h1>
<p>We use publicly available data published in Pop et al. to illustrate the use of this pipeline.</p>

<h1>Quantification of Translation Velocity from Ribo-seq</h1>

<h2>1. Ribo-seq Data and Reference Genome Preparation</h2>

<h3>1.1 Download Raw FASTQ Files from NCBI's Sequence Read Archive (SRA)</h3>

<p>To obtain the ribosome-protected footprints of Saccharomyces cerevisiae S288c, you can download the raw reads file using following link and then convert raw data to FASTQ format by the following command:</p>

<p>Access the data at this <a href="https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1688545/SRR1688545">link</a>.</p>

<pre><code>fastq-dump --gzip SRR1688545</code></pre>


<h3>1.2 Download Reference Genome File and Annotation File</h3>

<p>Download the reference genome and annotation files for Saccharomyces cerevisiae from the <a href="https://www.ncbi.nlm.nih.gov/genome/?term=Saccharomyces+cerevisiae">NCBI Genome website</a>.</p>

<p>Generate a transcriptome file using the following command:</p>

<pre><code>gffread yeast.genomic.gff -g yeast.genomic.fa -w yeast.transcriptome.fa</code></pre>

<p>These steps will prepare the necessary data for your Ribo-seq analysis.</p>

<h2>2. Ribo-seq Quality Check</h2>

<p>For quality checking the Ribo-seq data, you can use the following command:</p>

<pre><code>fastqc -o ./quality_check/ -t 6 ./SRR1688545.fastq.gz</code></pre>

<h2>3. Adapter Removal and Transcriptome Mapping</h2>

<p>To perform adapter removal and transcriptome mapping, execute the following R script:</p>

<pre><code>Rscript mapping.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/ /home/bbian/smORFer_test/MiMB_ribosome_profiling/out /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/genome_data/yeast.transcriptome.fa</code></pre>

<h2>4. Statistics and Filtering Ribosomal Reads</h2>

<p>Generate an overview of the occurrence of each read length and filter ribosomal reads by length using the following R script:</p>

<pre><code>Rscript read_length_distribution.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out</code></pre>

<h2>5. Calibration of RPFs</h2>

<h3>5.1 Split Reads by Length and Assign 5' End</h3>

<p>Split reads by length and perform 5' assignment for read lengths of 27 to 31 nucleotides:</p>

<pre><code>Rscript split_by_length.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/ 5prime 27-31</code></pre>

<h3>5.2 Generate Calibration Plots</h3>

<p>Generate plots to calibrate the RPFs:</p>

<pre><code>Rscript calibration_count_plot.R /home/bbian/smORFer_test/MiMB_ribosome_profiling/out/calibration/split_by_length/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/</code></pre>

<h3>5.3 Calibrate the Reads According to Manually Determined Offsets</h3>

<p>Calibrate the reads according to manually determined offsets using the calibration_5prime_config.csv file:</p>

<pre><code>Rscript calibrate_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibration_5prime_config.csv /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out</code></pre>

<h3>5.4 Merge Different Read Lengths After Calibration</h3>

<p>Merge different read lengths after calibration into one calibrated BAM format file:</p>

<pre><code>Rscript merge_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/</code></pre>

<h3>5.5 Generate a Coverage Plot Using Precisely Calibrated Reads</h3>

<p>Generate a coverage plot using precisely calibrated reads:</p>

<pre><code>Rscript coverage_start_stop.original.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibrated_coverage</code></pre>

<h2>6. Convert the Calibrated BAM Format File to Normalized Footprints</h2>

<h3>6.1 Convert the Calibrated BAM Format File to BED Format File</h3>

<p>Convert the calibrated BAM format file to a BED format file containing read counts at each position:</p>

<pre><code>perl countReads3Edge_calibrated.readscount.pl yeast.Ribo.seq_calibrated.bam > inforam_trans.bed positive_strand_info_trans.wig negative_strand_info_trans.wig</code></pre>

<h3>6.2 Convert the BED File to Nucleotide Wave in Transcriptome Mapping Manner</h3>

<p>Convert the BED file to nucleotide wave in transcriptome mapping manner:</p>

<pre><code>cat inforam_trans.bed | perl -anle '$F[5] eq "+" and print' | bedtools groupby -g 1 -c 2,5 -o collapse,collapse> Yeast_Riboseq_calbrited.nucle.wave</code></pre>

<h3>6.3 Convert the Transcriptome-Based Nucleotide-Wise Wave to Codon-Wise Wave</h3>

<p>Convert the transcriptome-based nucleotide-wise wave to codon-wise wave:</p>

<pre><code>python nuc2codon.transcriptome.mapping.py yeast.transcriptome.fa yeast_Riboseq_calbrited.nucle.wave yeast_Riboseq_calbrited.codon.wave</code></pre>

<h3>6.4 Select the Most Similar Protein Sequences Compared to AlphaFold</h3>

<p>Select the most similar protein sequences compared to the protein sequences from AlphaFold:</p>

<pre><code>python multivar_adjust_format_add_codon.py yeast_bowtie_Riboseq.codon.wave yeast.geneid2name.table.txt yeast.genename2uniprot.tab.txt yeast | tee yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon</code></pre>

<h3>6.5 Remove Different Lengths in Ribo-seq Data</h3>

<p>Remove different lengths in Ribo-seq data:</p>

<pre><code>python filter.different.length.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length</code></pre>

<h3>6.6 Filter the Low-Reads Ribo-seq Data</h3>

<p>Filter the low-reads Ribo-seq data:</p>

<pre><code>python filter_low_reads_multivarant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60</code></pre>

<h3>6.7 Scale the Ribo-seq Footprints</h3>

<p>Scale the Ribo-seq footprints:</p>

<pre><code>python scale_Riboseq_multivariant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60 yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h1>Determination of Protein Structure Features from AlphaFold</h1>

<h3>1. Download the PDB Format File from AlphaFold Database of Yeast</h3>

<h3>1.1 Unzip the file:</p>

<pre><code>tar –xvf UP000002311_559292_YEAST_v4.tar</code></pre>

<h3>1.2 Unzip All the PDB Files</h3>

<p>Unzip all the PDB files:</p>

<pre><code>for i in *.pdb.gz ; do
  n=`basename $i .gz`
  gzip -dc $i > $n
  echo $n
done > output.txt</code></pre>

<h3>2. Extract the Positive Charge and Proline Residue</h3>

<p>Extract the positive charge and proline residue using the EMBOSS pipinfor tool:</p>

<pre><code>pepinfo -sequence $filename -graph ps -outfile $filename.pepinfo</code></pre>

<h3>3. Calculate the IDR Score for All Proteins in Yeast</h3>

<p>Calculate the IDR score for all proteins in yeast:</p>

<pre><code>for pdb in *.pdb ; do
  fa=`basename $pdb .pdb`.fa
  python pdb2fa.py $pdb $fa
done</code></pre>

<pre><code>for f in /Your/work/path/Alpha-fold/yeast/*.fa ; do
        n=`basename $f`
        python iupred2a.py $f long >/Your/work/path/IDRs_calculate/iupred2a/yeast/$n.txt
done</code></pre>

<h3>4. Calculate the Local Relative Contact Order and Local Absolute Contact Order</h3>

<p>Calculate the local relative contact order and local absolute contact order:</p>

<pre><code>perl contactOrder_local_fast_use.pl $filename 1>$filename.contact.order.fast.use 2>$filename.contact.order.fast.use.progress.txt</code></pre>

<h1>Association Analysis Between Translation Velocity and Protein Structure Features</h1>

<h3>1. Metagene Analysis</h3>

<p>Perform metagene analysis:</p>

<h4>1.1 Merge Charge, Proline, and Reads Together</h4>

<p>Merge charge, proline, and reads together:</p>

<pre><code>python merge_charge_reads_each_gene.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>1.2 Metagene Analysis for Positive Charge with 95% Confidence Interval</h4>

<p>Perform metagene analysis for positive charge with a 95% confidence interval:</p>

<pre><code>python metagene_charge_plot.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>1.3 Metagene Analysis for Proline Residue with 95% Confidence Interval</h4>

<p>Perform metagene analysis for proline residue with a 95% confidence interval:</p>

<!-- Continue with the rest of your instructions -->

<h3>2. Comparison of Translation Velocity Between Different Secondary Structure Elements in Diverse Organisms</h3>

<p>Compare translation velocity between different secondary structure elements in diverse organisms:</p>

<h4>2.1 Extract Secondary Structure for All the Genes</h4>

<p>Extract secondary structure for all the genes:</p>

<pre><code>python collect_allresidue_ASA_ss_multvariant.new2.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all yeast.extract_SS_relativeASA.txt</code></pre>

<h4>2.2 Make the Boxplot</h4>

<p>Make the boxplot:</p>

<pre><code>Rscript </code></pre>

<h3>3. Partial Correlation Analysis</h3>

<p>Perform partial correlation analysis:</p>

<h4>3.1 Merge Codon Usage Frequency and Footprints Together for All the Genes</h4>

<p>Merge codon usage frequency and footprints together for all the genes:</p>

<pre><code>python codon_usage_analysis_withoutAMBIGUOUS_CODON.py yeast.codon.usage.table.txt yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.2 Merge rASA and Footprints Together for All the Genes</h4>

<p>Merge rASA and footprints together for all the genes:</p>

<pre><code>python ASA_single_gene_correlation_multivariant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.3 Merge IDRs and Footprints Together for All the Genes</h4>

<p>Merge IDRs (Intrinsically Disordered Regions) and footprints together for all the genes:</p>

<pre><code>python merge_idrs_reads_nig.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.4 Contact Order Normalization</h4>

<p>Perform contact order normalization:</p>

<pre><code>python normalized_CO_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all out 2> log.txt</code></pre>

<h4>3.5 Merge Normalized Contact Order Normalization and Footprints</h4>

<p>Merge normalized contact order normalization and footprints for all the genes:</p>

<pre><code>python normalized_CO_counts_merge.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.6 Merge Local Relative Contact Order and Footprints</h4>

<p>Merge local relative contact order and footprints for all the genes:</p>

<pre><code>python relative_contact_order_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.7 Merge Local Absolute Contact Order and Footprints</h4>

<p>Merge local absolute contact order and footprints for all the genes:</p>

<pre><code>python absolute_contact_order_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.8 Merge All the Data Together</h4>

<p>Merge all the data together:</p>

<pre><code>python partial_correlation_merge_batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.9 Partial Correlation for Codon Usage Frequency</h4>

<p>Perform partial correlation analysis for codon usage frequency:</p>

<pre><code>python spearman_partial_correlation_nig_codonusage.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.10 Partial Correlation for ASA (Accessible Surface Area)</h4>

<p>Perform partial correlation analysis for ASA:</p>

<pre><code>python spearman_partial_correlation_nig_rASA.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.11 Partial Correlation for IDRs (Intrinsically Disordered Regions)</h4>

<p>Perform partial correlation analysis for IDRs:</p>

<pre><code>python spearman_partial_correlation_nig_idr.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.12 Partial Correlation for Normalized Contact Order</h4>

<p>Perform partial correlation analysis for normalized contact order:</p>

<pre><code>python spearman_partial_correlation_nig_normalizedCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.13 Partial Correlation for Local Relative Contact Order</h4>

<p>Perform partial correlation analysis for local relative contact order:</p>

<pre><code>python spearman_partial_correlation_nig_relativeCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.14 Partial Correlation for Local Absolute Contact Order</h4>

<p>Perform partial correlation analysis for local absolute contact order:</p>

<pre><code>python spearman_partial_correlation_nig_absoluteCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all</code></pre>

<h4>3.15 One Sample t-test for All the Features</h4>

<p>Perform a one-sample t-test for all the features:</p>

<pre><code>python ttest_for_all.py</code></pre>










































   









