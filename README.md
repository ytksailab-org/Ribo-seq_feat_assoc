# VeloPro

#Here we developed a pipeline: VeloPro to decipher the association pattern between translation velocity and protein structure features in diverse organisms including bacteria, fungi, protozoa, nematode, plants, insect and mammals.

<h3>Introduction</h3>

Translation velocity plays an important role iin modulating co-translational protein folding and protein functional integrity.Here we provide an easy-to-use and modifiable pipeline that anable to conduct the association analysis between translation velocity and many protein structure features in diverse organisms. This pipeline contains translation velocity quantification procedure and protein structure deterimination procedure in detail with step-by-step explanations.

Software and Installation

fastq-dump
gffread
fastqc



The example procedure
We use publically available data published in pop et al. to illustrat the use of this pipeline.

Quantification of translation velocity from Ribo-seq dataset

1. Ribo-seq data and reference genome preparation\
1.1 Download raw FASTQ files from NCBI's sequence Read Archive (SRA) with accession number: SRR1688545 (https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1688545/SRR1688545) which contains ribosome protected footprints of Saccharomyces cerevisiae S288c.
fastq-dump --gzip SRR1688545

1.2 Download reference genome file and annotation file and generate transcriptome file
gffread yeast.genomic.gff -g yeast.genomic.fa -w yeast.transcriptome.fa

2. Ribo-seq quality check.
fastqc -o ./quality_check/ -t 6 ./SRR1688545.fastq.gz

3. For adapter removal and transcriptome mapping
Rscript mapping.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/sequencing_data/ /home/bbian/smORFer_test/MiMB_ribosome_profiling/out /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/genome_data/yeast.transcriptome.fa

4. Statistics and filtering ribosomal reads by length to get an overview of the occurence of each read length
Rscript read_length_distribution.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out

5. Calibration of RPFs\
5.1. Split reads by length, 5' assignment for read length of 27 to 31 nucleotides\
Rscript split_by_length.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/ 5prime 27-31\
5.2 Generate the plot to calibrate the RPFs\
Rscript calibration_count_plot.R /home/bbian/smORFer_test/MiMB_ribosome_profiling/out/calibration/split_by_length/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/\
5.3 Calibrate the reads according to manually determined offsets using the calibration_5prime_config.csv\
Rscript calibrate_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibration_5prime_config.csv /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out\
5.4 Merge different read length after calibration into one calibrated bam format file\
Rscript merge_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/\
5.5 Generate a coverage plot using precisely calibrated reads\
Rscript coverage_start_stop.original.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibrated_coverage

7. Convert the calibrated bam format file to normalized footprints\
6.1 Convert the calibrated bam format file to bed format file containing reads counts at each position\
perl countReads3Edge_calibrated.readscount.pl yeast.Ribo.seq_calibrated.bam > inforam_trans.bed positive_strand_info_trans.wig negative_strand_info_trans.wig\
6.2 Convert the bed file to nuclotide wave in transcriptome mapping manner\
cat inforam_trans.bed | perl -anle '$F[5] eq "+" and print' | bedtools groupby -g 1 -c 2,5 -o collapse,collapse> Yeast_Riboseq_calbrited.nucle.wave\
6.3 Convert the transcritpome based nulceotide wise wave to codon wise wave\
python nuc2codon.transcriptome.mapping.py yeast.transcriptome.fa yeast_Riboseq_calbrited.nucle.wave yeast_Riboseq_calbrited.codon.wave\
6.4 Select most similar protein seq compare to the protein seq from AlphaFold\
python multivar_adjust_format_add_codon.py yeast_bowtie_Riboseq.codon.wave yeast.geneid2name.table.txt yeast.genename2uniprot.tab.txt yeast | tee yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon\
6.5 Remove different length in Ribo-seq data\
python filter.different.length.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length\
6.6 Filter the low reads Ribo-seq data\
python filter_low_reads_multivarant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length\ yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60\
6.7 Make the scaled the Ribo-seq footprints\
python scale_Riboseq_multivariant.add.codon.py yeast.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60 yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all


Determination of codon usage frequence and protein structure features from AlphaFold

1. Download the pdb format file from AlphaFold database of yeast\
1.1 unzip the file\
tar –xvf UP000002311_559292_YEAST_v4.tar

1.2 unzip all the pdb file\
for i in *.pdb.gz ; do\
  n=`basename $i .gz`\
  gzip -dc $i > $n\
  echo $n\
done > output.txt

2. Extract the positive charge and proline residue by pipinfor tool

pepinfo -sequence $filename -graph ps -outfile $filename.pepinfo

3. Calculate the IDR score for all the proteins in yeast

for pdb in *.pdb ; do\
  fa=`basename $pdb .pdb`.fa\
  python pdb2fa.py $pdb $fa\
done

for f in /home/bbian/Data_all/raw_data/Alpha-fold/yeast/*.fa ; do\
        n=`basename $f`\
        python iupred2a.py $f long >/home/bbian/Data_all/result/IDRs_calculate/iupred2a/yeast/$n.txt\
done

4. Calculate the local relative contact order and local absolute contac order\
perl contactOrder_local_fast_use.pl $filename 1>$filename.contact.order.fast.use 2>$filename.contact.order.fast.use.progress.txt





Association analysis between translation velocity and protein structure features

1. Metagene analysis\
1.1. Merge charge,proline and reads together\
python merge_charge_reads_each_gene.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
1.2 Metagene analysis for positive charge with 95% confidence interval\
python metagene_charge_plot.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
1.3 Metagene analysis for proline residue with 95% confidence interval

2. Comparison of translation velocity between different secondary structure elements in diverse organisms\
2.1 Extract secondary structure for all the genes\
python collect_allresidue_ASA_ss_multvariant.new2.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all yeast.extract_SS_relativeASA.txt
2.2 Make the boxplot\
Rscript 

3. Partial correlation analysis\
3.1 Merge codon usage frequence and footprints together for all the genes\
python codon_usage_analysis_withoutAMBIGUOUS_CODON.py yeast.codon.usage.table.txt yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.2 Merge rASA and footprints together for all the genes\
python ASA_single_gene_correlation_multivariant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.3 Merge IDRs and footprints together for all the genes\
python merge_idrs_reads_nig.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.4 Contact order normalization\
python normalized_CO_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all out 2> log.txt\
3.5 Merge normalized Contact order normalization and footprints for all the genes\
python normalized_CO_counts_merge.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.6 Merge local relative contact order and footprints for all the genes\
python relative_contact_order_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.7 Merge local absoluate contact order and footprints for all the genes\
python absolute_contact_order_multivarant.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.8 Merge all the data together\
python partial_correlation_merge_batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.9 Partial correlation for codon usage frequency\
python spearman_partial_correlation_nig_codonusage.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.10 Partial correlation for ASA\
python spearman_partial_correlation_nig_rASA.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.11 Partial correlation for idrs\
python spearman_partial_correlation_nig_idr.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.12 Partial correlation for normalized contact order\
python spearman_partial_correlation_nig_normalizedCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.13 Partial correlation for local relative contact order\
python spearman_partial_correlation_nig_relativeCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.14 Partial correlation for local absolute contact order\
python spearman_partial_correlation_nig_absoluteCO.batch.py yeast.scaled.bowtie.Riboseq.codon.wave.adjusted.format.add.codon.filter.different.length.higher60.removed.all\
3.15 One sample t-test for all the features\
python ttest_for_all.py



























   









