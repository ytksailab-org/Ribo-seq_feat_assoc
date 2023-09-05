# VeloPro

#Here we developed a pipeline: VeloPro to decipher the association pattern between translation velocity and protein structure features in diverse organisms including bacteria, fungi, protozoa, nematode, plants, insect and mammals.

<h3>Introduction</h3>

Translation velocity plays an important role iin modulating co-translational protein folding and protein functional integrity.Here we provide an easy-to-use and modifiable pipeline that anable to conduct the association analysis between translation velocity and many protein structure features in diverse organisms. This pipeline contains translation velocity quantification procedure and protein structure deterimination procedure in detail with step-by-step explanations.

Software and Installation

fastq-dump


The example procedure
We use publically available data published in pop et al. to illustrat the use of this pipeline.

Quantification of translation velocity from Ribo-seq dataset

1. Download raw FASTQ files from NCBI's sequence Read Archive (SRA) with accession number: SRR1688545 (https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1688545/SRR1688545) which contains ribosome protected footprints of Saccharomyces cerevisiae S288c.
fastq-dump --gzip SRR1688545

2. Ribo-seq quality check.
fastqc -o ./quality_check/ -t 6 ./SRR1688545.fastq.gz

3. For adapter removal and transcriptome mapping
Rscript mapping.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/sequencing_data/ /home/bbian/smORFer_test/MiMB_ribosome_profiling/out /Your/work/path/smORFer_test/MiMB_ribosome_profiling/example_data/genome_data/yeast.transcriptome.fa

4. Statistics and filtering ribosomal reads by length to get an overview of the occurence of each read length
Rscript read_length_distribution.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out

5. Calibration of RPFs
5.1. Split reads by length, 5' assignment for read length of 27 to 31 nucleotides
Rscript split_by_length.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/ 5prime 27-31
5.2 Generate the plot to calibrate the RPFs
Rscript calibration_count_plot.R /home/bbian/smORFer_test/MiMB_ribosome_profiling/out/calibration/split_by_length/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/
5.3 Calibrate the reads according to manually determined offsets using the calibration_5prime_config.csv
Rscript calibrate_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/mapping/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibration_5prime_config.csv /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out
5.4 Merge different read length after calibration into one calibrated bam format file
Rscript merge_reads.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/
5.5 Generate a coverage plot using precisely calibrated reads
Rscript coverage_start_stop.original.R /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibration/calibrated/ /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed /Your/work/path/smORFer_test/MiMB_ribosome_profiling/out/calibrated_coverage

6. Convert the calibrated bam format file to normalized footprints
6.1 Convert the calibrated bam format file to bed format file containing reads counts at each position 
perl countReads3Edge_calibrated.readscount.pl yeast.Ribo.seq_calibrated.bam > inforam_trans.bed positive_strand_info_trans.wig negative_strand_info_trans.wig
6.2 Convert the bed file to nuclotide wave in transcriptome mapping manner
cat inforam_trans.bed | perl -anle '$F[5] eq "+" and print' | bedtools groupby -g 1 -c 2,5 -o collapse,collapse> Yeast_Riboseq_calbrited.nucle.wave
6.3 Convert the transcritpome based nulceotide wise wave to codon wise wave
python nuc2codon.transcriptome.mapping.py yeast.transcriptome.fa yeast_Riboseq_calbrited.nucle.wave yeast_Riboseq_calbrited.codon.wave
6.4 Select most similar protein seq compare to the protein seq from AlphaFold
python multivar_adjust_format_add_codon.py $codon_input/${F1}_bowtie_Riboseq.codon.wave $table_input/${F1}.geneid2name.table.txt $table_input/${F1}.genename2uniprot.tab.txt $input1 | tee $codon_input/${F1}.bowtie.Riboseq.codon.wave.adjusted.format.add.codon
6.5 Remove different length in Ribo-seq data




   









