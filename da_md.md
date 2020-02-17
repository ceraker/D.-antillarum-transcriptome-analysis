# __*D. antillarum* transcriptome analysis__
#### Author: Cassie Raker
#### Last updated: February 16, 2020

#### Data uploaded and analyzed on KITT

###### accessing KITT account
```
ssh -p -Y 2292 craker@kitt.uri.edu
```
###### Prepare a working directory with data
Make a directory
```
mkdir diadema
cd diadema
```
Create a conda environment
```
conda create -n diadema
conda activate diadema
```
Create symbolic links to data
```
ln /home/Shared_Data/Prada_Data/*.fastq.gz .
```

##### Quality Analysis using FastQC
FastQC is used to identify any problems with reads, such as low quality or contamination from TruSeq adapters.
```
conda install fastqc
mkdir fastqc_results
cd fastqc_results
fastqc /home/craker/diadema/*.fastq.gz
```
Generate a multiqc report
```
conda install multiqc
multiqc .
scp -P 2292 craker@kitt.uri.edu:/home/craker/diadema/fastqc_results/multiqc_report.html /Users/cassieraker/Desktop/Urchin2/fastqc
```
Mean Quality Scores from MultiQC report
![mean quality scores](https://github.com/jpuritz/BIO_594_2019/blob/master/Final_Assignment/Raker_final_assignment/Images/fastqc_meanqualityscores.png)

###### Summarize results
```
cat */summary.txt > /home/craker/diadema/fastqc_results/fastqc_summaries.txt
grep FAIL fastqc_summaries.txt | wc -l
```
47 FastQC tests failed

##### Quality control using fastp
fastp is used to trim low quality reads, remaining adapter sequences, etc.
```
conda install fastp
```
Run fastp on one file
```
fastp --in1 DA-HI-A_S79_L007_R1_001.fastq.gz \
--out1 DA-HI-A_S79_L007_R1_001.fastq.gz.clean \
--cut_front 20 --cut_tail 20 --cut_window_size 5 \
--cut_mean_quality 15 -q 15 -w 16 \
--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
Loop fastp to run on all files
```
for filename in *.fastq.gz
> do
> fastp --in1 ${filename} --out1 ${filename}.clean \
> --cut_front 20 --cut_tail 20 --cut_window_size 5 \
> --cut_mean_quality 15 -q 15 -w 16 --trim_front1 14 \
> --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
> done
```
Quality scores post fastp trimming

![fastp quality scores](https://github.com/jpuritz/BIO_594_2019/blob/master/Final_Assignment/Raker_final_assignment/Images/fastp_quality.png)

Base contents post fastp trimming

![fastp base contents](https://github.com/jpuritz/BIO_594_2019/blob/master/Final_Assignment/Raker_final_assignment/Images/fastp_base_contents.png)

##### de novo transcriptome with Trinity
Without a reference transcriptome, Trinity can be used to do *de novo* transcriptome analysis.

Install Trinity
```
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.4.zip
unzip Trinity-v2.8.4.zip
cd trinityrnaseq-Trinity-v2.8.4/
make
make plugins
export TRINITY_HOME=/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4
```

test Trinity
```
cd sample_data/test_Trinity_Assembly/
./runMe.sh
```

run Trinity on all clean files, using --full_cleanup option to delete uneeded files at the end of the analysis to save space
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Trinity \
--seqType fq --max_memory 50G \
--single DA-HI-A_S79_L007_R1_001_clean.fq,\
DA-HI-B_S80_L007_R1_001_clean.fq,\
DA-HI-C_S81_L007_R1_001_clean.fq,\
DA-LOW-A_S72_L007_R1_001_clean.fq,\
DA-LOW-B_S73_L007_R1_001_clean.fq,\
DA-LOW-C_S74_L007_R1_001_clean.fq,\
DA-LOW-D_S75_L007_R1_001_clean.fq,\
DA-MED-A_S76_L007_R1_001_clean.fq,\
DA-MED-B_S77_L007_R1_001_clean.fq,\
DA-MED-D_S78_L007_R1_001_clean.fq,\
OADA0006_S16_L002_R1_001_clean.fq,\
OADA0049_S25_L003_R1_001_clean.fq,\
OADA0058_S19_L002_R1_001_clean.fq,\
OADA0071_S53_L005_R1_001_clean.fq,\
OADA0081_S67_L006_R1_001_clean.fq,\
OADA0085_S24_L002_R1_001_clean.fq,\
OADA0101_S18_L002_R1_001_clean.fq,\
OADA0102_S33_L003_R1_001_clean.fq,\
OADA0116_S26_L003_R1_001_clean.fq,\
OADA0139_S17_L002_R1_001_clean.fq,\
OADA0174_S20_L002_R1_001_clean.fq \
--CPU 20 --full_cleanup
```

calculate number of clusters
```
cat da.trinity.fasta | grep '>' | wc -l
645371 clusters
```

##### Transdecoder and CD-Hit
Both of these programs are used to further group clusters.

###### Transdecoder
```
conda install transdecoder
TransDecoder.LongOrfs -t da.trinity.fasta
TransDecoder.Predict -t da.trinity.fasta
```
###### CD-hit
```
conda install cd-hit
cd-hit-est -i ./da.trinity.fasta.transdecoder_dir/longest_orfs.cds -o cdhit89 -n 8 -c 0.89 -M 96000 -T 20 -d 0
```
calculate number of clusters
```
cat cdhit89 | grep '>' | wc -l
65628 clusters
```
run TransDecoder on cd-hit output
```
TransDecoder.LongOrfs -t cdhit89
TransDecoder.Predict -t cdhit89
```
rename best fasta file (so it's easy to find later)
```
mv cdhit89.transdecoder.pep.fasta diadema.best.fasta
```

##### Blast files using Trinotate
Trinotate is used to annotate files.

Install Trinotate
```
wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.1.1.zip
unzip Trinotate-v3.1.1.zip
cd Trinotate-Trinotate-v3.1.1/
chmod 755 Trinotate
/home/craker/data/Trinotate-Trinotate-v3.1.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

BLAST
```
blastx -query /home/craker/data/Trinity.fasta -db /home/craker/data/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
blastp -query /home/craker/data/cdhit89.transdecoder.pep -db /home/craker/data/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
```

##### OrthoFinder
The *Diadema antillarum* data was compared to sequences from *Strongylocentrotous purpuratus*. Data for *S. purpuratus* was downloaded from EchinoBase.org.

Download *Diadema antillarum* data to personal computer, and then moved to KITT
```
scp -P 2292 ~/Desktop/Diadema/SPU_peptide.fasta.zip craker@kitt.uri.edu:/home/craker/diadema/
unzip SPU_peptide.fasta.zip
```

Install OrthoFinder
```
wget https://github.com/davidemms/OrthoFinder/releases/download/v2.3.1-beta/OrthoFinder-2.3.1_source.tar.gz
tar xvzf OrthoFinder-2.3.1_source.tar.gz
```
Test OrthoFinder
```
cd OrthoFinder-2.3.1_source/orthofinder
orthofinder -f ExampleData -S diamond
```
Copy fasta files to OrthoFinder directory
```
mkdir DAfastas
cp /home/craker/diadema/diadema.best.fasta ./DAfastas/
cp /home/craker/diadema/ortho_fastas/SPU_peptide.fasta ./DAfastas/
```
###### Run OrthoFinder
```
orthofinder -f ./DAfastas/ -t 20
```

Generate a fasta file for each orthogroup by running:
orthofinder -fg RESULTS_DIR -M msa -os
where RESULTS_DIR is the directory containing your orthogroups results files (e.g. Orthogroups.csv)
```
orthofinder -fg /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/DAfastas/Results_May02_1 -M msa -os
```
Combine individual fasta files into one fasta file
```
cat *.fa > ortho.fasta
```

##### Filtering out only overlapping orthogroups from original fasta file
In this analysis, we want to focus on orthogroups found in common between *D. antillarum* and *S. purpuratus*


Navigate to .csv files containing overlapping orthogroups and examine file
```
cd /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/DAfastas/Results_May02_1/Orthologues_May03/Orthologues/Orthologues_diadema.best
nano diadema.best__v__SPU_peptide.csv
```
headers from original Trinity fasta file are all located in the second column
Use cut to create a text file with only this column
```
cut -f2 diadema.best__v__SPU_peptide.csv > onecolumn.txt
```
Copy onecolumn.txt to main directory for ease of further analysis and navigate back up to that directory
```
cp onecolumn.txt /home/craker/diadema/
cd /home/craker/diadema/
```
###### Use bbtools script filterbyname.sh and onecolumn.txt to create a new fasta file from da.trinity.fasta containing only the desired transcripts
Use include=t so output includes instead of excludes files named in onecolumn.txt
Use substring=t since Orthofinder added .p1.p1 to the end of headers
```
conda install bbtools
filterbyname.sh in=da.trinity.fasta names=onecolumn.txt out=da.filtered.fasta \
include=t substring=t
```
Use the da.filtered.fasta file for further downstream analysis

##### Use Trinotate to annotate the filtered fasta file
Install other necessary programs
```
conda install hmmer
conda install perl-dbd-sqlite
```
Run TransDecoder on the filtered fasta file
```
TransDecoder.LongOrfs -t da.filtered.fasta
TransDecoder.Predict -t da.filtered.fasta
```
Blast the peptide file
```
blastp -query /home/craker/diadema/da.filtered.fasta.transdecoder.pep -db /home/craker/diadema/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
```
Run HMMER
```
hmmscan --cpu 20 --domtblout TrinotatePFAM.out Pfam-A.hmm da.filtered.fasta.transdecoder.pep > pfam.log
```
Load files into an SQLite database
```
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init --gene_trans_map da.filtered.fasta.gene_trans_map --transcript_fasta da.filtered.fasta --transdecoder_pep da.filtered.fasta.transdecoder.pep
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.filtered.outfmt6
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.filtered.outfmt6
```
Generate a gene/transcript relationship file using TRINITY
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl da.filtered.fasta > da.filtered.fasta.gene_trans_map
```
Generate a Trinotate annotation report
```
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > da.filtered.trinotate_annotation_report.xls
```

##### Trinity transcript quantification
Necessary to perform further downstream analysis.

Install necessary programs
```
conda install salmon
conda install rsem
conda install bowtie2
```
Created a tab delimited text file of samples (da_samples.txt)

###### Running the analysis with the original fasta file:
Prepare the reference and run alignment and abundance estimation
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl \
--transcripts diadema.best.fasta --est_method salmon --aln_method bowtie2 \
--prep_reference --trinity_mode --samples_file da_samples.txt --seqType fq
```

Create matrices
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl \
--est_method salmon --gene_trans_map none --out_prefix salmon --name_sample_by_basedir \
pH_high_rep1/quant.sf pH_high_rep2/quant.sf pH_high_rep3/quant.sf pH_high_rep4/quant.sf \
pH_high_rep5/quant.sf pH_high_rep6/quant.sf pH_high_rep7/quant.sf pH_medx_rep1/quant.sf \
pH_medx_rep2/quant.sf pH_medx_rep3/quant.sf pH_medx_rep4/quant.sf pH_medx_rep5/quant.sf \
pH_medx_rep6/quant.sf pH_lowx_rep1/quant.sf pH_lowx_rep2/quant.sf pH_lowx_rep3/quant.sf \
pH_lowx_rep4/quant.sf pH_lowx_rep5/quant.sf pH_lowx_rep6/quant.sf pH_lowx_rep7/quant.sf \
pH_lowx_rep8/quant.sf

```

These matrix files, specifically the salmon.isoforms.count.matrix file, will be used for gene expression analysis.


Filter the results to only include most differentially espressed isoforms
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/filter_low_expr_transcripts.pl \
--m salmon.isoform.TPM.not_cross_norm --t da.trinity.fasta \
--min_expr_any 0 --highest_iso_only --trinity_mode \
> filter_hio.fasta

/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/filter_low_expr_transcripts.pl \
--m salmon.isoform.TPM.not_cross_norm --t da.trinity.fasta \
--min_expr_any 1 --min_pct_dom_iso 1 --trinity_mode \
> filter_expr1.fasta
```
Can create matrices using these results as well, just re-run above scripts with filter_hio.fasta and filter_expr1.fasta in place of da.trinity.fasta


###### Create matrix from da.filtered.fasta
Prepare the reference and align
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl \
--transcripts da.filtered.fasta --est_method salmon --aln_method bowtie2 \
--prep_reference --trinity_mode --samples_file da_samples.txt --seqType fq
```

Create matrices
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl \
--est_method salmon --gene_trans_map none --out_prefix salmon --name_sample_by_basedir \
pH_high_rep1/quant.sf pH_high_rep2/quant.sf pH_high_rep3/quant.sf pH_high_rep4/quant.sf \
pH_high_rep5/quant.sf pH_high_rep6/quant.sf pH_high_rep7/quant.sf pH_medx_rep1/quant.sf \
pH_medx_rep2/quant.sf pH_medx_rep3/quant.sf pH_medx_rep4/quant.sf pH_medx_rep5/quant.sf \
pH_medx_rep6/quant.sf pH_lowx_rep1/quant.sf pH_lowx_rep2/quant.sf pH_lowx_rep3/quant.sf \
pH_lowx_rep4/quant.sf pH_lowx_rep5/quant.sf pH_lowx_rep6/quant.sf pH_lowx_rep7/quant.sf \
pH_lowx_rep8/quant.sf
```



##### DESeq2
DESeq2 can actually be run from Trinity scripts, as long as the proper R packages are installed. To install the proper R packages:
```
% R
 > source("http://bioconductor.org/biocLite.R")
 > biocLite('edgeR')
 > biocLite('limma')
 > biocLite('DESeq2')
 > biocLite('ctc')
 > biocLite('Biobase')
 > install.packages('gplots')
 > install.packages('ape')
 > biocLite('qvalue')
 > biocLite('fastcluster')
 ```
 DESeq2 is then run using this script:
 ```
 $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl
 ```
###### Using my files generated from original diadema.best.fasta:
 ```
 /home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/run_DE_analysis.pl \
 --matrix salmon.isoform.counts.matrix \
 --method DESeq2 --samples_file da_samples.txt
 ```
Results of this script:
```
salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.count_matrix
salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.DE_results
salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.DE_results.MA_n_Volcano.pdf
salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.Rscript
salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.count_matrix
salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.DE_results
salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.DE_results.MA_n_Volcano.pdf
salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.Rscript
salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.count_matrix
salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.DE_results
salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.DE_results.MA_n_Volcano.pdf
salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.Rscript
```

MA and Volcano plots of entire transcript examining high and low pH treatments
![full transcript volcano plot](https://github.com/ceraker/D.-antillarum-transcriptome-analysis/blob/master/figures/salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.DE_results.MA_n_Volcano.png)


To extract and cluster differentially expressed transcripts, I would use this script:
```
$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl
```
Copy sample text file into the DESeq2_outdir directory, move there, and run script like this:
```
cp da_samples.txt /home/craker/DESeq2_outdir
cd DESeq2_outdir

/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix /home/craker/diadema/salmon.isoform.TMM.EXPR.matrix \
--samples da_samples.txt

```
Results of this script include a correlation heatmap:

```
DE_feature_counts.P0.001_C2.matrix
diffExpr.P0.001_C2.matrix
diffExpr.P0.001_C2.matrix.log2.centered.dat
diffExpr.P0.001_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf
diffExpr.P0.001_C2.matrix.log2.centered.sample_cor.dat
diffExpr.P0.001_C2.matrix.log2.centered.sample_cor_matrix.pdf
diffExpr.P0.001_C2.matrix.R
diffExpr.P0.001_C2.matrix.RData
```
Correlation heatmap on unfiltered samples
![unfiltered heatmap](https://github.com/ceraker/D.-antillarum-transcriptome-analysis/blob/master/figures/diffExpr.P0.001_C2.matrix.log2.centered.genes_vs_samples_heatmap.png)

DESeq2 analysis was then run on the two filtered matrices.

Heatmaps for two filtered matrices:

![Expr1 heatmap](https://github.com/ceraker/D.-antillarum-transcriptome-analysis/blob/master/figures/filter_expr1.diffExpr.P0.001_C2.matrix.log2.centered.genes_vs_samples_heatmap.png)
![HIO heatmap](https://github.com/ceraker/D.-antillarum-transcriptome-analysis/blob/master/figures/filter_hio.diffExpr.P0.001_C2.matrix.log2.centered.genes_vs_samples_heatmap.png)

DESeq2 was attempted on te filtered fasta file, but no significantly differentially expressed genes were found.
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/run_DE_analysis.pl \
 --matrix salmon.isoform.counts.matrix \
 --method DESeq2 --samples_file da_samples.txt

 /home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/analyze_diff_expr.pl \
 --matrix /home/craker/diadema/salmon.isoform.TMM.EXPR.matrix \
 --samples da_samples.txt

 Error, no differentially expressed transcripts identified at cuttoffs: P:0.001, C:2 at $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl line 203.

```

From here we need to determine how we want to proceed, and if we want to pursue futher analysis on the unfiltered results. 
