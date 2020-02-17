## Gene Expression Analysis with only larval urchins
#### Last updated: February 16, 2020

run Trinity on all clean files, using --full_cleanup option to delete uneeded files at the end of the analysis to save space

```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Trinity \
--seqType fq --max_memory 50G \
--single /home/craker/diadema/DA-HI-A_S79_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-HI-C_S81_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-HI-D_S82_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-LOW-A_S72_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-LOW-B_S73_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-LOW-C_S74_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-LOW-D_S75_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-MED-A_S76_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-MED-B_S77_L007_R1_001_clean.fq,\
/home/craker/diadema/DA-MED-D_S78_L007_R1_001_clean.fq \
--CPU 20 --full_cleanup --output trinity_larvae_out_dir
```


rename fasta file
```
mv trinity_larvae_out_dir.Trinity.fasta da.larvae.fasta
mv trinity_larvae_out_dir/ da.larvae_out_dir/
mv trinity_larvae_out_dir.Trinity.fasta.gene_trans_map da.larvae.fasta.gene_trans_map
```

calculate number of clusters
```
cat da.larvae.fasta | grep '>' | wc -l
```
389705 clusters

#### TransDecoder and CD-hit
```
TransDecoder.LongOrfs -t da.larvae.fasta
TransDecoder.Predict -t da.larvae.fasta
cd-hit-est -i ./da.larvae.fasta.transdecoder_dir/longest_orfs.cds -o cdhit89 -n 8 -c 0.89 -M 96000 -T 20 -d 0
```

calculate number of clusters
```
cat cdhit89 | grep '>' | wc -l
```
45177 clusters

run TransDecoder on cd-hit output
```
TransDecoder.LongOrfs -t cdhit89
TransDecoder.Predict -t cdhit89
```

rename best fasta file (so it's easy to find later)
```
mv cdhit89.transdecoder.pep da.larvae.best.fasta
```

#### Trinotate

BLAST
```
blastx -query /home/craker/diadema/da.larvae.fasta -db /home/craker/diadema/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
blastp -query /home/craker/diadema/da.larvae.transdecoderCDhit_dir/cdhit89.transdecoder.pep -db /home/craker/diadema/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
```

#### Orthofinder

Copy fasta files to OrthoFinder directory
```
cd OrthoFinder-2.3.1_source/orthofinder/
mkdir da.larvae.fastas
cp /home/craker/diadema/da.larvae.transdecoderCDhit_dir/da.larvae.best.fasta /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.larvae.fastas
cp /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/DAfastas/SPU_peptide.fasta /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.larvae.fastas
```

run OrthoFinder (remember to deactivate conda environment)
```
orthofinder -f ./da.larvae.fastas/ -t 20
```

Generate a fasta file for each orthogroup by running:
orthofinder -fg RESULTS_DIR -M msa -os
where RESULTS_DIR is the directory containing your orthogroups results files (e.g. Orthogroups.csv)
```
orthofinder -fg /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.larvae.fastas/Results_Jan17_1 -M msa -os
```
Combine individual fasta files into one fasta file (located in /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.adult.fastas/Results_Jan17_1/Orthologues_Jan21/Sequences)
```
cat *.fa > da.larvae.ortho.fasta
```
Move fasta file to main directory
```
cp da.larvae.ortho.fasta /home/craker/diadema/
```

##### Filtering out only overlapping orthogroups from original fasta file
In this analysis, we want to focus on orthogroups found in common between *D. antillarum* and *S. purpuratus*


Navigate to .csv files containing overlapping orthogroups and examine file
```
cd /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.larvae.fastas/Results_Jan17_1/Orthologues_Jan17/Orthologues/Orthologues_da.larvae.best
nano da.larvae.best__v__SPU_peptide.csv
```
headers from original Trinity fasta file are all located in the second column
Use cut to create a text file with only this column
```
cut -f2 da.larvae.best__v__SPU_peptide.csv > da.larvae.onecolumn.txt
```
Copy onecolumn.txt to main directory for ease of further analysis and navigate back up to that directory
```
cp da.larvae.onecolumn.txt /home/craker/diadema/
cd /home/craker/diadema/
```

(reactivate conda environment)


###### Use bbtools script filterbyname.sh and da.larvae.onecolumn.txt to create a new fasta file from da.larvae.fasta containing only the desired transcripts
Use include=t so output includes instead of excludes files named in onecolumn.txt
Use substring=t since Orthofinder added .p1.p1 to the end of headers
```
conda install bbtools
filterbyname.sh in=da.larvae.fasta names=da.larvae.onecolumn.txt out=da.larvae.filtered.fasta include=t substring=t
```
Use the da.larvae.filtered.fasta file for further downstream analysis

##### Use Trinotate to annotate the filtered fasta file
Install other necessary programs
```
conda install hmmer
conda install perl-dbd-sqlite
```
Run TransDecoder on the filtered fasta file
```
TransDecoder.LongOrfs -t da.larvae.filtered.fasta
TransDecoder.Predict -t da.larvae.filtered.fasta
```

Blast the peptide file
```
blastp -query /home/craker/diadema/da.larvae.filtered.fasta.transdecoder.pep -db /home/craker/diadema/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > da.larvae.filtered.blastp.outfmt6
```
Run HMMER
```
hmmscan --cpu 20 --domtblout TrinotatePFAM.out Pfam-A.hmm da.larvae.filtered.fasta.transdecoder.pep > da.larvae.filtered.pfam.log
```

Load files into an SQLite database
```
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init --gene_trans_map da.larvae.filtered.fasta.gene_trans_map --transcript_fasta da.larvae.filtered.fasta --transdecoder_pep da.larvae.filtered.fasta.transdecoder.pep
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp da.larvae.filtered.blastp.outfmt6
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
```

# (didn't do because did not generate a blastx file for some reason)
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.filtered.outfmt6


Generate a gene/transcript relationship file using TRINITY
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl da.larvae.filtered.fasta > da.larvae.filtered.fasta.gene_trans_map
```
Generate a Trinotate annotation report
```
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > da.larvae.filtered.trinotate_annotation_report.xls
```

##### Trinity transcript quantification
Necessary to perform further downstream analysis.

Install necessary programs
```
conda install salmon
conda install rsem
conda install bowtie2
```
Created a tab delimited text file of samples (da_samples_adult.txt)


###### Running the analysis with the original fasta file:
Prepare the reference and run alignment and abundance estimation
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl \
--transcripts da.larvae.fasta --est_method salmon --aln_method bowtie2 --trinity_mode \
--prep_reference --samples_file da_samples_larvae.txt --seqType fq
```


Create matrices
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl \
--est_method salmon --gene_trans_map none --out_prefix da.larvae.salmon --name_sample_by_basedir \
da.larvae.fasta.pH/pH_high_rep1/quant.sf da.larvae.fasta.pH/pH_high_rep2/quant.sf da.larvae.fasta.pH/pH_high_rep3/quant.sf \
da.larvae.fasta.pH/pH_medx_rep1/quant.sf da.larvae.fasta.pH/pH_medx_rep2/quant.sf da.larvae.fasta.pH/pH_medx_rep3/quant.sf \
da.larvae.fasta.pH/pH_lowx_rep1/quant.sf da.larvae.fasta.pH/pH_lowx_rep2/quant.sf da.larvae.fasta.pH/pH_lowx_rep3/quant.sf \
da.larvae.fasta.pH/pH_lowx_rep4/quant.sf
```


These matrix files, specifically the salmon.isoforms.count.matrix file, will be used for gene expression analysis.


Filter the results to only include most differentially espressed isoforms
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/filter_low_expr_transcripts.pl \
--m da.larvae.salmon.isoform.TPM.not_cross_norm --t da.larvae.fasta \
--min_expr_any 0 --highest_iso_only --trinity_mode \
> da.larvae.filter_hio.fasta
```
Retained 235804 / 389705 = 60.51% of total transcripts.

```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/filter_low_expr_transcripts.pl \
--m da.larvae.salmon.isoform.TPM.not_cross_norm --t da.larvae.fasta \
--min_expr_any 1 --min_pct_dom_iso 1 --trinity_mode \
> da.larvae.filter_expr1.fasta
```
Retained 226005 / 389705 = 57.99% of total transcripts.


Can create matrices using these results as well, just re-run above scripts with da.larvae.filter_hio.fasta and da.larvae.filter_expr1.fasta in place of da.larvae.fasta


###### Create matrix from da.filtered.fasta
Prepare the reference and align
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl \
--transcripts da.larvae.filtered.fasta --est_method salmon --aln_method bowtie2 \
--prep_reference --trinity_mode --samples_file da_samples_larvae.txt --seqType fq
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
###### Using my files generated from original da.adult.fasta:
 ```
 /home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/run_DE_analysis.pl \
 --matrix da.larvae.salmon_out_dir/da.larvae.salmon.isoform.counts.matrix \
 --method DESeq2 --samples_file da_samples_larvae.txt
 ```
Results of this script:
```
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.count_matrix
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.DE_results
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.DE_results.MA_n_Volcano.pdf
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.Rscript
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.count_matrix
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.DE_results
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.DE_results.MA_n_Volcano.pdf
da.larvae.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.Rscript
da.larvae.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.count_matrix
da.larvae.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.DE_results
da.larvae.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.DE_results.MA_n_Volcano.pdf
da.larvae.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.Rscript
```

To extract and cluster differentially expressed transcripts, I would use this script:
```
$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl
```
Copy sample text file into the DESeq2_outdir directory, move there, and run script like this (MUST run from the DESeq2 directory!):
```
cp da_samples_larvae.txt /home/craker/diadema/da.larvae.DESeq2_outdir
cd da.larvae.DESeq2_outdir

/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix /home/craker/diadema/da.larvae.salmon_out_dir/da.larvae.salmon.isoform.TMM.EXPR.matrix --samples da_samples_larvae.txt

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












.
