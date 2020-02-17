## Gene Expression Analysis with only adult urchins
#### Last updated: February 16, 2020

run Trinity on all clean files, using --full_cleanup option to delete uneeded files at the end of the analysis to save space
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Trinity \
--seqType fq --max_memory 50G \
--single OADA0006_S16_L002_R1_001_clean.fq,\
OADA0006_S16_L002_R1_001_clean.fq,\
OADA0049_S25_L003_R1_001_clean.fq,\
OADA0058_S19_L002_R1_001_clean.fq,\
OADA0071_S53_L005_R1_001_clean.fq,\
OADA0081_S67_L006_R1_001_clean.fq,\
OADA0085_S24_L002_R1_0011_clean.fq,\
OADA0102_S33_L003_R1_001_clean.fq,\
OADA0116_S26_L003_R1_001_clean.fq,\
OADA0139_S17_L002_R1_001_clean.fq,\
OADA0174_S20_L002_R1_001_clean.fq \
--CPU 20 --full_cleanup
```

is running but maybe has errors:
```
^CError, cmd: seqtk-trinity seq -A /home/craker/diadema/OADA0006_S16_L002_R1_001_clean.fq >> single.fa died with ret 2 at /home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/insilico_read_normalization.pl line 762.
Error, cmd: /home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/insilico_read_normalization.pl --seqType fq --JM 50G  --max_cov 200 --min_cov 1 --CPU 20 --output /home/craker/diadema/trinity_out_dir/insilico_read_normalization   --max_pct_stdev 10000  --single /home/craker/diadema/OADA0006_S16_L002_R1_001_clean.fq,/home/craker/diadema/OADA0006_S16_L002_R1_001_clean.fq,/home/craker/diadema/OADA0049_S25_L003_R1_001_clean.fq,/home/craker/diadema/OADA0058_S19_L002_R1_001_clean.fq,/home/craker/diadema/OADA0071_S53_L005_R1_001_clean.fq,/home/craker/diadema/OADA0081_S67_L006_R1_001_clean.fq,/home/craker/diadema/OADA0085_S24_L002_R1_0011_clean.fq,/home/craker/diadema/OADA0102_S33_L003_R1_001_clean.fq,/home/craker/diadema/OADA0116_S26_L003_R1_001_clean.fq,/home/craker/diadema/OADA0139_S17_L002_R1_001_clean.fq,/home/craker/diadema/OADA0174_S20_L002_R1_001_clean.fq died with ret 512 at /home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Trinity line 2689.
```


rename fasta file
```
mv trinity_out_dir.Trinity.fasta da.adult.fasta
```

calculate number of clusters
```
cat da.adult.fasta | grep '>' | wc -l
```
436218 clusters

#### TransDecoder and CD-hit
```
TransDecoder.LongOrfs -t da.adult.fasta
TransDecoder.Predict -t da.adult.fasta
cd-hit-est -i ./da.adult.fasta.transdecoder_dir/longest_orfs.cds -o cdhit89 -n 8 -c 0.89 -M 96000 -T 20 -d 0
```

calculate number of clusters
```
cat cdhit89 | grep '>' | wc -l
```
53362 clusters

run TransDecoder on cd-hit output
```
TransDecoder.LongOrfs -t cdhit89
TransDecoder.Predict -t cdhit89
```

rename best fasta file (so it's easy to find later)
```
mv cdhit89.transdecoder.pep da.adult.best.fasta
```

#### Trinotate

BLAST
```
blastx -query /home/craker/diadema/da.adult.fasta -db /home/craker/diadema/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
blastp -query /home/craker/diadema/da.adult.transdecoderCDhit_dir/cdhit89.transdecoder.pep -db /home/craker/diadema/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > da.adult.blastp.outfmt6
```


#### Orthofinder

Copy fasta files to OrthoFinder directory
```
cd OrthoFinder-2.3.1_source/orthofinder/
mkdir da.adult.fastas
cp /home/craker/diadema/da.adult.transdecoderCDhit_dir/da.adult.best.fasta /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.adult.fastas
cp /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/DAfastas/SPU_peptide.fasta /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.adult.fastas
```

run OrthoFinder (remember to deactivate conda environment)
```
orthofinder -f ./da.adult.fastas/ -t 20
```

Generate a fasta file for each orthogroup by running:
orthofinder -fg RESULTS_DIR -M msa -os
where RESULTS_DIR is the directory containing your orthogroups results files (e.g. Orthogroups.csv)
```
orthofinder -fg /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.adult.fastas/Results_Jan17_1 -M msa -os
```
Combine individual fasta files into one fasta file (located in /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.adult.fastas/Results_Jan17_1/Orthologues_Jan21/Sequences)
```
cat *.fa > da.adult.ortho.fasta
```
Move fasta file to main directory
```
cp da.adult.ortho.fasta /home/craker/diadema/
```

(reactivate conda environment)

##### Filtering out only overlapping orthogroups from original fasta file
In this analysis, we want to focus on orthogroups found in common between *D. antillarum* and *S. purpuratus*


Navigate to .csv files containing overlapping orthogroups and examine file
```
cd /home/craker/diadema/OrthoFinder-2.3.1_source/orthofinder/da.adult.fastas/Results_Jan17_1/Orthologues_Jan17/Orthologues/Orthologues_da.adult.best
nano da.adult.best__v__SPU_peptide.csv
```
headers from original Trinity fasta file are all located in the second column
Use cut to create a text file with only this column
```
cut -f2 da.adult.best__v__SPU_peptide.csv > da.adult.onecolumn.txt
```
Copy onecolumn.txt to main directory for ease of further analysis and navigate back up to that directory
```
cp da.adult.onecolumn.txt /home/craker/diadema/
cd /home/craker/diadema/
```
###### Use bbtools script filterbyname.sh and onecolumn.txt to create a new fasta file from da.adult.fasta containing only the desired transcripts
Use include=t so output includes instead of excludes files named in onecolumn.txt
Use substring=t since Orthofinder added .p1.p1 to the end of headers
```
conda install bbtools
filterbyname.sh in=da.adult.fasta names=da.adult.onecolumn.txt out=da.adult.filtered.fasta include=t substring=t
```
Use the da.adult.filtered.fasta file for further downstream analysis


##### Use Trinotate to annotate the filtered fasta file
Install other necessary programs
```
conda install hmmer
conda install perl-dbd-sqlite
```
Run TransDecoder on the filtered fasta file
```
TransDecoder.LongOrfs -t da.adult.filtered.fasta
TransDecoder.Predict -t da.adult.filtered.fasta
```
Blast the peptide file
```
blastp -query /home/craker/diadema/da.adult.filtered.fasta.transdecoder.pep -db /home/craker/diadema/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > da.adult.filtered.blastp.outfmt6
```
Run HMMER
```
hmmscan --cpu 20 --domtblout TrinotatePFAM.out Pfam-A.hmm da.adult.filtered.fasta.transdecoder.pep > da.adult.filtered.pfam.log
```

Load files into an SQLite database
```
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init --gene_trans_map da.adult.fasta.gene_trans_map --transcript_fasta da.adult.filtered.fasta --transdecoder_pep da.adult.filtered.fasta.transdecoder.pep
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp da.adult.filtered.blastp.outfmt6
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
```

(didn't do because did not generate a blastx file for some reason)
```
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.filtered.outfmt6
```


Generate a gene/transcript relationship file using TRINITY
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl da.adult.filtered.fasta > da.adult.filtered.fasta.gene_trans_map
```
Generate a Trinotate annotation report
```
/home/craker/diadema/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > da.adult.filtered.trinotate_annotation_report.xls
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
--transcripts da.adult.fasta --est_method salmon --aln_method bowtie2 --trinity_mode \
--prep_reference --samples_file da_samples_adult.txt --seqType fq
```
have to use original fasta for some reason


Create matrices
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl \
--est_method salmon --gene_trans_map none --out_prefix da.adult.salmon --name_sample_by_basedir \
da.adult.fasta.pH/pH_high_rep4/quant.sf da.adult.fasta.pH/pH_high_rep5/quant.sf da.adult.fasta.pH/pH_high_rep6/quant.sf da.adult.fasta.pH/pH_high_rep7/quant.sf \
da.adult.fasta.pH/pH_medx_rep4/quant.sf da.adult.fasta.pH/pH_medx_rep5/quant.sf da.adult.fasta.pH/pH_medx_rep6/quant.sf \
da.adult.fasta.pH/pH_lowx_rep5/quant.sf da.adult.fasta.pH/pH_lowx_rep6/quant.sf da.adult.fasta.pH/pH_lowx_rep7/quant.sf da.adult.fasta.pH/pH_lowx_rep8/quant.sf
```


These matrix files, specifically the salmon.isoforms.count.matrix file, will be used for gene expression analysis.


Filter the results to only include most differentially espressed isoforms
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/filter_low_expr_transcripts.pl \
--m da.adult.salmon.isoform.TPM.not_cross_norm --t da.adult.fasta \
--min_expr_any 0 --highest_iso_only --trinity_mode \
> da.adult.filter_hio.fasta
```
Retained 267768 / 436218 = 61.38% of total transcripts.

```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/filter_low_expr_transcripts.pl \
--m da.adult.salmon.isoform.TPM.not_cross_norm --t da.adult.fasta \
--min_expr_any 1 --min_pct_dom_iso 1 --trinity_mode \
> da.adult.filter_expr1.fasta
```
Retained 295569 / 436218 = 67.76% of total transcripts.


Can create matrices using these results as well, just re-run above scripts with da.adult.filter_hio.fasta and da.adult.filter_expr1.fasta in place of da.adult.fasta


###### Create matrix from da.filtered.fasta
Prepare the reference and align
```
/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl \
--transcripts da.adult.filtered.fasta --est_method salmon --aln_method bowtie2 \
--prep_reference --trinity_mode --samples_file da_samples_adult.txt --seqType fq
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
 ??????????
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
 --matrix da.adult.salmon_out_dir/da.adult.salmon.isoform.counts.matrix \
 --method DESeq2 --samples_file da_samples_adult.txt
 ```
Results of this script:
```
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.count_matrix
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.DE_results
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.DE_results.MA_n_Volcano.pdf
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_lowx.DESeq2.Rscript
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.count_matrix
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.DE_results
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.DE_results.MA_n_Volcano.pdf
da.adult.salmon.isoform.counts.matrix.pH_high_vs_pH_medx.DESeq2.Rscript
da.adult.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.count_matrix
da.adult.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.DE_results
da.adult.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.DE_results.MA_n_Volcano.pdf
da.adult.salmon.isoform.counts.matrix.pH_lowx_vs_pH_medx.DESeq2.Rscript
```

To extract and cluster differentially expressed transcripts, I would use this script:
```
$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl
```
Copy sample text file into the DESeq2_outdir directory, move there, and run script like this:
```
cp da_samples_adult.txt /home/craker/diadema/da.adult.DESeq2_outdir
cd da.adult.DESeq2_outdir

/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix /home/craker/diadema/da.adult.salmon_out_dir/da.adult.salmon.isoform.TMM.EXPR.matrix --samples da_samples_adult.txt

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
