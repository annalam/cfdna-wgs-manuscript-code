# cfdna-wgs-manuscript-code
This repository contains the code used to carry out the analyses in our article "Somatic chronology of treatment-resistant prostate cancer via deep whole-genome ctDNA sequencing".




## Pre-requisites

Our analysis software and scripts are written using the Julia and Rust programming languages. To use our code, you will first need to:

1. Install Julia: https://julialang.org/downloads/
2. Install the Rust compiler and toolchain: https://www.rust-lang.org/tools/install
3. Install required Julia packages:
```
import Pkg;
Pkg.add(split("Distributions HypothesisTests PyCall KernelDensity Interpolations"))
```
4. Modify your JULIA_LOAD_PATH environment variable so that Julia can find the code modules found in the root directory of this repository:
```
export JULIA_LOAD_PATH=/path/to/this/repository:
```
5. Add the Julia scripts found in subfolder `scripts` of this repository to your PATH environment variable.
6. Install Mutato: https://www.github.com/annalam/mutato
7. Install Breakfast: https://www.github.com/annalam/breakfast
8. Install ANNOVAR: https://annovar.openbioinformatics.org/en/latest/user-guide/download/

Steps 5 - 8 are only required for some of the analysis steps, so you may be able to skip some steps depending on what analyses you wish to run.

Simulated patient data files for experimenting with our subclonal reconstruction algorithms is also available here:
https://www.dropbox.com/s/hofrv1vlnln1mau/example_patients.zip?dl=0

Patient whole-genome sequencing data have been deposited in the European Genome-Phenome Archive (EGA) database under the accession code EGAS00001005783 and is available under standard EGA controlled release.


## Read alignment

All plasma cfDNA and metastatic tissue biopsy samples were sequenced using Illumina X Ten instruments. The paired end sequencing reads were adapter-trimmed using cutadapt-1.11, quality-masked using seqkit-0.8, and then aligned against the human hg38 reference genome using Bowtie-2.3.0. Duplicate DNA fragments were identified and marked using samblaster-0.1.24:
```
fasta interleave sample_1.fq.gz sample_2.fq.gz | \
  cutadapt --interleaved -f fastq -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC - | \
  fasta trim by quality - 30 | fasta mask by quality - 20 | \
  bowtie2 -p20 -X 1000 --score-min L,0,-0.6 --ignore-quals -x hg38 --interleaved - | \
  samblaster | samtools view -u - | samtools sort -@ 8 -m 4G -o sample.bam
```





## Somatic mutation analysis

As the first step in mutation analysis, we generated a matrix containing read-level evidence for each potential variant (row) in each sample (column):
```
mutato call --alt-reads=5 --alt-frac=0.05 --min-mapq=0 hg38.fa *.bam > variants.tsv
```

Next we generated a file `tumor_normal_pairs.txt` listing all tumor and cfDNA samples in the cohort, together with their matched germline samples. Here is an example:
```
TEST	REF
AE-015-Baseline-cfDNA	AE-015-WBC
AE-015-Progression2-cfDNA		AE-015-WBC
AE-018-Baseline-cfDNA	AE-018-WBC
```

Cancer samples are listed in the first column, and matched germline samples in the second column. Additional negative control cfDNA samples were included by listing each of them on their own line, with an empty first column.

We searched for somatic mutations fulfilling the following criteria:
- At least 8 mutant allele reads
- Mutant allele fraction ≥ 10%
- Mutant allele fraction 10x higher than in the matched germline sample
- Mutant allele fraction 50x higher than the average of all cancer-negative samples (germline samples and known cancer-negative cfDNA samples)
- A minimum read depth of 10x in the matched germline sample
- Average mutant allele distance from nearest read end ≥ 15
- Average mapping quality of mutant allele reads ≥ 20

This was done using the following command:
```
variant somatic --alt-reads=8 --alt-frac=0.1 --test-ref-ratio=10 \
  --test-bg-ratio=50 --ref-reads=10 --min-sidedness=15 --min-mapq=20 \
  variants.tsv tumor_normal_pairs.txt > somatic.tsv
```

Next, we annotated all mutations with their predicted biological effect using `variant predict effect` (which internally uses ANNOVAR). We also annotated each mutation with information about its frequency in the COSMIC (version 77) and GNOMAD (version 3.0) databases:
```
variant predict effect somatic.tsv | \
  variant annotate - cosmic_77_hg38.jls | variant annotate gnomad_v3.jls | \
  variant discard if frequency above gnomad_v3.jls 0.005 > somatic_annotated.tsv
```



## Germline heterozygous SNP identification

We generated a list of germline heterozygous SNPs for each patient, using the following command:
```
variant heterozygous snps --min-depth=30 --max-depth=500 --min-mapq=30 \
  --min-sidedness=15 variants.tsv WBC | variant discard indels - | \
  variant annotate - gnomad_v3.jls | egrep 'CHROM|GNOMAD' > hetz_snps.tsv
```

Indels were omitted since they are associated with alignment artifacts that can compromise accurate allele fraction quantification.


## Deleterious germline variant identification

To search for germline variants, we looked for variants found in the matched leukocyte samples that fulfilled the following criteria:
- At least 5 reads supporting the variant
- Allele fraction ≥ 20%
- Allele fraction is at least 20 times higher than the background error rate (i.e. the average allele fraction of leukocyte samples that had an allele fraction < 20% for the variant)

This analysis was carried out using the `variant germline` tool:
```
variant germline --alt-reads=5 --alt-frac=0.2 --bg-ratio=20 \
  variants.tsv WBC > germline.tsv
```

We further narrowed this list down to protein altering germline variants with a population frequency below 0.5%, and annotated the variants with information from the COSMIC and ClinVar (dated 2020-07-06)  databases:
```
variant predict effect germline.tsv | variant protein altering - | \
  variant discard if frequency above - gnomad_v3.jls 0.005 | \
  variant annotate - gnomad_v3.jls | \
  variant annotate - cosmic_77_hg38.jls | \
  variant annotate - clinvar-2020-07-06_hg38.jls \
  > germline_annotated.tsv
```

The resulting list of rare germline variants was curated for deleterious alterations by looking for protein-truncating mutations in DNA repair genes and other critical cancer genes, and by searching for variants annotated as “pathogenic” in ClinVar.


## Copy number analysis

First we generated a BED file describing a grid of half-overlapping 1000 bp windows covering the entire genome:
```
copynum grid hg38.chrom.sizes 1000 > grid.bed
```

Then we calculated the GC nucleotide content of each window, for use in subsequent GC bias correction:
```
fasta gc content hg38.fa grid.bed > grid.gc
```

Next, we counted the number of concordantly aligned DNA fragments within each 1000 bp genomic window using the `sam count` tool:
```
sam count sample.bam grid.bed > sample.tsv
```

Finally we generated coverage logratio and heterozygous SNP allele fraction (HSAF) tracks using the `copynum call genomewide` tool:
```
copynum call genomewide --output-format=igv --logratio-dots=Inf \
  --snp-median-decimate=9 --discard-noisiest=0.05 --min-ref-depth=30 \ 
  --controls=cna_controls.txt --gc-fractions=grid.gc --report-dir=./ \
  --hetz-snps=hetz_snps.tsv grid.bed tumor_normal_pairs.txt
```

We used the same `tumor_normal_pairs.txt` file used in the mutation analysis step.

## Somatic chromosomal rearrangement analysis

We first generated a list of candidate rearrangements using `breakfast detect`, which carries out a split-read analysis to identify DNA sequencing reads that overlap a rearrangement junction:
```
breakfast detect --max-frag-len=1000 --anchor-len=30 sample.bam hg38 \
  > rearrangements/original/sample.sv
```

To remove germline rearrangements and false positives, we filtered out any rearrangements that were also found in the patient’s matched germline sample, or in any sample from another patient. We also required a minimum of 5 unique supporting reads for somatic rearrangements:
```
breakfast blacklist matched_leukocyte.sv other_patient_sample_1.sv … other_patient_sample_N.sv > blacklist.txt
breakfast filter --min-reads=5 --merge-duplicates --blacklist=blacklist.txt sample.sv
```

Finally we annotated all rearrangements with information regarding the location of their genomic breakpoints relative to known genes:
```
breakfast annotate sample.filtered.sv ensembl_84_genes.bed > sample.annotated.sv
```




## Inference of cancer fraction and truncal copy number profile

To infer a truncal copy number model for each patient sample, we started out with files `*_logratio.igv` and `*_hetz_snp.tsv` for each sample, and `*_mutations.tsv` for the patient. These files were created using the procedures described in earlier sections of this document. We have made [example data](https://www.dropbox.com/s/hofrv1vlnln1mau/example_patients.zip?dl=0) available for two simulated patient cases PatientA and PatientB (see section “Generation of simulated cancer patients”).

The first step in truncal copy number model was to identify the boundaries of genomic copy number alterations using a joint analysis of all same-patient samples:
```
copynum segment PatientA
```

This produces a `*_segments.tsv` file for each patient sample. We then used the segmented copy number profiles to infer the cancer fraction and truncal copy number profile of all patient samples:
```
subclones truncal model PatientA
```

The command creates PDF figures `*_truncal_model.pdf` showing the inferred truncal copy number model for each sample. Additionally, a file called `PatientA_populations.tsv` is generated, listing the cancer fraction and diploid level of each sample in a tabular format.


## Subclonal reconstruction

After identifying the truncal copy number model for each patient sample, we carried out somatic mutation clustering, subclone identification, and phylogenetic tree inference:
```
subclones cluster PatientA
```

This updates the file `PatientA_populations.tsv` with additional rows and columns listing the different cancer cell populations identified in the patient’s samples, as well as their subclonal fractions in each sample. A PDF figure `PatientA_subclones.pdf` is also generated, visualizing the mutation clusters.

After the subclone fractions and phylogenetic tree had been determined, we reconstructed the copy number profiles of individual subclones:
```
subclones copy number deconvolution PatientA
```

This analysis creates a file `PatientA_cn_deconvolution.pdf` showing the subclone copy number profiles.


## Generation of simulated cancer patients

To generate simulated cancer patient `PatientC` with two samples and three cancer cell subpopulations (subclones), we used our script ‘simulate wgs`:
```
simulate wgs PatientD 2
```


## GATK Mutect2 validation

We used GATK Mutect2 version 4.2.3.0 to validate our somatic mutation calls. The analysis was carried out according to the [recommended workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2) described in the GATK Mutect2 documentation.

First we aligned our paired end Illumina sequencing reads against the hg38 human reference genome using BWA version 0.7.17:
```
fasta interleave sample_1.fq.gz sample_2.fq.gz | \
  cutadapt --interleaved -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC - | \
  bwa mem -p -t 20 \
    -R "@RG\tSM:sample\tID:sample.L1\tLB:0\tPL:ILLUMINA\tPU:sample_1" hg38.fa - | \  
  samblaster | samtools sort -m 4G -o sample.bam -
```

We then carried out base quality recalibration on the aligned BAM files:
```
gatk BaseRecalibrator -I sample.bam -R hg38.fa \
  --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf \
  --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf \
  --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
  -O sample.bqsr_table
gatk ApplyBQSR --bqsr-recal-file sample.bqsr_table --input sample.bam \
  --output sample.bqrecalib.bam
```

Next we ran GATK Mutect2 on matched leukocyte and ctDNA-negative cfDNA samples to generate initial variant allele pileups. The following command was run on each control sample’s base quality recalibrated BAM file:
```
gatk Mutect2 -R hg38.fa --max-mnp-distance 0 \
  --intervals broad_wgs_calling_regions.hg38.interval_list \
  -I normal_sample.bqrecalib.bam -O normal_sample.vcf.gz
```

After generating a VCF file for each control sample, we combined them into a single panel-of-normals database:
```
gatk GenomicsDBImport -R hg38.fa \
  -L broad_wgs_calling_regions.hg38.interval_list \
  --genomicsdb-workspace-path pon_db \
  -V normal_sample_1.vcf.gz  ...  -V normal_sample_N.vcf.gz
```

Finally, we generated an initial list of somatic mutations for each tumor and cfDNA cancer sample using the following command:
```
gatk Mutect2 -R hg38.fa -L broad_wgs_calling_regions.hg38.interval_list \
  -I cancer_sample_1.bqrecalib.bam \
  -I cancer_sample_2.bqrecalib.bam \
  -I cancer_sample_3.bqrecalib.bam \
  -l matched_normal.bqrecalib.bam \
  -normal matched_normal -pon pon.vcf.gz \
  --f1r2-tar-gz patient.f1r2.tar.gz \
  -O patient.somatic.vcf.gz
```

The optional argument `--f1r2-tar-gz patient.f1r2.tar.gz` instructed GATK Mutect2 to collect statistics into file `patient.f1r2.tar.gz` that was then used to generate an orientation bias model:
```
gatk LearnReadOrientationModel -I patient.f1r2.tar.gz -O patient.orientation_model.tar.gz
```

Finally, we filtered the somatic mutation candidates using GATK Mutect2’s post-processing filters, including the orientation bias filter:
```
gatk IndexFeatureFile -I patient.somatic.vcf.gz
gatk FilterMutectCalls -R hg38.fa \
  --ob-priors patient.orientation_model.tar.gz \
  -V patient.somatic.vcf.gz -O patient.somatic_final.vcf.gz
```



## PyClone validation of subclonal reconstruction

As part of our effort to validate our subclonal reconstructions, we ran the third-party software PyClone on our mutation data. First, we generated sample-specific input files in PyClone’s input format which includes an allele-specific copy number for the genomic location of each mutation. Mutations located in genomic regions that displayed subclonal copy number heterogeneity in one or more patient samples were excluded from the input file.

We then ran PyClone on each patient’s sample input files using this command: 
```
PyClone run_analysis_pipeline --in_files sample_1.tsv sample_2.tsv \
  --working_dir ./ --tumour_contents 0.45 0.78
```

After PyClone had inferred per-sample CCF values for all somatic mutations provided as input, we read them from the output file `tables/loci.tsv`, and visualized them as a scatter plot.





## PhyloWGS validation

As part of our effort to validate our subclonal reconstructions, we ran the third-party software PhyloWGS on our mutation data. First, we converted our somatic mutation tables into input files in the `mutect_tcga` format supported by PhyloWGS. Then we converted our copy number segmentation files into the input files format required by PhyloWGS. Once this was done, we used the following PhyloWGS script to pre-process the data and merge the sample-specific files:
```
create_phylowgs_inputs.py --output-variants ssm.txt \
  --output-cnvs cnv.txt \
  --vcf-type AE-015-Baseline-cfDNA=mutect_tcga \
  --vcf-type AE-015-Progression2-cfDNA=mutect_tcga \
  --cnvs AE-015-Baseline-cfDNA=AE-015-Baseline-cfDNA.seg \
  --cnvs AE-015-Progression2-cfDNA=AE-015-Progression2-cfDNA.seg \
  AE-015-Baseline-cfDNA=AE-015-Baseline-cfDNA.vcf \ 
  AE-015-Progression2-cfDNA=AE-015-Progression2-cfDNA.vcf
```

The resulting `ssm.txt` and `cnv.txt` files were used as inputs for mutation clustering and subclone identification:
```
multievolve.py --num-chains 4 --ssms ssm.txt --cnvs cnv.txt --output-dir chains
```

After the analysis had completed, we ran the PhyloWGS `write_results.py` script that generates the final output files:
```
write_results.py AE-015 trees.zip AE-015.summ.json.gz \
  AE-015.muts.json.gz AE-015.mutass.zip
```

PhyloWGS does not report CCF values for individual mutations, but instead only reports mutation-cluster assignments and cluster CCFs. Furthermore, PhyloWGS actually reports thousands of possible subclonal reconstructions, which are scored based on a normalized likelihood value assigned to each reconstruction. To compare our subclonal reconstructions against the best PhyloWGS reconstruction, we used the normalized likelihood calculation formula in PhyloWGS tree_viewer.js to identify the numeric identifier of the best reconstruction (e.g. `1452`), using information available in the `AE-015.summ.json` file. We then read mutation-to-cluster assignments from the file `1452.json` found inside the `AE-015.mutass.zip` file.

Since PhyloWGS does not output any per-mutation CCF values, we could not directly generate mutation CCF scatterplots that would have allowed visual investigation of the validity of the cluster assignments. So we instead generated mutation CCF scatterplots using the mutation CCF values inferred by our pipeline, and colored the mutations in these scatterplots based on their PhyloWGS cluster assignments.






