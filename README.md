This repository contains code used to generate the results in study "Somatic chronology of treatment-resistant prostate cancer via deep whole-genome ctDNA sequencing", currently in review.

All analyses were run on a computational server running the Linux operating system. No non-standard hardware was utilized.

Sequenced paired end reads were aligned against the hg38 human reference genome using Bowtie 2.3.0:  
https://github.com/BenLangmead/bowtie2

Somatic and germline mutation analysis was carried out using the in-house software "mutato". The Rust source code for this software can be found under subfolder "mutato" of this repository.

Rearrangement analysis was carried out using the Breakfast software:  
https://github.com/annalam/breakfast

Read counting within genomic windows was carried out using BEDTOOLS:  
https://github.com/arq5x/bedtools2

FASTQ preprocessing and quality metrics calculation were performed using Seqkit:  
https://github.com/annalam/seqkit
