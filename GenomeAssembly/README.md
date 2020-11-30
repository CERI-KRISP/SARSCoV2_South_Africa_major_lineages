### Pre-processing Raw Reads from Illumina

Raw reads coming from Illumina sequencing were assembled using Genome Detective 1.126 [Genome Detective](https://www.genomedetective.com/) and the Coronavirus Typing Tool 7,8. 

The initial assembly obtained from Genome Detective was polished by aligning mapped reads to the references and filtering out low-quality mutations using bcftools 1.7-2 mpileup method. 

The bash script for this post-processing step is available here and can be used as follows. The files needed for this command can be download on the report page of the Genome Detective sequence assembly. 

`./makeconsensus-illumina2.sh <reference_sequence> <BAM sequence alignment> <output fasta file>` 
