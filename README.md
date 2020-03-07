# Redside Dace RNA project
A repository for scripts assembling and analyzing RNA seq data for redside dace. Brief folder descripts are below, and more detailed information on scripts are within the readme files of each folder.

- `shell_scripts` contains code for taking raw mRNA reads and:
  - Creating a *de novo* transcriptome with Trinity
  - Transcript and superTranscript quantification using Salmon and Corset
  
 - `transcriptome_annotation` contains code for taking the Trinity-assembled transcriptome and finding known transcript function, following the Trinotate pipeline. 
