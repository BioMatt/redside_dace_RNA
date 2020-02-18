# Shell scripts

These scripts were run on the [Cedar cluster](https://docs.computecanada.ca/wiki/Cedar) on [Westgrid](https://www.westgrid.ca/) and [Compute Canada](https://www.computecanada.ca/), unless otherwise stated below. 

The pipeline here follows the [superTranscripts](https://github.com/Oshlack/Lace/wiki/Example:-Differential-Transcript-Usage-on-a-non-model-organism) pipeline pretty closely. 


#### The scripts were used in this order:
  - `raw_fastqc.sh` used to take a look at quality before trimming.
    - the `raw_dace_reads.txt` file used to list the unpaired reads and file locations for trimming was made by navigating to the folder they were in and using `printf '%s\n' "$PWD"/* > raw_dace_reads.txt`.
    - MultiQC was used after this to collate fastqc reports into a single report, using Ubuntu on a Windows PC. 
  - `trim_raw_reads.sh` is where trimmomatic was called, using the `TruSeq3-SE.fa` file for adapters.
  - `trimmed_fastqc.sh` is trimmomatic as well, and the `trimmed_read_list.txt` file was made using the `printf` command above.
    - As before, MultiQC was used to collate quality reports, again on Ubuntu and not the cluster. 
  - `trinity_assembly_skylake.sh` is where Trinity is used to assemble a _de novo_ transcriptome. The job took ~11.5 days to run.
    - `trinity_sample_list.txt` includes the file path locations, copied from the trimmed read list text file and a little bit of metadata on each sample, with 3 groups in this experiment.
  - `salmon_index_pre_corset.sh` is where an index of the trascriptome is generated for using Salmon to quantify transcript abundance. 
  - `salmon_quant_pre_corset.sh` is where Salmon is actually used to quantify transcript abundance. Equivalence classes (`--dumpEq`) is used with Salmon because Corset requires that information for clustering the transcripts.
    - the `array_paired_trimmed_reads.txt` file was made with `printf` as in the `raw_dace_reads.txt` file, but with each forward and reverse file on the same line. Salmon quantification is then run in an array, so 30 jobs (one for each sample) are run at once on the cluster. Each job represents one sample, or one line from the array text file. 
     - the directory format for outputting data when Salmon finishes led to many extra direcories the way the script is written here.
  - `corset_dace.sh` takes the equivalence classes text file from each of the salmon outputs, and clusters transcripts based on Salmon's quasi alignments. 
