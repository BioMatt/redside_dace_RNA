# Transcriptome annotation steps

These lines of code and scripts were run with the [Trinity-assembled transcriptome](https://github.com/BioMatt/redside_dace_RNA/blob/master/shell_scripts/trinity_assembly_skylake.sh) for redside dace to annotate the transcriptome with gene/transcript functions. The code follows the [Trinotate pipeline](https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required), except RNAmmer is *not* run because it was too much trouble to get working. 

Unless otherwise stated (i.e., everything except the tmhmm and signalP scripts) code was run using [Bioconda](https://bioconda.github.io/user/install.html), installed on Ubuntu, on a 12-core desktop PC running Windows 10. 


#### Code was used in this order:
  - TransDecoder was used first to identify the longest open reading frames in transcripts, and output `longest_orfs.pep`, which identified protein sequences within the transcriptome. `conda install transdecoder` was used to install the program.  
    - `TransDecoder.LongOrfs -t dace_transcriptome/dace_transcriptome_Trinity.fasta -O transdecoder_out`  
  - Trinotate was then installed (`conda install trinotate`), and the three databases needed going forward were built. There is probably a *much* cleaner way to do this, but I went into the trinotate binaries (bin) directory and built the databases there.   
    - Building the database for trinotate to pull together annotation data at the end:  
      - `Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate`  
      - This makes the file `Trinotate.sqlite`  
    - Building the sprot database for the blastx and blastp searches. This seems to do the heavy-lifting with annotations:  
      - `makeblastdb -in uniprot_sprot.pep -dbtype prot`  
      - This makes the file `uniprot_sprot.pep`  
     - Building the pfam, or protein family database, for other annotations:  
        - `gunzip Pfam-A.hmm.gz`  
        - `hmmpress Pfam-A.hmm`  
        - This makes the file `Pfam-A.hmm.gz`
        
  - With the three databases created, and within the trinotate binaries directory, run blastx first on the Trinity fasta  
    - `blastx -query /home/fish_people/redside_dace/dace_transcriptome/dace_transcriptome_Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6`  
    - This step was the longest in the process, taking ~3 days of computer time.   
   - Repeat the process of searching the sprot database with blastp, using the transdecoder .pep output file  
    - `blastp -query /home/fish_people/redside_dace/transdecoder_out/longest_orfs.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6`  
   - Run hmmscan to find known protein residues  
    - `hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm /home/fish_people/redside_dace/transdecoder_out/longest_orfs.pep > pfam.log`
   - [`tmhmm.sh`](https://github.com/BioMatt/redside_dace_RNA/blob/master/transcriptome_annotation/tmhmm.sh) and [`signalp.sh`](https://github.com/BioMatt/redside_dace_RNA/blob/master/transcriptome_annotation/signalp.sh) were both run on the [Cedar cluster](https://docs.computecanada.ca/wiki/Cedar) on [Westgrid](https://www.westgrid.ca/) and [Compute Canada](https://www.computecanada.ca/).  
      - Admin help was needed to increase the maximum number of entries allowed for signalP from 10000 to 2000000 following the [guide](https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required#signalp-v4-free-academic-download)
      
   - With all outputs collected, collect them into an sqlite database for Trinotate   
    - First, the transcripts, gene and transcript map, and protein sequences  
      - `Trinotate Trinotate.sqlite init --gene_trans_map /home/fish_people/redside_dace/dace_transcriptome/dace_Trinity.fasta.gene_trans_map --transcript_fasta /home/fish_people/redside_dace/dace_transcriptome/dace_transcriptome_Trinity.fasta --transdecoder_pep /home/fish_people/redside_dace/transdecoder_out/longest_orfs.pep`  
    - Then load the BLAST protein hits  
         - `Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6`  
    - Load the BLAST transcript hits  
          - `Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6`  
    - Load PFAM domain entries  
          - `Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out`  
    - Load transmembrane domains from tmhmm  
         - `Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out`  
    - Load signal peptide predictions from signalP  
        - `Trinotate Trinotate.sqlite LOAD_signalp signalp.out`  
      
   - Output the annotation report! Unfiltered first
     - `Trinotate Trinotate.sqlite report > dace_annotation_report.xls`
   - Then filtering by a liberal threshold of E values below 0.001 
     -  `Trinotate Trinotate.sqlite report -E 1e-3 > dace_annotation_report_filtered.xls`
