# blobtools-light
Light version of the upcoming blobtools package

blobtools-light allows the visualisation of (draft) genome assemblies using TAGC (Taxon-annotated Gc-Coverage) plots (Kumar et al. 2012).

Software requirements:
- Python 2.7+
- Numpy 1.9.1 ('pip install numpy')
- matplotlib 1.4 ('pip install matplotlib')
- NCBI Taxdump ('wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz')

Data requirements:
- Assembly* (FASTA)
- Coverage information (one or more):
  - COV (Tab-separated: "header\tcoverage")
  - BAM/SAM
  - CAS format (CLC mapper)
  * if assembly was generated through Spades, Abyss or Velvet, coverage can be parsed from contig headers (use e.g. -spades SPADESASSEMBLY)
- BLAST result (using BLAST 2.2.29+):
  - suggested command : "blastn -task megablast -query ASSEMBLY.fa -db /blast_db/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -culling_limit 1 -num_threads 16 -evalue 1e-25 -out ASSEMBLY.vs.nt.25cul1.1e25.megablast.out"
  - The important part is "-outfmt '6 qseqid staxids bitscore'"

Making blobplot file:
~/blobtools-light/makeblobs.py -a assembly/assembly.fa -bam mapping/assembly.mapping.bam -taxdb /exports/blast_db/ -blast assembly.vs.nt.25cul1.1e25.megablast.out -o test

Plotting TAGC-file:
~/blobtools-light/plotblobs.py test.blobplot.txt

