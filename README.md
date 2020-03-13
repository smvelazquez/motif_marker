# motif_marker
-------------------------------------------------------------
By: Samantha Velazquez

The *motif_identification.py* script is capable of taking in a single-line fasta file alongisde a .txt file that contains a list of motifs one would like to investigate. 
Motifs are defined, translated into regular expressions, and used to search the fasta file. All found motifs are then placed into a list where they are marked onto an image later. 
The image is described as follows:

  1. Exons are denoted by large rectangular boxes, to scale. The rest of DNA is visualized as a single thin black line, also to scale. 
  2. Motifs are colored by their respective motif and are marked to scale. 
  3. A legend of motifs is underneath the image indicating which gene correlates to which color.
  4. Gene names are located next to their respective strand for clarity. 
  
This specific motif marker recognizes exon sequences as using capitalized characters (ATUCG) and intron sequences using lowercase (atucg), as per UCSC Genome Browser sequence download format. This script does not support any other characters aside from ATCG and U as the ambiguous characters. 

# Downloads and Requirements
------------------------------------------------------------

This script requires Python 3.7.1. It also requires a single-line fasta file. Please use the following command to convert multiline to single-line. 

```
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' interleaved.fasta > singleline.fasta
```
