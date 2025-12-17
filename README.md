# SARS-CoV-2_spike_insertion_sequence_analysis

Supplementary data and code for the analysis of the SARS-CoV-2 S1/S2 insertion sequence, with a summary of the contents below

**SARS-CoV-2_spike_insertion_sequence_analysis**
- README.md
- data
    - GARD définitif_11.fa - segment11 alignment from Temmam 2022
    - output_matches-to-inserts-sim.RData - Number of matches of 9 nucleotides or more to the insertion sequence in 1000 simulated genomes
    - pekar_MRCA_masked.fasta - Most recent common ancestor for SARS-CoV-2 from Pekar 2025
    - pekar-segment11-fasta-alignment-masked.fasta - Aligned segment11 from Pekar 2022 masked recCA to the equivalent section of Wuhan-1
    - sars-cov-2 genome.fasta - Wuhan-1 SARS-CoV-2 genome downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) on 1/8/25
    - segment11_v3 - Phylip-formatted version of the GARD definitif_11.fa file
    - segment11_v3.state - Ancestral reconstructions outputted by IQ-TREE 2
    - segment11_v3.treefile - Phylogenetic tree of segment11 sequences from Temmam 2022
- figures
    - Figure 1.ai
    - Figure 2.ai
    - Supplementary Figure 1.ai
    - Supplementary Figure 2.ai
    - Supplementary Figure 3.ai
- scripts
    - sars-cov-2_insertion-sequence-analysis.R - Analysis of the insertion sequence and mutation frequency in the SARS-CoV-2 MRCA
