library(Biostrings)
library(msa)
library(stringr)
library(viridis)

sars2genome <- readDNAStringSet("../data/sars-cov-2 genome.fasta") #Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/1798174254 on 1/8/25
reverseComplement(sars2genome) -> rc_sars2genome

sarstwelvemers <- sapply(1:29890, function (x) substr(as.character(sars2genome[[1]]), x, x + 11))
rcsarstwelvemers <- sapply(1:29890, function (x) substr(as.character(rc_sars2genome[[1]]), x, x + 11))

insert1 <- "CTCCTCGGCGGG"
insert2 <- "TCCTCGGCGGGC"

sort(sapply(sarstwelvemers, function (x) sum(strsplit(x, "")[[1]] == strsplit(insert1, "")[[1]])), decreasing = T)[1:10]
sort(sapply(sarstwelvemers, function (x) sum(strsplit(x, "")[[1]] == strsplit(insert2, "")[[1]])), decreasing = T)[1:10]
sort(sapply(rcsarstwelvemers, function (x) sum(strsplit(x, "")[[1]] == strsplit(insert1, "")[[1]])), decreasing = T)[1:10]
sort(sapply(rcsarstwelvemers, function (x) sum(strsplit(x, "")[[1]] == strsplit(insert2, "")[[1]])), decreasing = T)[1:10]
#There are 4 sars2 twelvemers with 9 matches to the insert sequence.
#The presence of 5 9-matching twelvemers suggests that
#matching 9 times can happen by chance (and it doesn't look like these are all repeats of each other; see below)

# CTCCTCGGCGGG  #insert 1
#  TCCTCGGCGGGC #insert 2
# CTACTTGGCGTG  #match 1 (SARS2 to insert 1) #At position 23454 - only ~150 nt away, which is curious. Two of the mutations are nonsynonymous in the context of the insertion site (ACT<->CCT and TGG<->CGG), but the third is synonymous and unlikely (CGT<->CGG, where T<->G is not seen from the MRCA reconstruction to SARS2 (using two sided arrows because it is possible for either the insertion site or the donor site to have mutated after the duplication; but I haven't checked whether they are synonymous/nonsynonymous in the context of the donor. Synonymous changes in the donor might happen, while nonsynonymouus would be assumed to be selected against (while at the insertion site, nonsynonymous might be assumed to be selected FOR)).)
# CTTCTTCGCGGG  #match 2 (SARS2 reverse complement to insert 1) #At position 11630. All three mutations are nonsynonymous (TCT<->CCT and TCG<->CGG)
# CTCCTCCACGGA  #match 3 (SARS2 reverse complement to insert 1) #At position 29528. A P<->R mutation involves both a synonymous and a nonsynonymous change (CCA<->CGG) and the last change is nonsynonymous (A[CA]<->G[CA]). The synonymous change is not super likely (A<->G; G->A is in 2 of 9 4-fold degenerate Gs: low N.)
#  TCTTCTGCAGGC #match 5 (SARS2 to insert 2) #At position 179. All three mutations are nonsynonymous (CTT<->CCT, CTG<->CGG, CAG<->CGG)


#Considering both insert1 and insert2 in a way that accounts for subsequent 12-mers being correlated in their amount of matching to insert1 and insert2
cumsum(round(digits = 6, table(sapply(1:(length(sarstwelvemers) - 1), function (x) max(sum(strsplit(sarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(sarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/29889)[11:1])

#Just within spike
cumsum(round(digits = 6, table(sapply(21563:25372, function (x) max(sum(strsplit(sarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(sarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/length(21563:25372))[10:1])

#Just 12-mers that align the 3rd bases of donor codons with the 3rd bases of insertion codons
cumsum(round(digits = 6, table(sapply(seq(21564, 25374, by = 3), function (x) max(sum(strsplit(sarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(sarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/length(seq(21564, 25374, by = 3)))[10:1])

#Plot of fractions of 12-mers matching
genome.12mer.match <- cumsum(round(digits = 6, table(sapply(1:(length(sarstwelvemers) - 1), function (x) max(sum(strsplit(sarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(sarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/29889)[11:1])
spike.12mer.match <- cumsum(round(digits = 6, table(sapply(21563:25372, function (x) max(sum(strsplit(sarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(sarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/length(21563:25372))[10:1])
spike.sameframe.12mer.match <- cumsum(round(digits = 6, table(sapply(seq(21564, 25374, by = 3), function (x) max(sum(strsplit(sarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(sarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/length(seq(21564, 25374, by = 3)))[10:1])

plot(c(0:12), genome.12mer.match[c(1, 1, 1, 2:11)], type = "l")
lines(c(0:12), spike.12mer.match[c(1, 1, 1, 2, 2, 3:10)], col = "orange")
lines(c(0:12), spike.sameframe.12mer.match[c(1, 1, 1, 2, 2, 3:10)], col = "blue")

#Zoom in to the high-matching low-frequency 12-mers
plot(ylim = c(0, 0.15), xlim = c(0, 8), c(0:12), genome.12mer.match[c(1, 1, 1, 2:11)], type = "l")
lines(c(0:12), spike.12mer.match[c(1, 1, 1, 2, 2, 3:10)], col = "orange")
lines(c(0:12), spike.sameframe.12mer.match[c(1, 1, 1, 2, 2, 3:10)], col = "blue")

#Looking at 12-mers on the negative strand

#Considering both insert1 and insert2 in a way that accounts for subsequent 12-mers being correlated in their amount of matching to insert1 and insert2
cumsum(round(digits = 6, table(sapply(1:(length(rcsarstwelvemers) - 1), function (x) max(sum(strsplit(rcsarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(rcsarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/29889)[10:1])
cumsum(round(digits = 6, table(sapply(4520:8330, function (x) max(sum(strsplit(rcsarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(rcsarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/length(4520:8330))[9:1])
cumsum(round(digits = 6, table(c(sapply(rcsarstwelvemers[seq(4519, 8329, by = 3)], function (x) sum(strsplit(x, "")[[1]] == strsplit(insert1, "")[[1]])),
                                 sapply(rcsarstwelvemers[seq(4520, 8330, by = 3)], function (x) sum(strsplit(x, "")[[1]] == strsplit(insert2, "")[[1]]))))/2542)[8:1])

#plots
plot(c(0:12), c(rep(0, 3), cumsum(round(digits = 6, table(sapply(1:(length(rcsarstwelvemers) - 1), function (x) max(sum(strsplit(rcsarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(rcsarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/29889)[10:1])), type = "l")
lines(c(0:12), c(rep(0, 4), cumsum(round(digits = 6, table(sapply(4520:8330, function (x) max(sum(strsplit(rcsarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(rcsarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/length(4520:8330))[9:1])), col = "orange")
lines(c(0:12), c(rep(0, 5), cumsum(round(digits = 6, table(c(sapply(rcsarstwelvemers[seq(4519, 8329, by = 3)], function (x) sum(strsplit(x, "")[[1]] == strsplit(insert1, "")[[1]])),
                                                             sapply(rcsarstwelvemers[seq(4520, 8330, by = 3)], function (x) sum(strsplit(x, "")[[1]] == strsplit(insert2, "")[[1]]))))/2542)[8:1])), col = "blue")

#zoom
plot(ylim = c(0, 0.01), xlim = c(0, 6), c(0:12), c(rep(0, 3), cumsum(round(digits = 6, table(sapply(1:(length(rcsarstwelvemers) - 1), function (x) max(sum(strsplit(rcsarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(rcsarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/29889)[10:1])), type = "l")
lines(c(0:12), c(rep(0, 4), cumsum(round(digits = 6, table(sapply(4520:8330, function (x) max(sum(strsplit(rcsarstwelvemers[x], "")[[1]] == strsplit(insert1, "")[[1]]), sum(strsplit(rcsarstwelvemers[x+1], "")[[1]] == strsplit(insert2, "")[[1]]))))/length(4520:8330))[9:1])), col = "orange")
lines(c(0:12), c(rep(0, 5), cumsum(round(digits = 6, table(c(sapply(rcsarstwelvemers[seq(4519, 8329, by = 3)], function (x) sum(strsplit(x, "")[[1]] == strsplit(insert1, "")[[1]])),
                                                             sapply(rcsarstwelvemers[seq(4520, 8330, by = 3)], function (x) sum(strsplit(x, "")[[1]] == strsplit(insert2, "")[[1]]))))/2542)[8:1])), col = "blue")


library(DescTools) #Confidence intervals for binomial proportions. https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Jeffreys_interval
BinomCI(6, 85, method = "jeffreys") # 6 A -> T mutations out of 85 4-fold degenerate As in the locally non-recombining region of spike
BinomCI(0, 150, method = "jeffreys") # 0 T -> G mutations out of 150 4-fold degenerate Ts in the locally non-recombining region of spike
BinomCI(28, 278, method = "jeffreys")

readDNAMultipleAlignment("../data/pekar-segment11-fasta-alignment-masked.fasta") -> pekaraln #aligned segment 11 from Pekar's 2022 masked recCA (segment 11 according to Temmam 2022) to the equivalent section of Wuhan-1 using Clustal omega.
table(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "A")])
table(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "T")])
table(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "C")])
table(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "G")])
BinomCI(4, 538, method = "jeffreys")
BinomCI(0, 608, method = "jeffreys")

barplot(ylim = c(0, 100),
        c(sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "A")] == x)))/length(which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "A"))* 100,
          sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "C")] == x)))/length(which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "C"))* 100,
          sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "G")] == x)))/length(which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "G"))* 100,
          sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "T")] == x)))/length(which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "T"))* 100))

basecounts <- c(sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "A")] == x))),
                sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "C")] == x))),
                sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "G")] == x))),
                sapply(c("A", "C", "G", "T"), function (x) length(which(strsplit(as.character(pekaraln@unmasked$`EPI_ISL_406798_Human_Wuhan/WH01/2019_China_2019`), "")[[1]][which(strsplit(as.character(pekaraln@unmasked$pekar), "")[[1]] == "T")] == x))))
sapply(1:16, function (x) text(x*1.2, 20, basecounts[x]))

#fourfold degenerate sites
pekarcodingseq <- "AGATTTCCTAATATTACAAACTTATGCCCTTTTGGTGAAGTTTTTAACGCCACCACATTCGCATCAGTTTATGCTTGGAACAGAAAGAGAATTAGCAACTGTGTTGCNGATTATTCTNTCCTATATAATTCCACTTCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTTTATGCAGATTCATTTGTAGTTAGAGGTGATGAAGTCAGACAAATTGCTCCAGGANAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGTGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCTGGCAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCTTATGGTTTCCACCCTACTAATGGTGTTGGTTACCAACCATATAGGGTAGTAGTACTATCTTTTGAGCTTCTAAATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAACTTGGTTAAAAACAAATGTGTCAATTTCAACTTTAATGGTTTAACTGGCACAGGTGTTCTTACAGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGTAGAGACATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACCTCTAACCAGGTTGCTGTTCTTTATCAGGATGTTAACTGCACAGAAGTCCCTGTTGCTATCCATGCAGACCAACTTACTCCTACTTGGCGTGTTTACTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAACATGTCAACAACTCGTATGAGTGTGACATACCTATTGGTGCAGGAATATGCGCCAGTTATCAGACTCAAACTAATTC------------ACGTAGTGTAGCCAGTCAATCCATTATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCTACAAATTTTACTATTAGTGTTACCACAGAAATTCTACCTGTGTCTATGACTAAGACATCGGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGAGTGCAGCAACCTTTTGTTGCAATATGGCAGTTTTTGCACACAACTAAATCGTGCTTTAACTGGAATAGCTGTTGAACAAGACAAAAACACACAAGAAGTTTTTGCTCAAGTCAAACAAATTTACAAGACACCACCAATTAAAGATTTTGGTGGTTTCAATTTTTCACAAATATTACCAGATCCATCAAAACCAAGCAAGAGGTCATTTATTGAGGATTTACTTTTCAACAAAGTGACACTTGCTGATGCTGGCTTCATCAAACAATATGGTGATTGCCTTGGTGATATTGCTGCTAGAGATCTTATTTGTGCACAAAAGTTTAATGGCCTTACTGTTCTGCCACCTTTGCTCACAGATGAAATGATTGCTCAATACACTTCTGCACTATTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAGGTGCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAGTTACACAGAATGTTCTCTATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTATTGGCAAAATTCAAGACTCACTT"
sars2codingseq <- "AGATTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGACATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGATGTTAACTGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCCACAAATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGAATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCTTTAACTGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAATTTACAAAACACCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGATCCATCAAAACCAAGCAAGAGGTCATTTATTGAAGATCTACTTTTCAACAAAGTGACACTTGCAGATGCTGGCTTCATCAAACAATATGGTGATTGCCTTGGTGATATTGCTGCTAGAGACCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGAAATGATTGCTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAGGTGCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAGTTACACAGAATGTTCTCTATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTATTGGCAAAATTCAAGACTCACTT"
table(sapply(which(sapply(1:611, function (x) substr(pekarcodingseq, x*3 - 2, x*3) %in% c("TCA", "CTA", "CCA", "CGA", "ACA", "GTA", "GCA", "GGA"))), function (y) if (substr(pekarcodingseq, y*3 - 2, y*3 - 1) == substr(sars2codingseq, y*3 - 2, y*3 - 1)) substr(sars2codingseq, y*3, y*3) else F))
BinomCI(4, 79, method = "jeffreys") #Codon 19 mutated from ACN to AGN, and codon 436 mutated from CTN to TTN: both are no longer 4-fold degenerate.
table(sapply(which(sapply(1:611, function (x) substr(pekarcodingseq, x*3 - 2, x*3) %in% c("TCC", "CTC", "CCC", "CGC", "ACC", "GTC", "GCC", "GGC"))), function (y) if (substr(pekarcodingseq, y*3 - 2, y*3 - 1) == substr(sars2codingseq, y*3 - 2, y*3 - 1)) substr(sars2codingseq, y*3, y*3) else F))
table(sapply(which(sapply(1:611, function (x) substr(pekarcodingseq, x*3 - 2, x*3) %in% c("TCG", "CTG", "CCG", "CGG", "ACG", "GTG", "GCG", "GGG"))), function (y) if (substr(pekarcodingseq, y*3 - 2, y*3 - 1) == substr(sars2codingseq, y*3 - 2, y*3 - 1)) substr(sars2codingseq, y*3, y*3) else F))
table(sapply(which(sapply(1:611, function (x) substr(pekarcodingseq, x*3 - 2, x*3) %in% c("TCT", "CTT", "CCT", "CGT", "ACT", "GTT", "GCT", "GGT"))), function (y) if (substr(pekarcodingseq, y*3 - 2, y*3 - 1) == substr(sars2codingseq, y*3 - 2, y*3 - 1)) substr(sars2codingseq, y*3, y*3) else F))
BinomCI(0, 152, method = "jeffreys") #Codon 45 mutated from ACT to GCA: GCN is still 4-fold degenerate, so it counts towards T->A. Codon 75 mutated from GTT to ATT: ATT is no longer 4-fold degenerate, so it doesn't count.

basecounts.4fg <- c(72, 1, 2, 4, 0, 28, 0, 4, 2, 0, 7, 0, 6, 8, 0, 138) #In addition to the resolutions above for codons 19, 436, 45, and 75, we resolve codon 534 as no longer 4-fold degenerate (CTG -> TTG)
barplot(ylim = c(0, 100), basecounts.4fg*100/rep(c(79, 32, 9, 152), each = 4))
sapply(1:16, function (x) text(x*1.2, 20, basecounts.4fg[x]))

#Possible insert1 synonymous codings
#SPRRA
#insert1: CN CCN CGN/AGR CGN/AGR G = 4*4*6*6 = 576 possible synonymous codings
#insert2: N CCN CGN/AGR CGN/AGR GC = 576 possible synonymous codings; no different from insert1, so I didn't analyze these.
sprra.syn.codings <- apply(expand.grid(c("CA", "CC", "CG", "CT"),
                                       c("CCA", "CCC", "CCG", "CCT"),
                                       c("CGA", "CGC", "CGG", "CGT", "AGA", "AGG"),
                                       c("CGA", "CGC", "CGG", "CGT", "AGA", "AGG"),
                                       "G"),
                           1, function (x) paste0(x, collapse = ""))
cumsum(table(sapply(sprra.syn.codings, function (x) sum(strsplit(x, "")[[1]] == strsplit("CACGTAGTGTAG", "")[[1]])))/576) #SARS2 adjoining sequence
cumsum(table(sapply(sprra.syn.codings, function (x) sum(strsplit(x, "")[[1]] == strsplit("CACGTAGTGTGG", "")[[1]])))/576) #ancestral adjoining sequence
#The median synonymous coding has five identical residues to the adjoining SARS sequence - the same as the actual SPRRA sequence has. It also have five identical residues to the
#adjoining ancestral sequence, while the actual SPRRA has 6. This happens 89.6% of the time with synonymous SPRRA codings - somewhat unusual, but not crazily so.

#To run iqtree, I first reformatted Temmam's segment11 alignment (GARD définitif_11.fa) as phylip format using the website http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi
#I then ran the following command on biowulf (in the /data/SBGE/meru/furin/ directory) after saving the phylip formatted alignment as segment11_v3: iqtree2 -s segment11_v3 -o L.R.cornut -m GTR+F+G+I -asr -redo
#The substitution model (GTR+F+G+I) was chosen to match what Pekar used to create RecCA, according to the Methods section of his 2022 RecCA paper.
#The outputted tree (segment11_v3.treefile) matched segment 11 from Temmam's supplementary figure 2. Node7 on the tree corresponds to the MRCA of the human and BANAL coronaviruses (verified by visualization on icytree.org).

read.table("../data/segment11_v3.state", stringsAsFactors = F, header = T) -> ancrec
s2anc <- split(ancrec, ancrec$Node)$Node7
t(apply(s2anc[c(1049:1061, 1074:1084, 1088:1090),4:7], 1, cumsum)) -> s1s2ancrec
plot(ylim=c(0,1), xlim=c(1,30), 1, type = "n")
for (x in 1:27) polygon(c(0 + x, 0 + x, 1 + x, 1 + x), c(1, 0, 0, 1), col = "red", border = NA) #These are Sanger sequencing .ab1 file colors (A = green, C = blue, G = black, T = red), but they're obviously not good for colorblind people. Find better colors.
for (x in 1:27) polygon(c(0 + x, 0 + x, 1 + x, 1 + x), c(s1s2ancrec[x,3], 0, 0, s1s2ancrec[x,3]), col = "black", border = NA)
for (x in 1:27) polygon(c(0 + x, 0 + x, 1 + x, 1 + x), c(s1s2ancrec[x,2], 0, 0, s1s2ancrec[x,2]), col = "blue", border = NA)
for (x in 1:27) polygon(c(0 + x, 0 + x, 1 + x, 1 + x), c(s1s2ancrec[x,1], 0, 0, s1s2ancrec[x,1]), col = "green", border = NA)

#Simulate 100 genomes with the nucleotide frequency seen in the SARS-CoV-2 genome; see how often there are 12-nt stretches with N matches to the insert.
nucdist <- table(strsplit(as.character(sars2genome), "")[[1]])
table(sapply(1:2989000, function (x) sum(sample(c("A", "C", "G", "T"), prob = nucdist, replace = T, size = 12) == strsplit(insert1, "")[[1]])))

#The above only looks for insert1. I wanted to look for both insert1 and insert2 while not double-counting sites in the genome that match both. This is what I came up with:

strsplit(insert1, "")[[1]] -> insert1.str
strsplit(insert2, "")[[1]] -> insert2.str

count.matches.to.insert.in.simulated.genomes <- function (replicates, acgt.dist = c(8954, 5492, 5863, 9594), match.minimum = 9) {
  num.matches.vector <- numeric(replicates)
  for (x in 1:replicates) {
    simgenome <- sample(c("A", "C", "G", "T"), prob = acgt.dist, replace = T, size = 29903)
    simgenome.rc <- sapply(simgenome, function (x) c("A", "C", "G", "T")[x == c("T", "G", "C", "A")])[29903:1]
    match.insert1 <- which(sapply(1:29890, function (x) length(which(simgenome[x:(x+11)] == insert1.str))) > (match.minimum - 1))
    match.insert2 <- which(sapply(1:29890, function (x) length(which(simgenome[x:(x+11)] == insert2.str))) > (match.minimum - 1))
    match.insert1.rc <- which(sapply(1:29890, function (x) length(which(simgenome.rc[x:(x+11)] == insert1.str))) > (match.minimum - 1))
    match.insert2.rc <- which(sapply(1:29890, function (x) length(which(simgenome.rc[x:(x+11)] == insert2.str))) > (match.minimum - 1))
    num.matches.vector[x] <- length(unique(c(match.insert1 + 1, match.insert2))) + length(unique(c(match.insert1.rc + 1, match.insert2.rc)))
  }
  num.matches.vector
}

matches.to.inserts <- count.matches.to.insert.in.simulated.genomes(1000)
mean(matches.to.inserts)
save(matches.to.inserts, file = "../data/output_matches-to-inserts-sim.RData")
