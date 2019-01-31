# Pallavi Gupta
# Biol 432
# January 30 2019

# Extract sequences from .ab1 files and create a FASTA file

library(sangerseqR) # load library
library(dplyr)

sequences <- list.files("data/ab1_sequences", full.names = T) # file pathways into a list
  seqname <- gsub("data/ab1_sequences/", "", sequences) # create a list of file names
quality <- read.csv("data/BarcodePlateStats.csv") # load quality score dataset
BADquality <- subset(quality, Ok == "FALSE") # subset poor quality sequences

FASTAseqs <- character(length(1:33)) # create an empty character vector for sequences
position <- 1

for (i in sequences) { # extract sequences from .ab1 files and put into FASTA format
  header <- gsub(".*data/ab1_sequences/", "", i)
  if (header %in% BADquality) {
    next
  }
  ITS <- read.abif(i)
  ITSseq <- sangerseq(ITS)
  SeqX <- makeBaseCalls(ITSseq)
  PSeqX <- SeqX@primarySeq
  FASTAformat <- gsub("^", "\n", PSeqX)
  FASTAformat1 <- gsub("^", header, FASTAformat)
  FASTAformat2 <- gsub("^", "> ", FASTAformat1)
  FASTAformat3 <- gsub("$", "\n", FASTAformat2)
  FASTAseqs [position] <- FASTAformat3
  position = position + 1
}

FinalFormat <- cat(FASTAseqs) # final formating
FinalFormat # print sequences
