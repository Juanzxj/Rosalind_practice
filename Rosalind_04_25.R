#Rosalind 2024-04-24
# Counting DNA Nucleotides
# Problem
# A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.
# An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."
# Given: A DNA string s
# of length at most 1000 nt.
# Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s

countNucleotides <- function(DNAstring){
  nucleotides <- c('A', 'C','G','T')
  counts <- sapply(nucleotides, function(nuc){
    sum(DNAstring==nuc)
  })
  names(counts) <- nucleotides
  return(counts)
}
DNAexample <- readLines("rosalind_dna.txt")
dnastring <- unlist(strsplit(DNAexample, ""))
countNucleotides <- countNucleotides(dnastring)
countNucleotides
################################################
# Transcribing DNA into RNA 
# Problem
# An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
# Given a DNA string t
# corresponding to a coding strand, its transcribed RNA string u
# is formed by replacing all occurrences of 'T' in t
# with 'U' in u
# Given: A DNA string t
# having length at most 1000 nt.
# Return: The transcribed RNA string of t

transcribeDNAtoRNA <- function(dnaString){
  RNAsring <- chartr("T", "U", dnaString)
  return(RNAsring)
}
dnaExample <- readLines("rosalind_rna.txt")
RNAtranscription <- transcribeDNAtoRNA(dnaExample)
RNAtranscription

################################################
# Complementing a Strand of DNA-
# revc
# Problem
# In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.
# The reverse complement of a DNA string s
# is the string sc
# formed by reversing the symbols of s
# , then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").
# Given: A DNA string s
# of length at most 1000 bp.
# Return: The reverse complement sc of s

reverseComplementDNA <- function(dnaString) {
  # Reverse the DNA string
  reversedString <- rev(unlist(strsplit(dnaString, "")))
  
  # Map each nucleotide to its complement
  complements <- setNames(c("T", "G", "A", "C"), c("A", "C", "T", "G"))
  complementedString <- sapply(reversedString, function(nuc) {
    complements[nuc]
  })
  
  # Collapse the list of characters into a single string
  reversedComplement <- paste0(complementedString, collapse = "")
  
  return(reversedComplement)
}

dnaExample <- readLines("rosalind_revc.txt")
dnaExample <-paste0(dnaExample, collapse = "")
reverseComplementDNA <- reverseComplementDNA(dnaExample)
reverseComplementDNA
################################################
# Computing GC Content
# Identifying Unknown DNA Quicklyclick to expand
# Problem
# The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.
# DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.
# In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.
# Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
# Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.

# Function to calculate GC-content
calculateGCContent <- function(dnaString) {
  if (nchar(dnaString) == 0) return(0)  # Return 0 if string is empty to avoid division by zero
  dnaVector <- unlist(strsplit(dnaString, ""))
  gcCount <- sum(dnaVector %in% c('G', 'C'))
  gcContent <- (gcCount / nchar(dnaString)) * 100
  return(gcContent)
}

#  find the sequence with the highest GC content in a datase
findHighestGCContent <- function(processed_data) {
  maxGCContent <- 0  # Initialize with a very low GC content
  maxGCID <- ""       # Initialize the ID of the sequence with the highest GC content
  
  # Loop through each sequence in the processed data
  for (id in names(processed_data)) {
    currentGCContent <- calculateGCContent(processed_data[[id]])
    
    # Check if the current sequence has higher GC content than the max found so far
    if (currentGCContent > maxGCContent) {
      maxGCContent <- currentGCContent
      maxGCID <- id
    }
  }
  
  # Return the ID of the sequence with the highest GC content and the content value
  return(list(ID = maxGCID, GCContent = maxGCContent))
}

# Read the FASTA data from a file
fastaData <- readLines("rosalind_gc.txt")
# Initialize variables to store the current header and sequences
sequences <- list()
currentHeader <- NULL
currentSequence <- ""
# Process the data to format as requested as list
processed_data <- list()
current_header <- ""
for (line in fastaData) {
  if (startsWith(line, ">")) {
    current_header <- substr(line, 2, nchar(line))  #extracts the header text (excluding the > character) using the substr function. It starts from the second character (index 2) to the end of the line
    processed_data[[current_header]] <- ""   #empty string if it is the line of header
  } else {
    processed_data[[current_header]] <- paste0(processed_data[[current_header]], line)
  }
}
# Display the processed data
processed_data
results <- findHighestGCContent(processed_data) 

results

################################################################################################################################################  
# Counting Point Mutations 
# Problem
# hamm
# Figure 2. The Hamming distance between these two strings is 7. Mismatched symbols are colored red.
# Given two strings s and t of equal length, the Hamming distance between s and t denoted dH(s,t) is the number of corresponding symbols that differ in s
# and t. See Figure 2.
# Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
# Return: The Hamming distance dH(s,t)

dnaString <- readLines("rosalind_hamm.txt")

dnaString1 <- unlist(strsplit(dnaString[1],split = ""))
dnaString2 <- unlist(strsplit(dnaString[2],split = ""))
dnaString1
dnaString2
num <-0
for (i in 1:length(dnaString1)){
  if (dnaString1[i] != dnaString2[i]){
    num <-num +1
  }
}

print(num)
################################################
# Translating RNA into Protein
# prot
# Problem
# The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.
# The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.
# Given: An RNA string s
# corresponding to a strand of mRNA (of length at most 10 kbp).
# Return: The protein string encoded by s

library(readxl)
# Read the data from the Excel file
RNAcodomtable <- read_excel("RNAcodomtable.xlsx")

RNAStringlines <- readLines("rosalind_prot.txt")
RNAString <- paste(RNAStringlines, collapse = "")
RNAString

# Define the genetic code as a named vector
genetic_code <- c(
  UUU = "F", CUU = "L", AUU = "I", GUU = "V",
  UUC = "F", CUC = "L", AUC = "I", GUC = "V",
  UUA = "L", CUA = "L", AUA = "I", GUA = "V",
  UUG = "L", CUG = "L", AUG = "M", GUG = "V",
  UCU = "S", CCU = "P", ACU = "T", GCU = "A",
  UCC = "S", CCC = "P", ACC = "T", GCC = "A",
  UCA = "S", CCA = "P", ACA = "T", GCA = "A",
  UCG = "S", CCG = "P", ACG = "T", GCG = "A",
  UAU = "Y", CAU = "H", AAU = "N", GAU = "D",
  UAC = "Y", CAC = "H", AAC = "N", GAC = "D",
  UAA = "*", CAA = "Q", AAA = "K", GAA = "E",
  UAG = "*", CAG = "Q", AAG = "K", GAG = "E",
  UGU = "C", CGU = "R", AGU = "S", GGU = "G",
  UGC = "C", CGC = "R", AGC = "S", GGC = "G",
  UGA = "*", CGA = "R", AGA = "R", GGA = "G",
  UGG = "W", CGG = "R", AGG = "R", GGG = "G"
)

# Function to translate RNA into a protein string
translateRNA <- function(rna) {
  # Split the RNA into codons of three nucleotides
  codons <- sapply(seq(1, nchar(rna), by=3), function(i) substr(rna, i, i+2))
  # extracts a substring from the RNA string that starts at i and ends at i+2. This substring is exactly one codon long.
  # Convert codons to amino acids, stop at the first stop codon
  amino_acids <- sapply(codons, function(codon) genetic_code[codon])
  protein_string <- paste0(amino_acids[amino_acids != "*"], collapse = "")
  
  return(protein_string)
}

# Sample RNA string
RNAStringlines <- readLines("rosalind_prot.txt")
rna_string <- paste(RNAStringlines, collapse = "")

# Translate RNA to protein
protein_string <- translateRNA(rna_string)
print(protein_string)
############################################################################################################# 
# Calculating Protein Mass
# Problem
# In a weighted alphabet, every symbol is assigned a positive real number called a weight. A string formed from a weighted alphabet is called a weighted string, and its weight is equal to the sum of the weights of its symbols.
# The standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic mass of the corresponding amino acid.
# Given: A protein string P
# of length at most 1000 aa.
# Return: The total weight of P
# . Consult the monoisotopic mass table.

# Define the monoisotopic mass table for amino acids
monoisotopic_masses <- list(
  A = 71.03711,   # Alanine
  C = 103.00919,  # Cysteine
  D = 115.02694,  # Aspartic acid
  E = 129.04259,  # Glutamic acid
  F = 147.06841,  # Phenylalanine
  G = 57.02146,   # Glycine
  H = 137.05891,  # Histidine
  I = 113.08406,  # Isoleucine
  K = 128.09496,  # Lysine
  L = 113.08406,  # Leucine
  M = 131.04049,  # Methionine
  N = 114.04293,  # Asparagine
  P = 97.05276,   # Proline
  Q = 128.05858,  # Glutamine
  R = 156.10111,  # Arginine
  S = 87.03203,   # Serine
  T = 101.04768,  # Threonine
  V = 99.06841,   # Valine
  W = 186.07931,  # Tryptophan
  Y = 163.06333   # Tyrosine
)
# monoisotopic_masses[["V"]]
# Function to calculate the total weight of a protein string
calculate_protein_weight <- function(protein_string) {
  total_weight <- 0
  for (i in 1:nchar(protein_string)) {
    amino_acid <- substr(protein_string, i, i)
    total_weight <- total_weight + monoisotopic_masses[[amino_acid]]
  }
  return(total_weight)
}
# monoisotopic_masses[[amino_acid]] use key (amino_acid) to get the value in the list
# substr() is a function in R that extracts part of a string. It takes three arguments:
#   string: The string from which to extract a substring.
# start: The starting position of the substring.
# stop: The ending position of the substring.
# Example usage:
# Example protein string
protein_string <- readLines("rosalind_prtm.txt")
protein_string <- paste0(protein_string, collapse = "")
total_weight <- calculate_protein_weight(protein_string)

formatted_weight <- sprintf("%.3f", total_weight)
print(formatted_weight)

###########################################################################################################
# Introduction to the Bioinformatics Armory 
#   Problem
#   This initial problem is aimed at familiarizing you with Rosalind's task-solving pipeline. To solve it, you merely have to take a given DNA sequence and find its nucleotide counts; this problem is equivalent to “Counting DNA Nucleotides” in the Stronghold.
# Of the many tools for DNA sequence analysis, one of the most popular is the Sequence Manipulation Suite. Commonly known as SMS 2, it comprises a collection of programs for generating, formatting, and analyzing short strands of DNA and polypeptides.
# One of the simplest SMS 2 programs, called DNA stats, counts the number of occurrences of each nucleotide in a given strand of DNA. An online interface for DNA stats can be found here.
# Given: A DNA string s
#  of length at most 1000 bp.
# Return: Four integers (separated by spaces) representing the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s
# . Note: You must provide your answer in the format shown in the sample output below.

library(stringr)
dna_sequence <- readLines("rosalind_ini.txt")
dna_sequence <- paste0(dna_sequence, collapse = "")
# Count occurrences of each nucleotide
count_A <- str_count(dna_sequence, "A")
count_C <- str_count(dna_sequence, "C")
count_G <- str_count(dna_sequence, "G")
count_T <- str_count(dna_sequence, "T")

# Print the counts
print(paste(count_A, count_C, count_G, count_T))


########################################################################################################### 
# perm
# Load necessary library
library(combinat)
# Generates all permutations of numbers from 1 to 𝑛  using the permn function.
# Function to generate and print all permutations of numbers 1 to n
generate_permutations <- function(n) {
  # Generate all permutations
  permutations <- permn(n)
  
  # Print the number of permutations
  cat("Number of permutations:", length(permutations), "\n")
  
  # Print each permutation on a new line
  for (perm in permutations) {
    cat(perm, "\n")
  }
}
# The cat() function in R is used for concatenating and printing objects to the console or to a file.
# Example usage with n = 3
generate_permutations(5)

#######################################################################################################
# Finding a Motif in DNA
# subs

find_substring_positions <- function(s, t) {
  # Initialize a vector to store the starting positions of each match
  positions <- numeric()
  
  # Compute the length of string t
  length_t <- nchar(t)
  
  # Loop through string s up to the point where t can still fit
  for (i in 1:(nchar(s) - length_t + 1)) {
    # Extract the substring of length t starting from position i
    substring_s <- substr(s, i, i + length_t - 1)
    
    # If the substring matches t, store the position
    if (substring_s == t) {
      positions <- c(positions, i)
    }
  }
  
  # Return the positions
  return(positions)
}

# Sample dataset
s <- "AAAGAGATTGAAACAGCGCTTGCACAGCGCAACAGCGCCACAGCGCACAGCGCACAGCGCAGACAGCGCTACAGCGCACAGCGCACAGCGCCTGGACAGCGCACACAGCGCTTCCTACAGACGGCTGGGACAGCGCACAGCGCCGCTAGACAGCGCAAACAGCGCCACAGCGCACACAAACTACACAGCGCACACAGCGCACAGCGCACAGCGCCTACACAGCGCTGAAACAGCGCCAGATGACAGCGCAGACAGCGCAGAGGACAGCGCACAGCGCAACAGCGCGACAGCGCGGACAGCGCACAGCGCACAGCGCATATACAGCGCCACAGCGCTACAGCGCTAATACGTTCACAGCGCGACAGCGCACAGCGCACAGCGCACAGCGCTGCACAGCGCACAGCGCAACAGCGCGGGCACAGCGCCACAGCGCGTCCCAAACAGCGCTGCCGCACAGCGCACAGCGCACAGCGCGACAGCGCGACAGCGCTACGACAGCGCCATGGGACAGCGCAACAGCGCTACAGCGCCAAAGACAGCGCTTCACAGCGCCACAGCGCAGTGTGGTTTCACAGCGCCACAGCGCAAAGGACAGCGCAACACAGCGCGGATACAGCGCCACAGCGCCCACAGCGCAGATTGCGACAGCGCAGTACAGCGCACAGCGCAGAAGGGACAGCGCACAGCGCGGGACACAGCGCTACAGCGCCCTTGGACAGCGCTCACAGCGCTACAGCGCACAGCGCGGAGTACAGCGCAATATATAGAACAGCGCACAGCGCAATCGCACAGCGCAGATCGGCCACAGCGCGGGACAGCGCCGACAGCGCACAGCGCGGTGCAACAGCGCTGACAGCGCAGACAGCGCACAGCGCTGGACCCAACAGCGCGCGTGACAGCGCCGTAACAGCGCAACAGCGCCTACAGCGCACAGCGCAATTAGACAGCGCACAGCGCGATACAGCGCACAGCGC"
t <- "ACAGCGCAC"

# Call the function and print the output
result_positions <- find_substring_positions(s, t)
cat(result_positions)

#######################################################################################################
#######################################################################################################
# grph-
# Overlap Graphs
# Function to parse FASTA data
parse_fasta <- function(fasta) {
  # Split the data into entries
  entries <- unlist(strsplit(fasta, ">"))
  entries <- entries[nchar(entries) > 0]
  
  # Create lists for labels and sequences
  labels <- vector("character", length(entries))
  sequences <- vector("character", length(entries))
  
  # Process each entry
  for (i in seq_along(entries)) {
    lines <- unlist(strsplit(entries[i], "\n")) #
    # splits the string by newline characters to separate the label (first line) from the sequence (subsequent lines).
    labels[i] <- trimws(lines[1])
    #The label is extracted directly (after trimming whitespace).
    sequences[i] <- paste(trimws(lines[-1]), collapse = "")
    #The sequence lines are joined together into a single string, removing any internal whitespace.
  }
  
  return(list(labels = labels, sequences = sequences))
}
#parse_fasta to get the labels and sequences from the input data.
# Function to find overlaps
find_overlaps <- function(fasta_data) {
  parsed_data <- parse_fasta(fasta_data)
  #takes a single string containing FASTA formatted data as input
  labels <- parsed_data$labels
  sequences <- parsed_data$sequences
  n <- length(labels)
  #For each suffix, it checks against the prefixes (first three characters) of all other sequences.
  for (i in 1:n) {
    suffix <- substr(sequences[i], nchar(sequences[i]) - 2, nchar(sequences[i]))
    for (j in 1:n) {
      if (i != j) {
        prefix <- substr(sequences[j], 1, 3)
        if (suffix == prefix) {
          cat(labels[i], labels[j], "\n")
        }
      }
    }
  }
}

# Example FASTA data (input)
fasta_data <- "
>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG

"

# Find and print overlaps
find_overlaps(fasta_data)

####################################################################################################### 
####################################################################################################### 
# Enumerating k-mers Lexicographically
# lexf
#   Problem
#   Assume that an alphabet 𝒜 has a predetermined order; that is, we write the alphabet as a permutation 𝒜=(a1,a2,…,ak)  where a1<a2<⋯<ak
#   For instance, the English alphabet is organized as (A,B,…,Z)
#   Given two strings s and t having the same length n we say that s
#   precedes t in the lexicographic order (and write s<Lext) if the first symbol s[j] that doesn't match t[j] satisfies sj<tj in 𝒜
# Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer n (n≤10 ).
# Return: All strings of length n that can be formed from the alphabet, ordered lexicographically (use the standard order of symbols in the English alphabet).

generate_strings <- function(alphabet, n) {
  # Convert the alphabet string into a vector of characters
  symbols <- strsplit(alphabet, " ")[[1]] # splits the alphabet string wherever there is a space.
  
  # Use expand.grid to generate all possible combinations of length n
  combinations <- expand.grid(rep(list(symbols), n))
  
  # Convert the data frame rows to single strings
  strings <- apply(combinations, 1, paste, collapse = "")
  #1: This argument tells apply to operate across rows. If you were to use 2 instead, the function would operate across columns.
  # Print the results lexicographically
  # Sort the strings lexicographically
  # Sort the strings based on the first symbol, then lexicographically
  sorted_strings <- strings[order(substr(strings, 1, 1), strings)]
  
  # Print the results lexicographically sorted by the first symbol
  cat(sorted_strings, sep = "\n")
}

# Example input
alphabet <- "A B C D E"
n <- 4

# Generate and print all strings of length n
generate_strings(alphabet, n)

####################################################################################################### 
####################################################################################################### 
#SPLC
# RNA Splicing
# Problem
# After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.
# Given: A DNA string s
# (of length at most 1 kbp) and a collection of substrings of s
# acting as introns. All strings are given in FASTA format.
# Return: A protein string resulting from transcribing and translating the exons of s
# . (Note: Only one solution will exist for the dataset provided.)


#parse_fasta Function:
parse_fasta <- function(fasta) {
  entries <- strsplit(fasta, "\n", fixed = TRUE)  # Split the data by newlines
  entries <- entries[[1]][nchar(entries[[1]]) > 0]  # Remove any empty entries
  names <- vector("character")
  sequences <- vector("character")
  current_name <- NULL
  current_sequence <- ""
  
  for (line in entries) {
    if (startsWith(line, ">")) {
      if (!is.null(current_name)) {
        names <- c(names, current_name)
        sequences <- c(sequences, current_sequence)
        current_sequence <- ""
      }
      current_name <- substring(line, 2)  # Get the name without '>'
    } else {
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  if (!is.null(current_name)) {
    names <- c(names, current_name)
    sequences <- c(sequences, current_sequence)
  }
  
  names(sequences) <- names
  return(sequences)
}



#remove_introns 
remove_introns <- function(dna, introns) {
  for (intron in introns) {
    dna <- gsub(intron, "", dna)
  }
  return(dna)
}

transcribe_dna_to_rna <- function(dna) {
  return(gsub("T", "U", dna))
}

# Function to translate RNA to protein, correcting function name and logic
translate_rna_to_protein <- function(rna, genetic_code) {
  protein <- ""
  if (nchar(rna) >= 3) {  # Check if RNA is long enough to have at least one codon
    for (i in seq(1, nchar(rna) - 2, by = 3)) {
      codon <- substr(rna, i, i + 2)
      amino_acid <- genetic_code[codon]
      if (amino_acid == "*") break  # Stops at a stop codon
      protein <- paste0(protein, amino_acid)
    }
  } else {
    return(NULL)  # Return NULL or an appropriate value if RNA is too short
  }
  return(protein)
}
# Sample dataset
genetic_code <- c("AAA"="K", "AAC"="N", "AAG"="K", "AAU"="N",
                  "ACA"="T", "ACC"="T", "ACG"="T", "ACU"="T",
                  "AGA"="R", "AGC"="S", "AGG"="R", "AGU"="S",
                  "AUA"="I", "AUC"="I", "AUG"="M", "AUU"="I",
                  "CAA"="Q", "CAC"="H", "CAG"="Q", "CAU"="H",
                  "CCA"="P", "CCC"="P", "CCG"="P", "CCU"="P",
                  "CGA"="R", "CGC"="R", "CGG"="R", "CGU"="R",
                  "CUA"="L", "CUC"="L", "CUG"="L", "CUU"="L",
                  "GAA"="E", "GAC"="D", "GAG"="E", "GAU"="D",
                  "GCA"="A", "GCC"="A", "GCG"="A", "GCU"="A",
                  "GGA"="G", "GGC"="G", "GGG"="G", "GGU"="G",
                  "GUA"="V", "GUC"="V", "GUG"="V", "GUU"="V",
                  "UAA"="*", "UAC"="Y", "UAG"="*", "UAU"="Y",
                  "UCA"="S", "UCC"="S", "UCG"="S", "UCU"="S",
                  "UGA"="*", "UGC"="C", "UGG"="W", "UGU"="C",
                  "UUA"="L", "UUC"="F", "UUG"="L", "UUU"="F")

# Define functions: parse_fasta, remove_introns, transcribe_dna_to_rna, translate_rna_to_protein
# fasta_input <- readLines("rosalind_splc.txt")
# fasta_input <- paste0(fasta_input, collapse = "")
fasta_input <- ">Rosalind_6825
ATGTCACATTTATGGCTGAGTAGACCTTCTTCTCATCCTAGACCAGGACACGACAGTTAT
ACTCCATTATCACGTCTAGACCGGCTCTTCCTAGGTGGTGGCCGGCAGGTATCGGGCGTG
TTAACCGCTGATGTCAAAGTCACTAGCAGTCCACTCACATGAGACAACTCCACTCAGAAC
CAAGTGTGATATCACGGTGGGCATTAGATCCACGGATGGCGCCTCGCAACGATTTTTCAT
ATATGTATATTTTACAGTTACAAACTATAATACCATGCTAGGTGTGATGTATGTCGCCTC
TAGGGTGTCCTTCTCTGCTGCCATGCTATATTACGCTCCGCGTAGGACTATCAGAGTCAT
CGACGACGCTTCACAGATAGATTCACCACACTGAGCGACGAGAAAGGGTCATCAATGATT
GTGTCTTTGAAATGCTTGGAGTAGAGTGTCCAAGCAACGAACAGCTCTGTTCTGAGCCCC
TCGAACGCGCACGGCTTATCAAGGCGGGTTATTTAAAAATCACAATGTATTTCCAATGGG
CGCACGACCCCTCCGGCCGTGGACGAACGAGACAGGCGAAAGTCATTGAGTTATCCCGCG
CAGTAGATCAATCTTAGGGGCCCCCACCTTTCCCGATTGATCAGTAGCAATTCTGTCTTT
CTACAAACGCCCTTCTGGAGCTCTGTGGCCGTTGCCGCGATCCGGCTCCCTCGATGGAAC
ACCCATAAAAATCAATAGTAATTATGACTCGATAGGAGCTGACCAGTCCTCCACACACTT
AGATTCGCGTGTGGTCCGATGGCCATGGAAAAAACAGCCGAGACAGGCCACGAATAACGA
GGAAGATTTTCACGTAGGTGCATGCCAGGCTCGTTGCCGAAAGGGGCCACCCGTCCGTGT
TAAGTTCACGGCCAATGGTGGGTTTGAAACTCAGGATCAGCAAACAGTATCGAGCAGGTG
ACTTGCTTTAATGCTACATGCGTGACACTGTAG
>Rosalind_9261
CAATAGTAATTATGACTCGATAGGAGCTGACCAG
>Rosalind_7663
TGGCCATGGAAAAAACAGCCGAGACAGGCCACGAATAACGAGGAAGAT
>Rosalind_5858
GCATTAGATCCACGGATGGCGCCTCGCAACGATTTTTCATATATGTA
>Rosalind_4718
AGTAGCAATTCTGTCTTTCTACAAACGCCCTTCTGGAGCTCT
>Rosalind_0544
AACTATAATACCATGCTAGGT
>Rosalind_2605
TATCGGGCGTGTTAACCGCTGATGTCAAAGTCACTAGCAGTCC
>Rosalind_0647
GGCGAAAGTCATTGAGTTATCCCGCGCAGTAG
>Rosalind_7584
ACTCCATTATCACGTCTAGACCGGCTCTTCCTAG
>Rosalind_5536
CCCGTCCGTGTTAAGTTCA
>Rosalind_3549
ATCACAATGTATTTCCAATG
>Rosalind_7760
TGCTTGGAGTAGAGTGTCCAAGCAACGAA
>Rosalind_7599
AGGATCAGCAAACAGT
>Rosalind_8844
GTGTCCTTCTCTGCTGCCATGCTATATTACGCTC
>Rosalind_5139
TCGACGACGCTTCACAGATAGATTCACCACACTG
"


# Process the input
sequences <- parse_fasta(fasta_input)
print("Sequences extracted:")
print(sequences)

main_dna <- sequences[[1]]  # Assuming the first sequence is the main DNA sequence
introns <- sequences[-1]  # The rest are introns

print("Main DNA:")
print(main_dna)
print("Introns:")
print(introns)

# Remove introns and transcribe to RNA
exons <- remove_introns(main_dna, introns)
print("Exons after removing introns:")
print(exons)

rna_sequence <- transcribe_dna_to_rna(exons)  # Transcribe DNA to RNA
print("RNA Sequence:")
print(rna_sequence)

# Translate RNA to protein
protein <- translate_rna_to_protein(rna_sequence, genetic_code)
print("Translated Protein:")
print(protein)


# ## read FASTA file-# splc- read FASTA file from the txt file, and cleaning
# #To read a FASTA formatted text file and parse it into a structured format in R
# 
# fasta_data <- readLines("rosalind_splc.txt")
# parse_fasta <- function(fasta_data) {
#   sequences <- list()  # Initialize an empty list to store sequences
#   current_name <- NULL
#   
#   for (line in fasta_data) {
#     if (startsWith(line, ">")) {  # Check if the line is a label
#       current_name <- substring(line, 2)  # Remove the '>'
#       sequences[[current_name]] <- ""  # Initialize an empty sequence for this label
#     } else if (!is.null(current_name)) {
#       sequences[[current_name]] <- paste0(sequences[[current_name]], line)  # Append the sequence
#     }
#   }
#   
#   return(sequences)
# }
# 
# # Call the function with the data read from the file
# parsed_sequences <- parse_fasta(fasta_data)
# for (name in names(parsed_sequences)) {
#   cat(name, "\n", "\"", parsed_sequences[[name]], "\"", "\n\n")
# }

####################################################################################################### 
####################################################################################################### 
# Finding a Shared Motif
# Problem
# A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".
# Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".
# Given: A collection of k
# (k≤100
# ) DNA strings of length at most 1 kbp each in FASTA format.
# Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)

##read FASTA file
#To read a FASTA formatted text file and parse it into a structured format in R

fasta_data <- readLines("rosalind_lcsm.txt")
parse_fasta <- function(fasta_data) {
  sequences <- list()  # Initialize an empty list to store sequences
  current_name <- NULL
  
  for (line in fasta_data) {
    if (startsWith(line, ">")) {  # Check if the line is a label
      current_name <- substring(line, 2)  # Remove the '>'
      sequences[[current_name]] <- ""  # Initialize an empty sequence for this label
    } else if (!is.null(current_name)) {
      sequences[[current_name]] <- paste0(sequences[[current_name]], line)  # Append the sequence
    }
  }
  
  return(sequences)
}

# Call the function with the data read from the file
parsed_sequences <- parse_fasta(fasta_data)
parsed_sequences
for (name in names(parsed_sequences)) {
  cat(name, "\n", "\"", parsed_sequences[[name]], "\"", "\n\n")
}
#"\n": This is a newline character. It's used here to move to the next line after printing name.
#"\"": This is an escaped double quote character. 
#In R, you need to use the backslash (\) to escape certain characters like quotes within strings. 
#The purpose here is to include literal quotation marks in the output around the sequence data.
##########Finding the longest common substring
longest_common_substring <- function(sequences) {
  if (length(sequences) == 0) return("")
  
  # Start with the first sequence
  first_seq <- sequences[[1]]
  max_length <- nchar(first_seq)
  longest_sub <- ""
  
  # Check all substrings of the first sequence
  #The outer loop (start) begins at the first character of first_seq and moves to the end
  for (start in 1:max_length) {
    #the inner loop (length) iteratively increases the length of the substring starting from start. 
    for (length in 1:(max_length - start + 1)) {
      substring <- substr(first_seq, start, start + length - 1)
      found_in_all <- TRUE
      
      for (seq in sequences[-1]) {
        if (!grepl(substring, seq)) {
          found_in_all <- FALSE
          break
        }
      }
      
      if (found_in_all && nchar(substring) > nchar(longest_sub)) {
        longest_sub <- substring
      }
    }
  }
  
  return(longest_sub)
}


# Output the result
result <- longest_common_substring(parsed_sequences)
cat("The longest common substring is:", result, "\n")

####################################################################################################### 
####################################################################################################### 
#Parse the FASTA Data
fasta_input <- "
>Rosalind_14
ACGTACGTGACG
>Rosalind_18
GTA
"

# Function to parse FASTA data
parse_fasta <- function(fasta_data) {
  sequences <- strsplit(fasta_data, "\n")[[1]]
  sequences <- sequences[nchar(sequences) > 0]  # Remove any empty lines
  list_of_sequences <- list()
  current_label <- NULL
  
  for (line in sequences) {
    if (startsWith(line, ">")) {
      current_label <- substring(line, 2)
      list_of_sequences[[current_label]] <- ""
    } else if (!is.null(current_label)) {
      list_of_sequences[[current_label]] <- paste0(list_of_sequences[[current_label]], line)
    }
  }
  
  return(list_of_sequences)
}

# Calling the parsing function
parsed_sequences <- parse_fasta(fasta_input)
s <- parsed_sequences[[1]]
t <- parsed_sequences[[2]]
s
t
find_subsequence_indices <- function(s, t) {
  positions <- integer(length(t))  # To store positions
  current_index <- 1              # Initialize current_index
  
  for (i in seq_along(t)) {
    char_to_find <- t[i]
    found <- FALSE
    
    # Loop through remaining characters in s starting from current_index
    while (current_index <= nchar(s)) {
      if (substr(s, current_index, current_index) == char_to_find) {
        positions[i] <- current_index
        current_index <- current_index + 1  # Move to next character
        found <- TRUE
        break  # Found the character, break the loop
      }
      
    }
    
    # If no character is found and while loop completes
    if (!found) {
      return(integer(0))  # Return empty if any character is not found
    }
  }
  current_index <- current_index + 1
  return(positions)
}


#  

# Example call to the function
positions <- find_subsequence_indices(s, t)

# Print the positions
if (length(positions) > 0) {
  cat("Positions:", positions, sep=" ")
} else {
  cat("No valid subsequence found.")
}


####################################################################################################### 
####################################################################################################### 
# sign-Enumerating Oriented Gene Orderings 

library(combinat)

# Function to generate all signed permutations of a given number n
generate_signed_permutations <- function(n) {
  perms <- permn(n)  # Generate all permutations of numbers 1 to n
  
  # Function to generate all sign combinations for a single permutation
  generate_signs <- function(perm) {
    signs <- expand.grid(rep(list(c(-1, 1)), length(perm)))
    signed_perms <- apply(signs, 1, function(sign) {
      return(sign * perm)  # Apply sign directly here
    })
    return(matrix(signed_perms, ncol = length(perm), byrow = TRUE))
  }
  
  # Apply the sign generation to all permutations
  signed_perms_list <- do.call(rbind, lapply(perms, generate_signs))  # Use rbind to correctly merge
  return(split(signed_perms_list, seq(nrow(signed_perms_list))))  # Split into a list by rows
}

# Example usage for n = 2
n <- 3
signed_perms <- generate_signed_permutations(n)


# Generate signed permutations for n = 2 (example usage)
n <-3
signed_perms <- generate_signed_permutations(n)
signed_perms
as.matrix(signed_perms, byrow=T)

# Assuming signed_perms is a list of vectors
perms_matrix <- do.call(rbind, signed_perms)
#View(perms_matrix)
# Optionally, convert it to a data frame if needed (for better formatting or additional handling)
perms_df <- as.data.frame(perms_matrix)
# Write the data frame to a text file

write.table(perms_df, "signed_permutations.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#check other methold
####################################################################################################### 
#######################################################################################################   
# lexv-Ordering Strings of Varying Length Lexicographically

generate_strings <- function(alphabet, max_length) {
  # Helper function for recursive generation of strings
  recursive_generate <- function(current_string, current_length) {
    if (current_length == max_length) {
      return(current_string)
    }
    
    # Generate strings by appending each alphabet symbol to the current string
    results <- character(0)
    for (symbol in alphabet) {
      new_string <- paste0(current_string, symbol)
      results <- c(results, recursive_generate(new_string, current_length + 1))
    }
    
    return(c(current_string, results))
  }
  
  # Initialize the process with an empty string and length 0
  results <- character(0)
  for (symbol in alphabet) {
    results <- c(results, recursive_generate(symbol, 1))
  }
  
  return(results)
}

# Sample input
alphabet <- c("D", "Q", "A", "R", "Z","U", "J", "H", "C", "F")
n <- 4
# Generate and print the strings
strings <- generate_strings(alphabet, n)
cat(strings, sep="\n")
length(strings)
write.table(strings, "strings.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

####################################################################################################### 
#######################################################################################################   
# Longest Increasing Subsequence -lgis
# Problem
# A subsequence of a permutation is a collection of elements of the permutation in the order that they appear. For example, (5, 3, 4) is a subsequence of (5, 1, 3, 4, 2).
# A subsequence is increasing if the elements of the subsequence increase, and decreasing if the elements decrease. For example, given the permutation (8, 2, 1, 6, 5, 7, 4, 3, 9), an increasing subsequence is (2, 6, 7, 9), and a decreasing subsequence is (8, 6, 5, 4, 3). You may verify that these two subsequences are as long as possible.
# Given: A positive integer n≤10000 followed by a permutation π of length n
# Return: A longest increasing subsequence of π followed by a longest decreasing subsequence of π
#lgis
# Longest Increasing Subsequence
 
# Function to find longest increasing subsequence
find_lis <- function(sequence) {
  n <- length(sequence)
  lengths <- rep(1, n)
  prev <- rep(NA, n)
  
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      if (sequence[i] > sequence[j] && lengths[i] < lengths[j] + 1) {
        lengths[i] <- lengths[j] + 1
        prev[i] <- j
      }
    }
  }
  
  # Find the maximum length and backtrack to get the subsequence
  max_length <- max(lengths)
  idx <- which.max(lengths)
  lis <- numeric(max_length)
  current_length <- max_length
  
  while (!is.na(idx)) {
    lis[current_length] <- sequence[idx]
    current_length <- current_length - 1
    idx <- prev[idx]
  }
  
  return(lis)
}

# Function to find longest decreasing subsequence
find_lds <- function(sequence) {
  # Reverse the sequence to use the same function for LIS
  reversed_sequence <- rev(sequence)
  lds_reversed <- find_lis(reversed_sequence)
  return(rev(lds_reversed))
}

# # Sample input
# n <- 5
# permutation <- c(5, 1, 4, 2, 3)


# Reading the file
file_content <- readLines("rosalind_lgis.txt")
# Extracting the first line for 'n' and the rest for 'm'
n <- as.integer(file_content[1])
# Assuming all other numbers are in the second line, split and convert to integer
m <- as.integer(unlist(strsplit(file_content[2], " ")))
# Printing the values to check
print(n)
print(m)

# Sample input
n <- n
permutation <- m

# Finding LIS and LDS
lis <- find_lis(permutation)
lds <- find_lds(permutation)

# Output
cat("LIS:", lis, "\n")
cat("LDS:", lds, "\n")

###################################################################################################### 
#######################################################################################################  

########################################################################################################### 
# Enumerating Gene Orders-perm

# Problem
# A permutation of length n is an ordering of the positive integers {1,2,…,n}. For example, π=(5,3,2,1,4)
# is a permutation of length 5
# Given: A positive integer n≤7
# Return: The total number of permutations of length n , followed by a list of all such permutations (in any order). 

# Load necessary library
library(combinat)
# Generates all permutations of numbers from 1 to 𝑛  using the permn function.
# Function to generate and print all permutations of numbers 1 to n
generate_permutations <- function(n) {
  # Generate all permutations
  permutations <- permn(n)
  
  # Print the number of permutations
  cat("Number of permutations:", length(permutations), "\n")
  
  # Print each permutation on a new line
  for (perm in permutations) {
    cat(perm, "\n")
  }
}
# The cat() function in R is used for concatenating and printing objects to the console or to a file.
# Example usage with n = 3
generate_permutations(5)

###################################################################################################### 
# Variables and Some Arithmetic-ini2
# Function to calculate the square of the hypotenuse
calculate_hypotenuse_square <- function(a, b) {
  # Calculate the square of the hypotenuse using Pythagorean theorem
  hypotenuse_square <- a^2 + b^2
  
  # Return the result
  return(hypotenuse_square)
}

# Example usage:
# a <- 3
# b <- 5
a <- 983
b <- 884

result <- calculate_hypotenuse_square(a, b)
cat(result)

###################################################################################################### 
# Function to extract and concatenate the slices-ini3
extract_slices <- function(s, a, b, c, d) {
  # Extract the first substring from index a to b (inclusive)
  slice1 <- substr(s, a + 1, b + 1)  # R uses 1-based indexing
  
  # Extract the second substring from index c to d (inclusive)
  slice2 <- substr(s, c + 1, d + 1)
  
  # Concatenate the two slices with a space in between
  result <- paste(slice1, slice2)
  
  # Return the result
  return(result)
}

# Example usage:
s <- "AMvKlwgeawqZ5nGYjdmvCxmYqjJHk8kZFBgVH3DmoS9COLZLcfs31z8LCypselus18fcyormKfmbF2Z1P8N4w89TZZuyTQr73Kquinquestriatus9CU68YQPspxxwUFm4GWqASfRZCVoJhN5K7T6pzfPbw19OFw1AQaxSPZR9Ztyc2."
a <- 56
b <- 63
c <- 98
d <- 112
result <- extract_slices(s, a, b, c, d)
cat(result)

#######################################################################################################   
# Function to calculate the sum of all odd integers between a and b-ini4
sum_odd_integers <- function(a, b) {
  # Generate a sequence from a to b
  numbers <- seq(a, b)
  
  # Filter only the odd numbers
  odd_numbers <- numbers[numbers %% 2 == 1]
  
  # Calculate the sum of the odd numbers
  sum_of_odds <- sum(odd_numbers)
  
  # Return the sum
  return(sum_of_odds)
}

# Example usage:
a <- 4184
b <- 8942
result <- sum_odd_integers(a, b)
cat(result)

###################################################################################################### 
# Function to extract even-numbered lines from a file-ini5
extract_even_lines <- function(input_file, output_file) {
  # Read the lines from the input file
  lines <- readLines(input_file)
  
  # Extract even-numbered lines (1-based indexing)
  even_lines <- lines[seq(2, length(lines), by = 2)]
  
  # Write the even-numbered lines to the output file
  writeLines(even_lines, output_file)
}

# Example usage:
input_file <- "rosalind_ini5.txt"   # Specify the input file name
output_file <- "rosalind_ini5_output.txt" # Specify the output file name
extract_even_lines(input_file, output_file)

#######################################################################################################  
# Function to count occurrences of each word in the string-ini6
# Function to count word occurrences from a text file
count_word_occurrences_from_file <- function(input_file) {
  # Read the contents of the file
  s <- readLines(input_file)  # Convert to lowercase to handle case insensitivity
  s <- paste(s, collapse = " ")  # Combine lines into a single string
  
  # Remove punctuation and split the string into words
  s <- gsub("[[:punct:]]", "", s)  # Remove punctuation
  words <- unlist(strsplit(s, "\\s+"))
  
  # Create a table of word counts
  word_counts <- table(words)
  
  # Print the words and their counts in the specified format
  for (word in names(word_counts)) {
    cat(word, word_counts[[word]], "\n")
  }
}

# Example usage:
input_file <- "rosalind_ini6.txt"  # Replace with your file path
count_word_occurrences_from_file(input_file)

