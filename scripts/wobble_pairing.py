#!/usr/bin/env python
"""
List of codons: [(codon_watson, codon_wobble)]
Wobble pairing rule follows the parsimony method:
1) codons with cognate tRNAs are assigned canonical decoding
2) codons without cognate tRNA are assigned foloowing restricted wobbling
G:U, A:C
3) extended paring (A:A, U:G) are assumed for codons that remain unassigned

Wobble pairs of Saccharomyces cerevisiae from tRNA database
"""
# anti-codon base : codon base
# G:C Watson-Crick 
# G:U Wobble
GCU_pairs = [ ("GGC", "GGT"), # Gly
              ("AGC", "AGT"), # Ser
              ("CTC", "CTT"), # Leu
              ("TTC", "TTT"), # Phe
              ("AAC", "AAT"), # Asn
              ("GAC", "GAT"), # Asp
              ("CAC", "CAT"), # His
              ("TAC", "TAT"), # Tyr
              ("TGC", "TGT") # Cys
          ]

# anti-codon base : codon base
# I:C Wobble
# I:U Wobble (less favorable geometry and two hydrogen bonds)
ICU_pairs = [ ("GCC", "GCT"), # Ala
              ("CCC", "CCT"), # Pro
              ("ACC", "ACT"), # Thr
              ("GTC", "GTT"), # Val
              ("TCC", "TCT"), # Ser
              ("CGC", "CGT"), # Arg
              ("ATC", "ATT") # Ile note: ATC has 1 canonical tRNA
          ]

# anti-codon base : codon base
# U:A Watson-Crick
# U:G Wobble
UAG_pairs = [ ("GCA","GCG"), # Ala
              ("CCA","CCG"), # Pro
              ("CTA","CTG") # Leu
          ]
