import pandas as pd
import sys


genetic_code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
       "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}


def get_second_largest_base (li, cons):
    consensus = cons.replace(" ","")
    if consensus=="A":
        li[0] = -1
    elif consensus=="C":
        li[1] = -1
    elif consensus=="G":
        li[2] = -1
    elif consensus=="T":
        li[3] = -1
    if li[0]==li[1]==li[2]==li[3]:
        return "N"
    if max(li) == li[0]:
        return "A"
    elif max(li) == li[1]:
        return "C"
    elif max(li) == li[2]:
        return "G"
    else :
        return "T"


fin = sys.argv[1]
fout = sys.argv[2]

d = pd.read_csv (fin)
d['consensus_aa'] = len(d['majorsequence'])*[None]
d['secondbase_aa'] = len(d['majorsequence'])*[None]
alternates=["","",""]

for codon_beg in range(107, 10376, 3):
    #print(d['majorsequence'][codon_beg])
    codon = d['majorsequence'][codon_beg]+d['majorsequence'][codon_beg+1]+d['majorsequence'][codon_beg+2]
    codon=codon.replace(" ","")
    alternates[0] = get_second_largest_base([d['qAs'][codon_beg],d['qCs'][codon_beg],d['qGs'][codon_beg],d['qTs'][codon_beg]], d['majorsequence'][codon_beg])
    alternates[1] = get_second_largest_base([d['qAs'][codon_beg+1],d['qCs'][codon_beg+1],d['qGs'][codon_beg+1],d['qTs'][codon_beg+1]], d['majorsequence'][codon_beg+1])
    alternates[2] = get_second_largest_base([d['qAs'][codon_beg+2],d['qCs'][codon_beg+2],d['qGs'][codon_beg+2],d['qTs'][codon_beg+2]], d['majorsequence'][codon_beg+2])
    codon_alt1 = alternates[0]+d['majorsequence'][codon_beg+1]+d['majorsequence'][codon_beg+2]
    codon_alt1=codon_alt1.replace(" ","")
    codon_alt2 = d['majorsequence'][codon_beg]+alternates[1]+d['majorsequence'][codon_beg+2]
    codon_alt2=codon_alt2.replace(" ","")
    codon_alt3 = d['majorsequence'][codon_beg]+d['majorsequence'][codon_beg+1]+alternates[2]
    codon_alt3=codon_alt3.replace(" ","")
    d['consensus_aa'][codon_beg] = genetic_code[codon]
    d['consensus_aa'][codon_beg+1] = genetic_code[codon]
    d['consensus_aa'][codon_beg+2] = genetic_code[codon]
    d['secondbase_aa'][codon_beg] = genetic_code[codon_alt1]
    d['secondbase_aa'][codon_beg+1] = genetic_code[codon_alt2]
    d['secondbase_aa'][codon_beg+2] = genetic_code[codon_alt3]

d.to_csv( fout )
