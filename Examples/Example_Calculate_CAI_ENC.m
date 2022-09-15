%Example of how to calculate CAI and ENC from a given sequence:

%Use CFP as an example:
CFP = 'ATGGGCAAGGGCGAAGAGCTTTTTACCGGTGTTGTGCCGATTTTAGTAGAACTGGACGGAGACGTGAACGGTCATAAGTTCTCTGTTCGTGGCGAAGGAGAGGGAGATGCCACCAATGGTAAGCTGACCCTGAAGTTCATCTGTACCACCGGTAAGCTGCCCGTGCCTTGGCCGACGCTGGTCACAACGTTGACGTGGGGCGTCCAATGCTTTTCACGCTATCCAGATCACATGAAACGCCACGACTTTTTTAAAAGCGCAATGCCTGAAGGTTATGTGCAGGAACGGACTATTAGCTTCAAAGACGATGGGACGTATAAGACCCGCGCGGAAGTGAAATTTGAAGGCGATACCTTAGTTAACCGCATTGAATTAAAAGGTATCGATTTCAAAGAGGATGGGAATATCCTGGGGCACAAATTGGAATACAACTTTAATTCGCACAACGTATACATTACAGCGGATAAACAGAAAAATGGCATCAAAGCCAACTTTAAAATCCGTCATAACGTAGAGGACGGTTCCGTGCAGCTGGCTGATCATTACCAGCAGAATACTCCGATTGGCGATGGCCCCGTTCTGCTCCCGGATAATCATTACCTGTCTACACAAAGCGTTCTTAGTAAAGACCCAAACGAGAAGCGTGACCATATGGTCCTGTTGGAATTCGTCACGGCAGCGGGGATTACTCATGGCATGGATGAACTCTATAAGTAA';

%Create an RSCU data struct of the new test sequence:
CFP_struct = RSCUstruct('CFP',CFP);

%Create an RSCU data struct to use as a reference, importing RSCU values
%from a set of highly expressed genes in E. coli originally downloaded from
%Genscript:

%First load in a cell array of every codon and corresponding RSCU value:
load Genscript_RSCU_High

%Next create the reference struct from the values:
Genscript_RSCU_High_Struct = RSCUarray2struct(Genscript_RSCU_High,'Genscript_RSCU_High');

%Now calculate CAI and ENC
CAI = CAI(CFP_struct,Genscript_RSCU_High_Struct)
ENC = ENC(CFP_struct)