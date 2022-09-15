%% Example_Recode_Target_CAIorCHI_Then_Analyze.m example
%***Note, takes ~30 seconds to run***

%Example of how to re-code a given sequence to a desired CAI or CHI and
%then analyze using clustering and PCA:

%Use CFP as a starting sequence:
CFP = 'ATGGGCAAGGGCGAAGAGCTTTTTACCGGTGTTGTGCCGATTTTAGTAGAACTGGACGGAGACGTGAACGGTCATAAGTTCTCTGTTCGTGGCGAAGGAGAGGGAGATGCCACCAATGGTAAGCTGACCCTGAAGTTCATCTGTACCACCGGTAAGCTGCCCGTGCCTTGGCCGACGCTGGTCACAACGTTGACGTGGGGCGTCCAATGCTTTTCACGCTATCCAGATCACATGAAACGCCACGACTTTTTTAAAAGCGCAATGCCTGAAGGTTATGTGCAGGAACGGACTATTAGCTTCAAAGACGATGGGACGTATAAGACCCGCGCGGAAGTGAAATTTGAAGGCGATACCTTAGTTAACCGCATTGAATTAAAAGGTATCGATTTCAAAGAGGATGGGAATATCCTGGGGCACAAATTGGAATACAACTTTAATTCGCACAACGTATACATTACAGCGGATAAACAGAAAAATGGCATCAAAGCCAACTTTAAAATCCGTCATAACGTAGAGGACGGTTCCGTGCAGCTGGCTGATCATTACCAGCAGAATACTCCGATTGGCGATGGCCCCGTTCTGCTCCCGGATAATCATTACCTGTCTACACAAAGCGTTCTTAGTAAAGACCCAAACGAGAAGCGTGACCATATGGTCCTGTTGGAATTCGTCACGGCAGCGGGGATTACTCATGGCATGGATGAACTCTATAAGTAA';

%Fill out all inputs:
exclusions = {}; %we won't exclude any codons from re-coding

%For CAI reference, first load in a cell array of every codon and
%corresponding RSCU value:
load Genscript_RSCU_High

%Next create the reference struct from the values:
Genscript_RSCU_High_Struct = RSCUarray2struct(Genscript_RSCU_High,'Genscript_RSCU_High');

CAI_Target = 0.92; %The target CAI value

GC_Target = 50; %The target GC% Value

ENC_Target = 40; %The target ENC Value

GCdiff_thresh = 50; %The maximum value tolerated of the absolute value of the difference between target GC% and measured GC%...here it is high so this is not a factor.

ENCdiff_thresh = 50; %The maximum value tolerated of the absolute value of the difference between target ENC and measured ENC...here it is high so this is not a factor.

[output_seq_CAI, CAI_Output, ENC_output, GC_output] = RecodeTarget_CAI(CFP,exclusions,...
    Genscript_RSCU_High_Struct,CAI_Target,GC_Target,ENC_Target,GCdiff_thresh,ENCdiff_thresh);

CAI_Output

%Now repeat the same operation optimizing towards CHI instead of CAI:

load CHI_RSCU_High

%Next create the reference struct from the values:
CHI_RSCU_High_Struct = RSCUarray2struct(CHI_RSCU_High,'CHI_RSCU_High');

CHI_Target = 0.96; %The target CHI value

[output_seq_CHI, CHI_Output, ENC_output, GC_output] = RecodeTarget_CAI(CFP,exclusions,...
    CHI_RSCU_High_Struct,CHI_Target,GC_Target,ENC_Target,GCdiff_thresh,ENCdiff_thresh);

CHI_Output
