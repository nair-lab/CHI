%% RecodeDNA_singleAA_DesiredCodon.m example

%Example of how to re-code a given sequence and replace all instances of a
%particular amino acid to a defined codon, and then randomize the sequence:

%Use CFP as a starting sequence:
CFP = 'ATGGGCAAGGGCGAAGAGCTTTTTACCGGTGTTGTGCCGATTTTAGTAGAACTGGACGGAGACGTGAACGGTCATAAGTTCTCTGTTCGTGGCGAAGGAGAGGGAGATGCCACCAATGGTAAGCTGACCCTGAAGTTCATCTGTACCACCGGTAAGCTGCCCGTGCCTTGGCCGACGCTGGTCACAACGTTGACGTGGGGCGTCCAATGCTTTTCACGCTATCCAGATCACATGAAACGCCACGACTTTTTTAAAAGCGCAATGCCTGAAGGTTATGTGCAGGAACGGACTATTAGCTTCAAAGACGATGGGACGTATAAGACCCGCGCGGAAGTGAAATTTGAAGGCGATACCTTAGTTAACCGCATTGAATTAAAAGGTATCGATTTCAAAGAGGATGGGAATATCCTGGGGCACAAATTGGAATACAACTTTAATTCGCACAACGTATACATTACAGCGGATAAACAGAAAAATGGCATCAAAGCCAACTTTAAAATCCGTCATAACGTAGAGGACGGTTCCGTGCAGCTGGCTGATCATTACCAGCAGAATACTCCGATTGGCGATGGCCCCGTTCTGCTCCCGGATAATCATTACCTGTCTACACAAAGCGTTCTTAGTAAAGACCCAAACGAGAAGCGTGACCATATGGTCCTGTTGGAATTCGTCACGGCAGCGGGGATTACTCATGGCATGGATGAACTCTATAAGTAA';

%In a loop, re-code every amino acid to a user defined codon.  This should
%be very biased, i.e. ENC = 20.  In this example, we will create an ENC =
%20, CAI = 1 sequence by replacing every codon with the one that is most
%represented in highly expressed genes:

%make a list of all amino acids to re-code
AAs = ['ala'
'arg'
'asn'
'asp'
'cys'
'gln'
'glu'
'gly'
'his'
'ile'
'leu'
'lys'
'phe'
'pro'
'ser'
'thr'
'tyr'
'val'];

%make a list of codons to swap to for each AA (preferred CAI codons):
codons = ['GCG'
'CGT'
'AAC'
'GAC'
'TGC'
'CAG'
'GAA'
'GGT'
'CAC'
'ATC'
'CTG'
'AAA'
'TTC'
'CCG'
'TCT'
'ACC'
'TAC'
'GTT'
'TAA'];

for i = 1:length(AAs)
    curr_AA = AAs(i,:);
    curr_codon = codons(i,:);
    CFP = RecodeDNA_singleAA_DesiredCodon(CFP,curr_AA,curr_codon);
end

CFP1 = CFP;

%Now calculate new sequence parameters:

%Create an RSCU data struct of the new test sequence:
CFP1_struct = RSCUstruct('CFP',CFP1);

%Create an RSCU data struct to use as a reference, importing RSCU values
%from a set of highly expressed genes in E. coli originally downloaded from
%Genscript:

%First load in a cell array of every codon and corresponding RSCU value:
load Genscript_RSCU_High

%Next create the reference struct from the values:
Genscript_RSCU_High_Struct = RSCUarray2struct(Genscript_RSCU_High,'Genscript_RSCU_High');

%Now calculate new CAI and ENC
CAI_1 = CAI(CFP1_struct,Genscript_RSCU_High_Struct)
ENC_1 = ENC(CFP1_struct)

% CAI_1 =
% 
%      1
% 
% 
% ENC_1 =
% 
%     20

%% RecodeRand.m example
%Now randomize the sequence with RecodeRand.m, excluding all leucine codons:

%first create an array of excluded leucine codons:
exclusions = {'CTA'; 'CTC'; 'CTG'; 'CTT'; 'TTA'; 'TTG'};
CFP2 = RecodeRand(CFP1,exclusions);

%Now calculate new sequence parameters:
CFP2_struct = RSCUstruct('CFP',CFP2);

%Compare sequences, note that CTGs are unchanged:
Align = localalign(CFP1, CFP2);
Align.Alignment{1}

% 'ATGGGTAAAGGTGAAGAACTGTTCACCGGTGTTGTTCCGATCCTGGTTGAACTGGACGGTGACGTTAACGGTCACAAATTCTCTGTTCGTGGTGAAGGTGAAGGTGACGCGACCAACGGTAAACTGACCCTGAAATTCATCTGCACCACCGGTAAACTGCCGGTTCCGTGGCCGACCCTGGTTACCACCCTGACCTGGGGTGTTCAGTGCTTCTCTCGTTACCCGGACCACATGAAACGTCACGACTTCTTCAAATCTGCGATGCCGGAAGGTTACGTTCAGGAACGTACCATCTCTTTCAAAGACGACGGTACCTACAAAACCCGTGCGGAAGTTAAATTCGAAGGTGACACCCTGGTTAACCGTATCGAACTGAAAGGTATCGACTTCAAAGAAGACGGTAACATCCTGGGTCACAAACTGGAATACAACTTCAACTCTCACAACGTTTACATCACCGCGGACAAACAGAAAAACGGTATCAAAGCGAACTTCAAAATCCGTCACAACGTTGAAGACGGTTCTGTTCAGCTGGCGGACCACTACCAGCAGAACACCCCGATCGGTGACGGTCCGGTTCTGCTGCCGGACAACCACTACCTGT-CTACCCAGTCTGTTCTGTCTAAAGACCCGAACGAAAAACGTGACCACATGGTTCTGCTGGAATTCGTTACCGCGGCGGGTATCACCCACGGTATGGACGAACTGTACAAATAA'
% '|||||:||:|| |||||||||||||| ||:|||||:|| || |||||||||||||| ||||| ||:|| ||||| |||||||| ||: |:|| ||:|||||||| || || || || || |||||||| |||||||| || ||||| ||||| |||||||| || || ||||| || |||||:|| || ||||| |||||:||:||:|| ||||| || ||||| || ||||||||||| || || ||||||||:||:|| ||||||||:|| |||||||||||| |:|| || || |||||:|| |||||:|| |||||:||  |:|| ||||| ||||||||:|||||||| |||||:||| |:|||||||||||:|| || || || ||:||:|||||||| || |||||||| ||:|||||||| || |||||||| || ||||| ||||| ||||| |||||:||:||:||||| ||||||||:|| ||||||||  | ||||| ||:||:|| ||:||:|||||:|||||||| || ||||||||:||||||||||||||||| ||:||||| |||||||| |||||||| || |||: | || ||||||||:|||: |||||| || |||||||||||:|| || ||||||||||||||:|| ||:|| || || || || || ||||| |||||||||||||| ||:|||'
% 'ATGGGAAAGGGGGAAGAACTGTTCACAGGAGTTGTACCTATACTGGTTGAACTGGATGGTGATGTAAATGGTCATAAATTCTCGGTAAGAGGCGAGGGTGAAGGGGATGCTACTAATGGGAAACTGACGCTGAAATTTATTTGCACTACCGGCAAACTGCCTGTGCCCTGGCCTACACTGGTAACGACTCTGACTTGGGGAGTACAATGTTTCTCGCGCTACCCTGATCACATGAAACGGCATGATTTCTTCAAGTCAGCCATGCCGGAGGGCTACGTTCAGGAAAGAACAATTTCCTTCAAGGATGACGGAACTTACAAGACTAGAGCCGAAGTGAAATTCGAGGGTGACACGCTGGTAAACAGAATCGAACTGAAGGGCATTGATTTTAAGGAGGACGGTAATATACTGGGTCATAAGCTGGAATATAATTTCAACTCGCATAACGTCTACATTACCGCTGACAAGCAAAAGAACGGGATCAAAGCAAATTTCAAAATAAGGCACAATGTAGAGGATGGATCAGTTCAACTGGCGGATCATTACCAGCAAAACACCCCGATCGGTGATGGACCGGTGCTGCTGCCCGACAACCATTATCTGAGC-ACTCAGTCTGTACTGAGTAAAGATCCTAACGAAAAACGAGATCATATGGTTCTGCTGGAGTTTGTAACGGCCGCTGGCATTACGCACGGGATGGACGAACTGTATAAGTAA'


%Now calculate new CAI and ENC (will be different every time)
CAI_2 = CAI(CFP2_struct,Genscript_RSCU_High_Struct)
ENC_2 = ENC(CFP2_struct)

% CAI_2 =
% 
%     0.4519
% 
% 
% ENC_2 =
% 
%    51.6718

%Now randomize the sequence entirely:

%first create an empty array of excluded codons:
exclusions = {};
CFP3 = RecodeRand(CFP1,exclusions);

%Now calculate new sequence parameters:
CFP3_struct = RSCUstruct('CFP',CFP3);

%Compare sequences:

seqs = {CFP1, CFP2, CFP3};
Align = multialign(seqs);
seqalignviewer(Align)

%Now calculate new CAI and ENC
CAI_3 = CAI(CFP3_struct,Genscript_RSCU_High_Struct)
ENC_3 = ENC(CFP3_struct)

%(will be different every time), e.g.:
% CAI_3 =
% 
%     0.3225
% 
% 
% ENC_3 =
% 
%    56.6430
