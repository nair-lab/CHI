function [output_seq] = RecodeRand(input_seq, exclusions)
%Recodes an input DNA sequence to random codons that are still acceptable
%for that amino acid.

%checks to make sure the sequence input is type char
if ischar(input_seq) == 0
    error('Input sequence must be of type char')
    return
end

%added in an exclusions filter so that certain codons won't be re-coded.
%This is passed in as a cell array input variable.
if iscellstr(exclusions) == 0
    error('Exclusions Input must be of type cellstr')
    return
end
exclusions = exclusions;

%Assigns codons to each amino acid for use in struct building
ala = ["GCA"; "GCC"; "GCG"; "GCT"];
arg = ["AGA"; "AGG"; "CGA"; "CGC"; "CGG"; "CGT"];
asn = ["AAC"; "AAT"];
asp = ["GAC"; "GAT"];
cys = ["TGC"; "TGT"];
gln = ["CAA"; "CAG"];
glu = ["GAA"; "GAG"];
gly = ["GGA"; "GGC"; "GGG"; "GGT"];
his = ["CAC"; "CAT"];
ile = ["ATA"; "ATC"; "ATT"];
leu = ["CTA"; "CTC"; "CTG"; "CTT"; "TTA"; "TTG"];
lys = ["AAA"; "AAG"];
met = ["ATG"];
phe = ["TTC"; "TTT"];
pro = ["CCA"; "CCC"; "CCG"; "CCT"];
ser = ["AGC"; "AGT"; "TCA"; "TCC"; "TCG"; "TCT"];
thr = ["ACA"; "ACC"; "ACG"; "ACT"];
trp = ["TGG"];
tyr = ["TAC"; "TAT"];
val = ["GTA"; "GTC"; "GTG"; "GTT"];
stop = ["TAA"; "TAG"; "TGA"];

%builds a struct with each amino acid and each codon for calculating RSCU
%and having appropriate indexing for a given amino acid and the codons
%assigned to it.
AAs.ala = ala;
AAs.arg = arg;
AAs.asn = asn;
AAs.asp = asp;
AAs.cys = cys;
AAs.gln = gln;
AAs.glu = glu;
AAs.gly = gly;
AAs.his = his;
AAs.ile = ile;
AAs.leu = leu;
AAs.lys = lys;
AAs.met = met;
AAs.phe = phe;
AAs.pro = pro;
AAs.ser = ser;
AAs.thr = thr;
AAs.trp = trp;
AAs.tyr = tyr;
AAs.val = val;
AAs.stop = stop;

%Assigns length of the number of codons in AAs to "codons".
codons = length(fieldnames(AAs));

%Gets a vector of "names" for indexing the AAs struct.
names = fieldnames(AAs);

len_seq = length(input_seq);

codon_index = len_seq/3; %indexing the number of codons
seq_index = 1; %indexing the first position of each codon
for i = 1:codon_index %for the number of codons in the sequence
    current_codon = cellstr(input_seq(seq_index:i*3));%gets the current codon from the input sequence as cell for matching
    if ismember(current_codon, exclusions) == 1
        seq_index = seq_index + 3;
        continue %skips over codons that are on the excluded list. 
    else
        for j = 1:length(names)%for the number of possible amino acids
            current_codon_index = AAs.(char(names(j))); %grabs the current AA
            current_codon_index = cellstr(current_codon_index);%for matching the index codon to current sequence codon
            if ismember(current_codon, current_codon_index) == 0
                %if it matches then it will randomize to another codon,
                %otherwise keep searching
                continue
            else
                for k = 1:length(AAs.(char(names(j))))%for the nubmer of AA specific codons
                    current_aa_codons = length(AAs.(char(names(j))));%number of codons for the current AA
                    random = randi(current_aa_codons); %random number constrained by number of codons
                    new_codon = AAs.(char(names(j)))(random); %assign new random codon
                    input_seq(seq_index:i*3) = new_codon;%replace the current index with the new random codon. 
                    seq_index = seq_index + 3;
                    break%terminates execution of the current loop once current codon is randomized
                end
            end
        end
    end
    
    output_seq = input_seq;
    
end

