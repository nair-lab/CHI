function [output_seq] = RecodeDNA_singleAA_DesiredCodon(input_seq, target_AA, new_codon)
%Takes an input DNA sequence and recodes a given amino acid entirely to a
%given codon based on inputs. AA should be passed in 3 letter lower case
%terminology.  Codon should be passed in 3 letter lower case. Start codon
%should be passed in as the 1st index divisible by 3 of the codon you wish
%to start at.  E.g. if you want to leave the 1st 51 bp alone type in 52. 

%% First some housekeeping tasks

%checks to make sure the sequence inputs are type char.
if ischar(input_seq) == 0
    error('Input sequence must be of type char')
    return
elseif ischar(target_AA) == 0
    error('Input AA must be of type char')
    return
elseif ischar(new_codon) == 0
    error('New codon must be of type char')
    return
end


%%
%This section just creates a reference struct for indexing and keeping the
%format of information associated with each codon consistent.  The
%reference struct is called AAs.

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

%builds a struct with each amino acid and each codon for having
%appropriate indexing for a given amino acid and the codonsassigned to it.

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

%Creates a vector of names for indexing the AAs struct.
names = fieldnames(AAs);

%% This part re-codes the sequence for only the given AA to the new codon:

%Determine the length of the input sequence and number of codons.
len_seq = length(input_seq);
codon_index = len_seq./3; %indexing the number of codons
seq_index = 1; %indexing the first position of each codon

%Now step through the sequence, and each time the input AA is encountered,
%recode it to the desired new codon.
for i = 1:codon_index %for the number of codons in the sequence
    current_codon = cellstr(input_seq(seq_index:i*3));%gets the current codon from the input sequence as cell for matching
    %now match the current codon to a corresponding AA using the struct
    %built earlier.
    for j = 1:length(names)%for the number of possible amino acids
        current_codon_index = AAs.(char(names(j))); %grabs the current AA
        current_codon_index = cellstr(current_codon_index);%for matching the index codon to current sequence codon
        if ismember(current_codon, current_codon_index) == 1 %if the current codon is a match to the set of codons for the current jth AA.
            current_AA = (names(j)); %retrieves the current AA to match to input AA to be recoded
        else
            continue
        end
        %if there is a match then change it to the desired codon.
        if ismember(current_AA,target_AA)==1
            input_seq(seq_index:i*3) = new_codon;
            break %stop iterating through AAs. 
        end
    end
    seq_index = seq_index + 3;
end

output_seq = input_seq;
        