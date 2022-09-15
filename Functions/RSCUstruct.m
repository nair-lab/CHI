function [output_struct] = RSCUstruct(Seq_Name,Seq)
%Takes an input sequence, which can also be a string of concatenated
%sequences, and generates an output struct of several sequence features and
%RSCU values for each codon

%checks to make sure the sequence input is type char
if ischar(Seq) == 0
    error('Input sequence must be of type char')
    return
end

%name of the incoming sequence becomes character for putting into struct
Seq_Name = char(Seq_Name);

%creates a struct where the first field is the input sequence name
output_struct(1).name = Seq_Name;

%Adds input sequence to the output struct
output_struct(1).Seq = Seq;

%adds the sum of all codons in input sequence to the output struct for each
%AA
output_struct(1).codon_use = codoncount(output_struct(1).Seq);

%adds codon frequency to the output struct
output_struct(1).codon_freq = codonbias(output_struct(1).Seq);

%adds total codon count field to the output struct
output_struct(1).total_codons = (length(Seq)/3);

%% Make AAs struct

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

%Assigns length of the codon count output to 64 for indexing.
codons = length(fieldnames((output_struct(1).codon_use)));
%Gets a vector of names for indexing the AAs struct.
names = fieldnames(AAs);
%%
%Adds total number of codons for each amino acid to a new struct field for
%RSCU calc.
for i = 1:length(names)
    total_AA_codon_sum = 0;
    currentAA = char(names(i));
    for j = 1:length(AAs.(currentAA))%indexes the number of codons for a given AA to increment the appropriate # of times.
        current_codon = AAs.(currentAA)(j);
        current_codon_sum = output_struct(1).codon_use.(current_codon);
        total_AA_codon_sum = total_AA_codon_sum + current_codon_sum;
    end
    output_struct(1).AA_codoncount.(currentAA) = total_AA_codon_sum;
end

%Adds RSCU field to the output struct
for i = 1:length(names)
    currentAA = char(names(i));
    total_AA_codon_sum = output_struct(1).AA_codoncount.(currentAA);%references the total sum of codons for a given AA in the genome. 
    for j = 1:length(AAs.(currentAA))%indexes the number of codons for a given AA to increment the appropriate # of times.
        current_codon = AAs.(currentAA)(j);
        current_codon_sum = output_struct(1).codon_use.(current_codon);
        RSCU = current_codon_sum./(total_AA_codon_sum./length(AAs.(currentAA))); %calculates RSCU, the sum of codons over the expected sum if equally represented.
        output_struct(1).RSCU.(currentAA).(current_codon) = RSCU;
    end
end

end

