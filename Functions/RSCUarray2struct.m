function [output_struct] = RSCUarray2struct(input_array, name)
%Converts an array of codons and RSCU values to a struct for input into
%other functions. The struct matches the format given by the RSCUstruct
%function. The main purpose is to easily import a new array of codons and
%their RSCU values from a gene set to calculate CAI.

%input array should be in the form of a cell array e.g.

%codon1 RSCU1;
%codon2 RSCU2 ...etc

%CODONS SHOULD BE UPPER CASE!!!

%Name is just a name of the array, can be a string or char. 

%%
%Starts off with some housekeeping tasks

%checks to make sure the sequence input is type array
if iscell(input_array) == 0
    error('Input array must be of type cell')
    return
end

%name of the incoming array becomes character for putting into struct
name = char(name);

%creates a struct where the first field is the input name
output_struct(1).name = name;

%%

%This section just creates a reference struct for indexing and keeping the
%format of information associated with each codon consistent.  The
%reference struct is called AAs.  There is also the variable "codons",
%which provides the length of iterations to get through all the codons, and
%names, which provides a vector of all the names of each AA.

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

%Gets a vector of names for indexing the AAs struct.
names = fieldnames(AAs);

%%

%Next we convert the array to a struct in the format where an AA can be
%called, then see the RSCU values associated with each in a format that
%matches the order of codons in the AA field.

for i = 1:length(names) %traverse each named field
    currentAA = char(names(i));%gets the ith amino acid
    for j = 1:length(AAs.(currentAA))%indexes the number of codons for a given AA to increment the appropriate # of times.
        current_codon = AAs.(currentAA)(j);%gets the jth codon for the ith amino acid
        for k = 1:length(input_array) %traverse each line of the input array
            current_array_codon = input_array(k,1); %for the kth codon in the array
            if ismember(current_codon, current_array_codon) == 1 %if the codon matches the jth codon
                output_struct.RSCU.(currentAA).(current_codon) = cell2mat(input_array(k,2)); %assign RSCU value to that codon field and convert to double
            end
        end
    end
end

%the output struct is the final delivered variable.
