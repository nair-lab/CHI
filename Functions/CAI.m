function [CAI_Out] = CAI(target_RSCUstruct,reference_RSCUstruct)
%%
%Calculates the CAI of a sequence in reference to a given struct of RSCU
%values for a reference set of genes. Requires first that you generate
%RSCUstruct for a reference set and sequence of interest. The "RSCUstruct"
%function does this. Function inputs are 2 structs generated by the
%function RSCUstruct for the sequence in question and the reference
%sequence (can be concatenated sequences of ORFs for instance if
%referencing a whole genome or could be a subset of highly expressed
%genes). Alternatively, if you already have RSCU values for each codon, the
%reference struct can be generated by the RSCUarray2struct function, which
%converts an array of codons and RSCU values to a struct that can be used
%as an input in the CAI function.

%***Note this version does not include stop codons when calculating CAI***

%% Input Checks
%checks to make sure the inputs are structs

if isstruct(target_RSCUstruct) == 0
    error('Target input must be of type struct')
    return
elseif isstruct(reference_RSCUstruct) == 0
    error('Reference input must be of type struct')
    return
end

%% Reference struct builder for amino acid/codon table. 
%This section just creates a reference struct for indexing each amino
%acid/codon, and keeping the format of information associated with each
%codon consistent. The output struct to use as a reference is called "AAs".
%The RSCUstruct function creates the same formatted data struct, so input
%structs are all indexed the same way as the AA struct built below.

%Assigns codons to each amino acid for use in struct building, note we do
%not include the stop codon here for CAI calculation!
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
%note we do not include the stop codon here for CAI calculation!

%builds a struct with each amino acid and each set of codons for
%calculating RSCU and having appropriate indexing for a given amino acid
%and the codons assigned to it. NO STOP CODON HERE!
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

%Make a list of names for indexing the AAs struct.
names = fieldnames(AAs);

%%
%Now calculating the maximum RSCU value for each AA in the reference struct
%to use in CAI calculation.  The output is a new struct used by the
%function called "reference_RSCUstruct_MAX" that looks the same as the
%input "reference_RSCUstruct ", only the AA field holds the maximum RSCU
%value for that AAs most preferred codon instead of RSCU values for each
%codon.

for i = 1:length(names) %traverse each AA
    currentAA = char(names(i));%gets the ith amino acid
    RSCUvalues = {}; %initialize a cell array to store max RSCU values
    count = 1;%for indexing RSCUvalues
    for j = 1:length(AAs.(currentAA))%indexes the number of codons for a given AA to increment the appropriate # of times.
        current_codon = AAs.(currentAA)(j);%gets the jth codon for the ith amino acid
        current_codon_RSCU = reference_RSCUstruct.RSCU.(currentAA).(current_codon);%gets the current RSCU
        RSCUvalues{count} = current_codon_RSCU; %add to growing list of the ith AAs RSCU values
        count = count+1;
    end
    %after the J loop finishes, determine the maximum RSCU value and add to
    %new struct:
    [M,~] = max(cell2mat(RSCUvalues)); %get the maximum RSCU value
    reference_RSCUstruct_MAX.(currentAA)= M; %add the max RSCU value to the growing new struct
end


%%
%Now the function calculates the CAI.  This is done by iterating through
%each codon, and calculating the individual codon relative adaptiveness
%(weights) by taking the ratio of the codon's RSCU value in the REFERENCE
%RSCU struct over the max codon RSCU for a given AA.  CAI is simply the
%geometric mean of all the weights, where each respective codon weight is
%multiplied every time that codon appears in a sequence, then the nth root
%(where n is the total number of codons) is taken.  Given the input
%structs, we can just multiply each calculated weight by the number of
%times that codon appears in the sequence using the codonuse input struct
%field, then multiply all of the products together and take the nth root
%(or use the geomean function).

weight_vec = []; %initialize a vector of weights.
for i = 1:length(names) %traverse each named field
    currentAA = char(names(i));%gets the ith amino acid
    for j = 1:length(AAs.(currentAA))%indexes the number of codons for a given AA to increment the appropriate # of times.
        current_codon = AAs.(currentAA)(j);%gets the jth codon for the ith amino acid
        %now grab both RSCU values from the input and reference struct
        curr_RSCU_ref = reference_RSCUstruct.RSCU.(currentAA).(current_codon);
        max_RSCU_ref = reference_RSCUstruct_MAX.(currentAA); %maximum RSCU value for a given AA
      
        curr_weight = curr_RSCU_ref./max_RSCU_ref; %calculate the weight for the current codon
        curr_codoncount = target_RSCUstruct.codon_use.(current_codon); %get the number of instances the current codon appears in the target sequence
        %now add the codon weight however many number of times it appears
        %in the sequence.  I will make a new vector each time and
        %concatenate it to the growing weight vector.
        curr_weights = [];
        curr_weights(1:curr_codoncount) = curr_weight; %just a vector of the weight n number of times where n is the codon count
        %Now concatenate to the growing weight vector:
        weight_vec = [weight_vec,curr_weights];
    end
end

%Now calculate the CAI by taking the geometric mean of the weight vector.
CAI_Out = geomean(weight_vec);
