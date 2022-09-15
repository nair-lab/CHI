function [output_seq, CAI_Output, ENC_output, GC_output] = RecodeTarget_CAI(input_seq, exclusions, CAI_reference_struct, CAI_Target, GC_Target, ENC_Target, GCdiff_thresh, ENCdiff_thresh)
%Recodes an input DNA sequence to random codons that are still acceptable
%for that amino acid.  Randomizes a randomly selected position, then checks
%to see if that made the sequence better in terms of getting closet to the
%target CAI, while maintaining GC content, and ENC. The target CAI will be
%given, and CAI is calculated through the CAI_reference struct. Note that
%CAI of 0 or 1 and extremes in general cannot be practically achieved. 

%% First some housekeeping tasks

%checks to make sure the sequence inputs are type char, exclusions is a
%cell array, and the reference struct is a struct.
if ischar(input_seq) == 0
    error('Input sequence must be of type char')
    return
elseif iscellstr(exclusions) == 0
    error('Exclusions Input must be of type cell')
elseif isstruct(CAI_reference_struct) == 0
    error('CAI reference struct must be of type struct')
    return
end

%INPUT SEQUENCE MUST BE UPPER CASE!

input_seq = upper(input_seq);

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

%% This part re-codes the sequence towards objectives defined below:

%exclusions input is fine as is, algorithm will skip re-coding of any
%excluded codons.

%input_seq for example is the CFP gene, or the gene you wish to re-code. It
%will be manipulated in the loops.

%CAI_reference_struct is a struct generated using RSCUarray2strct, and is a
%struct of RSCU values to reference when calling the CAI function and
%calculating CAI.  It is what you need to optimize towards.

%CAI target is the desired CAI of the recoded sequence.

%GCdiff_thresh, ENCdiff_thresh, and targets are
%all passed in as inputs and define distances in the respective parameters
%to maintain within a range of the accepted value as defined in the
%function.

%Here we go:

%calculate the starting CAI, GC content, and ENC for our gene to be
%randomized for a gene like cfp.
RSCUstructInput = RSCUstruct('target',input_seq);%get the sequence data struct for target seq (input seq)
CAI_Input = CAI(RSCUstructInput, CAI_reference_struct); %calculate CAI
properties_input = oligoprop(input_seq); %calculate DNA properties, output is a struct with GC content and other goodies.
GC_input = properties_input.GC; %get GC content for the input sequence
ENC_input = ENC(RSCUstructInput); %calculate ENC of the input gene.

%I will use the difference between the CAI, GC content, and ENC of the
%target and reference sequence as an indicator of how far apart they are
%and a metric by which to reject sequences that aren't meeting the criteria
%of staying similar in space. All of the differencecs are passed in through
%the input.
CAIdiff_Input = abs(CAI_Target - CAI_Input);%need to know difference of input to compare to output
GCdiff = abs(GC_Target - GC_input);
ENCdiff = abs(ENC_Target - ENC_input);

%checks to make sure the input sequence is already within acceptable
%difference range for GC%, and ENC.

if GCdiff >= GCdiff_thresh
    error('Starting GC difference too far appart')
    return
elseif ENCdiff >= ENCdiff_thresh
    error('Starting ENC difference too far appart')
    return
end

%Now determine the length of the input sequence and number of codons.
len_input_seq = length(input_seq);
num_codons = len_input_seq./3;

%loop will work by generating a random sequence index, randomizing the
%codon at that position, then determining if the new sequence is closer to
%the desired criteria or not.  If it isn't then it will repeat, otherwise
%it will perform the same operation on the new sequence until the desired
%sequence parameters are met.

%All of the thresholds/targets are passed in for the parameters being held
%constant.

isloop = 1; %initialize the while loop with a logical operator
while isloop == 1
    
    %if the target CAI is acquired, then stop the while loop.
    if CAIdiff_Input <= 0.01
        isloop = 0;
    end
    
    %get the indices for a random codon
    random_codon_end_index = randi(num_codons)*3; %random number constrained by number of codons, giving the index for the last base
    random_codon_start_index = random_codon_end_index - 2; %this will be the first position
    
    %get the random codon
    random_codon = cellstr(input_seq(random_codon_start_index:random_codon_end_index));
    if ismember(random_codon, exclusions) == 1
        continue %skips over codons that are on the excluded list.
    else
        %now will randomize the random codon to another acceptable
        %synonymous codon for the given AA.
        for i = 1:length(names)%for the number of possible amino acids
            current_codon_index = AAs.(char(names(i))); %grabs the current codon in the index struct
            current_codon_index = cellstr(current_codon_index);%for matching the index codon to current sequence codon
            if ismember(current_codon_index, random_codon) == 0
                %if it matches then it will randomize to another codon,
                %otherwise keep searching
                continue
            else
                current_aa_codons = length(AAs.(char(names(i))));%number of codons for the current AA
                random = randi(current_aa_codons); %random number constrained by number of codons
                new_codon = AAs.(char(names(i)))(random); %assign new random codon
                output_seq = input_seq; %need to initialize a sequence that will become the output seq for comparison.
                output_seq(random_codon_start_index:random_codon_end_index) = new_codon;%replace the current index with the new random codon.
            end
        end
    end
    
    
    %now need to check to see if the new output sequence is more desirable
    %in terms of reducing CAI distance to target while not exceeding GC, or
    %ENC difference thresholds. I will assign new vairables to characterize
    %the output sequence, and if it meets criteria, will re-assign them as
    %the new input variables.
    
    RSCUstructOutput = RSCUstruct('target',output_seq);%get the sequence data struct for output seq
    CAI_Output = CAI(RSCUstructOutput, CAI_reference_struct); %***remove semicolon to track progress*** calculate CAI for output seq
    CAIdiff_Output = abs(CAI_Target - CAI_Output);
    %determine if CAI became closer to the target CAI
    if CAIdiff_Output >= CAIdiff_Input
        continue %restart while loop if CAI got farther away from the target.
    else
        %determine if GC% exceeds threshod
        properties_output = oligoprop(output_seq); %calculate DNA properties, output is a struct
        GC_output = properties_output.GC; %get GC content for the input sequence
        GCdiff = abs(GC_Target - GC_output);
        if GCdiff >= GCdiff_thresh == 1
            continue %restart while loop
        else
            ENC_output = ENC(RSCUstructOutput); %calculate ENC of the reference gene.
            ENCdiff = abs(ENC_Target - ENC_output);
            if ENCdiff >= ENCdiff_thresh == 1
                continue %restart while loop
            else
                %if criteria are met, make the output sequence and its
                %CAI distance the inputs for the next iteration.
                CAIdiff_Input = CAIdiff_Output;
                input_seq = output_seq;
            end
        end
    end
end
end


