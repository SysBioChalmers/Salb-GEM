% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R6-B: Curation of KEGG Model reactions to add to draft of Salb-GEM

% Updated 2019-10-21
% Cheewin Kittikunapong

% KEGG IDs
% sco = S. coelicolor
% salb = S. albus J1074

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference
clear;
load('scrap/r6a_Salb_KEGGModels.mat');
load('scrap/r6a_Salb_KEGG2BiGG');
load('scrap/r5_draftSalb_curateAnnotations');

%% This approach uses RAVEN addRxnsGenesMets (directly from KEGGModel to draft)

% I did all the manual manipulation and filtering in Sublime and Excel and
% used Matlab to match indices of fields to arrange in order

% To convert KEGG Rxn/Met ID -> BiGG Rxn/Met ID:

% use BiGG database .txt file to sort reactions/metabolites so that we are
% looking at those occurring only in the cytoplasm (or extracellular space
% with manual inspection). 

% This is based on reasoning that other BiGG IDs associated with other 
% compartments would pertain to other multi-compartmental organisms 
% (e.g. possessing mitochondria).

% files used:
% (KEGGModel .xml|.xlsx):  KEGGModel.rxns|.mets (KEGG ID); can use modelStruct directly
% KEGG2BiGG.rxns           List of matches for BiGG Rxn IDs to be matched
%                          with the BiGG database file to extract info
% BiGG_rxns.csv:           List should be filtered for relevant compartments,
%                          which will most likely eliminate multiple ID matches
%                          See compartment information in rxn equations.
% (KEGG2BiGG.mets)         List of matches for BiGG met IDs. 
% BiGG_mets.csv:           Used to match IDs so that we can modify the
%                          fields mets and metNames before we add to draft
%                          model from KEGGModel structure
%
% BiGG database tab-delimited text files taken from: 
% http://bigg.ucsd.edu/data_access

% The BiGG db files have been filtered beforehand in Excel

%% Conceptual overview of procedure

% for each METABOLITE, scan the draft model OR template model to see if one
% of the multiple IDs mapped matches to an existing metabolite ID. The met
% ID already used in the model will be given preference.

% for each REACTION, scan the BiGG database file filtered for the relevant
% compartment information to see if one of the rxn IDs matches. This will
% most likely filter down multiple candidates to a manageable subset.

% Multiple matches may likely be encountered anyway but due to lower occurrence
% will allow ease in manual curation.

%% Sorting the metabolites
mets = KEGG2BiGG_Mets; %clear('KEGG2BiGG_Mets');

% use modelSalb or modelSco as reference
% modelSalb selected as met 'oh1' exists in modelSco but not in modelSalb
% none of the potential rxns to be added use 'oh1' specifically
ref = modelSalb.mets;

% Because the MNXid met outputs use universal BiGG ID's (no compartment
% information), we will have to remove compartments in the ref IDs for this
% to work:

% ref = regexprep(ref, '_e','');
ref = regexprep(ref, '_c','');
ref = unique(ref);
ref(isempty(ref)) = '';

done    = [];       % one ID matched and automatically sorted
check   = [];       % multiple ID candidates to be inspected manually
noMatch = [];       % no match at all (does not exist in model)

for i = 1:numel(mets)
    if contains(mets{i}, ';')
        manyIDs = strsplit(mets{i}, ';');
        match = ismember(manyIDs, ref);
        if numel(find(match)) == 1
            mets(i) = manyIDs(find(match));
            done    = [done i];
        elseif numel(find(match)) > 1
            check   = [check i];
        else
            check   = [check i];
            noMatch = [noMatch i];
        end
    end
end
fprintf('%1.0f met IDs have been selected automatically\n', numel(done));
fprintf('%1.0f met IDs from KEGGModel have to be manually checked:\n', numel(check));

%metsIdx         = check;
metsToCheck     = mets(check);
disp(metsToCheck);

fprintf('The following met IDs do not exist in the draft model:\n');
disp(mets(noMatch));


%% Finding the reactions the metabolites with multiple candidates participate in
% To help with deciding which met ID is most pertinent, we will sort using
% the list of candidate BiGG rxn IDs and the BiGG database:

% For each rxn that is included in the KEGG Model, we will check which
% variant of the participating metabolite's met ID occurs. This will
% provide us with recommendations of what the replacement BiGG met ID
% should be.

clear('ref');
rxns = KEGG2BiGG_Rxns; %clear('KEGG2BiGG_Rxns');
fid         = fopen(['../../ComplementaryData/reconstruction/BiGG_rxns_filtered.csv']);
loadedData  = textscan(fid,'%s %s %s %s %s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
ref.rxns        = loadedData{1};
ref.rxnNames    = loadedData{2};
ref.rxnEqns     = loadedData{3};

ref.rxnMets     = regexprep(ref.rxnEqns, ' \+ ',';');
ref.rxnMets     = regexprep(ref.rxnMets, ' <-> ',';');
% remove coefficients?? how to use \n for this??
ref.rxnMets     = regexprep(ref.rxnMets, '2.0 ',''); % goes upwards to 10.0

fprintf('In the reactions that involve the above metabolites, the following metabolite IDs are used:\n\n');

for i = 1:numel(check)
    idx = find(KEGGModel.S(check(i),:));
    typeRxn = extractBefore(rxns(idx),';');
    idx = find(ismember(ref.rxns, typeRxn));
    disp(ref.rxnMets(idx));
end
%% Manually correcting the met IDs

mets{check(1)} = 'glc__D';
mets{check(2)} = 'Glc_aD';
mets{check(3)} = '3hbcoa';
mets{check(4)} = 'f6p_c';
mets{check(5)} = '3hibutcoa';

% add compartment info back
mets = strcat(mets, '_c');

%% Sorting the reactions
% Similar procedure as above is applied to finding the reactions

done    = [];       % one ID matched and sorted automatically
check   = [];       % multiple ID hits and will require further sorting

for i = 1:numel(rxns)
     if contains(rxns{i}, ';')
        manyIDs = strsplit(rxns{i}, ';');
        match = ismember(manyIDs, ref.rxns);
        if numel(find(match)) == 1
            rxns(i) = manyIDs(find(match));
            done    = [done i];
        elseif numel(find(match)) > 1
            check   = [check i];
        end
     end 
end

fprintf('%1.0f rxns IDs have been selected automatically\n', numel(done));
fprintf('%1.0f rxns IDs from KEGGModel have to be manually checked:\n', numel(check));
disp(rxns(check));

%% Sorting the term candidates before manual curation
% The main sorting strategy is to see whether the IDs of metabolites involved in
% each of the candidate rxn IDs occur in their chemical equations. If a
% reaction's associated equation is missing at least one of the chemical
% species from the list of metabolite IDs we previously sorted, it is
% eliminated.

fprintf('These rxns in the suggested ID set have BOTH:\nrelevant compartment and metabolites contained in the KEGG Model\n\n');

for i = 1:numel(check)
    manyIDs = strsplit(rxns{check(i)}, ';');
    fprintf('For rxn %.0f (%s):\n\n', check(i), rxns{check(i)});
    for j = 1:numel(manyIDs)
        if ismember(manyIDs(j), ref.rxns)
            idx = find(ismember(ref.rxns, manyIDs{j}));
            match = strsplit(ref.rxnMets{idx}, ';');
            if numel(find(ismember(match, mets))) == numel(match)
                fprintf('%s\n\n',manyIDs{j});
            end
        end
    end
end

%% Manual correction

% If there are competing candidates suggested, sort by:
% 1) if relevant metabolites make sense; if chemical equation is identical,
% 2) then see how many and which models have adopted each ID term
% 3) 'ir' or 'i' usually designates irreversibility, which we cannot
% confirm without further information, so these will be sorted out

rxns{check(1)}		= 'SERD_L';     
rxns{check(2)}		= 'BGLA';
rxns{check(3)}		= 'ALCD2x';
rxns{check(4)}		= 'LEUTA';
rxns{check(5)}		= 'GLYCL';
rxns{check(6)}		= 'P5CRx';
rxns{check(7)}		= 'P5CR';
rxns{check(8)}		= 'GUAD';
rxns{check(9)}		= 'HACD1';
rxns{check(10)}		= 'ECOAH1';
rxns{check(11)}		= 'HPROa';
rxns{check(12)}		= 'SELNPS';
rxns{check(13)}		= 'ECOAH5';
rxns{check(14)}		= 'HACD7';      
rxns{check(15)}		= 'HACD6';
rxns{check(16)}		= 'HACD5';
rxns{check(17)}		= 'HACD3';
rxns{check(18)}		= 'HACD2';
rxns{check(19)}		= 'TREH';
rxns{check(20)}		= 'HKNDDH';
rxns{check(21)}		= 'HKNTDH';


%% Updating KEGG Model structure

% Update model structure with our sorted list of met and rxn IDs
model1 = KEGGModel;
model1.rxns = rxns;
model1.mets = mets;
% correct compartment information in model structure
model1.comps = {'c'};

% The remaining reactions with no match (in this case, they occur in the
% mitochondria or chloroplast, which for S. albus is irrelevant).
idx = find(contains(rxns, ';'));
model2 = removeReactions(model1, model1.rxns(idx)); %mets not removed
rxns(idx) = '';

model2.rxns = rxns;
    
%% Addressing duplicate met IDs
% Before we can merge our KEGG model with the draft, there are a couple of
% fixes to be made.

% First, an issue with duplicate met IDs was discovered, caused by 3
% distinct factors:

% 1) Some species are represented as both 'C' and 'G' in KEGG for compound
% and glycan, respectively.
% 2) Some species has multiple KEGG IDs due to each one describing them
% either generically or with stereochemistry context (L-amino acids)
% 3) One species (2hmc_c) has a tautomer denoted 'oxalc_c', but both
% species were erroneously assigned the met ID '2hmc_c'.

% list of duplicate IDs
dMets = {'glu__L_c' 'ala__L_c' 'sucr_c' '2hmc_c' 'lcts_c'};

% We will 'merge' the rows of the S matrix for all the reactions that each
% metabolite participates in. We will do this by summing the rows.
% The resulting sparse matrix should show that the one metabolite involves
% all the reactions from the case of the previous two distinct metabolites.

for i = 1:numel(dMets)
    i
    metIdx = find(ismember(model2.mets, dMets{i}))
    model2.S(metIdx(1),:)
    model2.S(metIdx(end),:)
    merge = model2.S(metIdx(1),:) + model2.S(metIdx(end),:);
    model2.S(metIdx(1),:) = merge;
    model2.S(metIdx(1),:)
    % make each ID to be removed 'unique' so that we can use RAVEN
    % removeMets to correctly remove the metabolite and updating the entire
    % model structure
    model2.mets(metIdx(end)) = strcat(model2.mets(metIdx(end)),'_DEL');
end

idx = find(contains(model2.mets, '_DEL'));
model2 = removeMets(model2, model2.mets(idx));

%% Addressing duplicate rxn IDs
% In the case of rxn IDs, it was demonstrated that each pair of duplicate
% reactions differed by the met IDs they used: in one case using KEGG 'G'
% metabolite IDs and another 'C' IDs.

model2.rxns{154} = 'SAGH_glycan';

% To avoid redundancy, we will remove rxns that already exist in the draft
% model. We can check what differs between their occurrences in the draft
% and the KEGG Model.
sharedRxns = intersect(model2.rxns, modelSalb.rxns);
model2 = removeReactions(model2, sharedRxns);

% This reaction is the interconversion of '2hmc_c' to its tautomer
% 'oxalc_c'. This reaction is redundant; since we merged the associated
% rows of the S matrix, the net chemical equation is null. Also, this
% reaction is the only occurrence of '2hmc_c', so we will remove the
% metabolite as well.
model2 = removeReactions(model2, '4OT',true); % remove unused metabolites

%% Adding the KEGG Model reactions to the draft model

% Before proceeding, we have to verify metNames are the same for shared met
% IDs, otherwise the addRxnsGenesMets function will not work.

mets = intersect(model2.mets,modelSalb.mets);

idx1 = getIndexes(model2, mets, 'mets');
idx2 = getIndexes(modelSalb, mets, 'mets');    
 
model2.metNames(idx1) = modelSalb.metNames(idx2);


% import reversibility information computed using eQuilibrator
fid         = fopen(['../../ComplementaryData/curation/reversibility/Salb_eQuilibrator_rev_filtered.csv']);
loadedData  = textscan(fid,'%s %s %f %f %f %f %f %f','delimiter','\t', 'HeaderLines',1); fclose(fid);
rxns        = loadedData{1};
dGm         = loadedData{5};
dGmStd      = loadedData{6};

% change bounds of any reactions that satisfy our threshold condition of
% -30 kJ/mol, referencing Gibbs free energy in physiological conditions
for i = 1:length(rxns)
    idx = find(ismember(model2.rxns, rxns(i)));
    if ~isempty(idx)
        if dGm(i) - dGmStd(i) < -30
            model2.lb(idx) = 0;
            %disp('forward only')
        elseif dGm(i) + dGmStd(i) > -30
            model2.lb(idx) = -1000;
            %disp('reversible')
        end
    end
end

modelSalb = addRxnsGenesMets(modelSalb, model2, model2.rxns, true, 'Additional reactions based on new genes from KEGG');

% There is an increase in the objective function after addition of the
% reactions from KEGG
[solution, hsSolOut] = solveLP(modelSalb, 0); -1*solution.f

%% save model to scrap folder

save('scrap/r6b_draftSalb_addKEGGRxns.mat', 'modelSalb');
