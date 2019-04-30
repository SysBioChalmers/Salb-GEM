% Cheewin Kittikunapong
% Genome-scale metabolic model reconstruction of Streptomyces 'albus' J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Script used to compare the use of different parameters including ,ax E value,
% minimum alignment length and minimum identities on BlastStructure for
% model generation by homology

% Comparison of:
%   -   # of Rxns
%   -   # of Mets
%   -   # of Genes
%   -   # of Dead End Metabolites
%   -   # of Blocked Rxns

% Project started 2019-03-11
% Script written 2019-04-16

%% Load template model
% We will build S. albus off its homology with S. coelicolor

% scoGEM grRules standardized
load('template_models/temp/modelSco_NADH17b8.mat');
modelSco

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% load FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = 'seq/GenBank/J1074_protein_fasta_sortLoci.fasta';

%% Build Salb draft model from Sco homology

% draftModel = getModelFromHomology(models,blastStructure,getModelFor,...
%  preferredOrder,strictness,onlyGenesInModels,maxE,minLen,minIde,mapNewGenesToOld);

% cell array for reference organisms structures and IDs

refModels = modelSco;
refModelIDs = modelSco.id;

% cell array for reference FASTA files may be appended with more templates
refFastaFiles = 'template_models/Sco_all_protein.faa';

% blastStructure has been exported manually into an Excel file for closer
% inspection
blastStructure = getBlast(organismID, fastaFile, refModelIDs, refFastaFiles);


%% Use blastStructure to generate a range of models with different parameters

% default parameters

% draftModel = getModelFromHomology(models,blastStructure,getModelFor,...
%  preferredOrder,strictness,onlyGenesInModels,maxE,minLen,minIde,mapNewGenesToOld);

maxE = 10^-30;
minLen = 200;
minIde = 40; %(in percent%)

%% different maximum E values

for i = 1:11
    %maxE = 10^(-36 + i);
    %model_E(i) = getModelFromHomology(refModels, blastStructure, organismID, ...
    %    {}, 1, 0, maxE);
    outputMets = detectDeadEnds(model_E(i));
    model_E(i).deadEnds = length(outputMets);
    DeadEnds = model_E(i).mets(outputMets);
    [rxnList, rxnFormulaList] = findRxnsFromMets(model_E(i), DeadEnds);
    model_E(i).blockedRxns = length(rxnList);
end

%% different minimum alignment lengths

for i = 1:5
%    minLen = 80+20*i;
%    model_L(i) = getModelFromHomology(refModels, blastStructure, organismID, ...
%        {}, 1, 0, maxE, minLen);
    outputMets = detectDeadEnds(model_L(i));
    model_L(i).deadEnds = length(outputMets);
    DeadEnds = model_L(i).mets(outputMets);
    [rxnList, rxnFormulaList] = findRxnsFromMets(model_L(i), DeadEnds);
    model_L(i).blockedRxns = length(rxnList);
end 

%% different minimum identities

for i = 1:5
%    minIde = 20*i;
%    model_I(i) = getModelFromHomology(refModels, blastStructure, organismID, ...
%        {}, 1, 0, maxE, minLen, minIde);
    outputMets = detectDeadEnds(model_I(i));
    model_I(i).deadEnds = length(outputMets);
    DeadEnds = model_I(i).mets(outputMets);
    [rxnList, rxnFormulaList] = findRxnsFromMets(model_I(i), DeadEnds);
    model_I(i).blockedRxns = length(rxnList);
end 
