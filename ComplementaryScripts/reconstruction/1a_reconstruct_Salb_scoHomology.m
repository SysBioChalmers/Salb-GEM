% Cheewin Kittikunapong
% Genome-scale metabolic model reconstruction of Streptomyces 'albus' J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% KEGG IDs
% sco = S. coelicolor
% salb = S. albus J1074

% Project started 2019-03-11
% Script written 2019-04-15

%% Loading Sco (S. coelicolor) as template model
% We will build S. albus off its homology with S. coelicolor

% choose to import new template from consensus scoGEM
% modelSco = importModel('template_models/ScoGEM.xml', false);

% OR

% choose to load template model with standardized grRules
% NOTE: 3 rxns, 1 met lost via this option

% scoGEM grRules standardized
modelSco = importModel...
    ('../../ComplementaryData/reconstruction/templateModels/scoGEM_stdGrRules.xml',false);

% If any, search for  reactions with incorrectly formatted grRules and substitute
% with manually re-formatted new strings
% old_str = {}; new_str = {};
% [modelSco, indicesChanged] = changeGrRules(modelSco, old_str, new_str);

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% load FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = '../../ComplementaryData/genome/J1074_protein_fasta_sortLoci.txt';

%% Build Salb draft model from Sco homology

% draftModel = getModelFromHomology(models,blastStructure,getModelFor,...
%  preferredOrder,strictness,onlyGenesInModels,maxE,minLen,minIde,mapNewGenesToOld);

% cell array for reference organisms structures and IDs

refModels = modelSco;
refModelIDs = modelSco.id;

% cell array for reference FASTA files may be appended with more templates
refFastaFiles = '../../ComplementaryData/genome/Sco_all_protein.faa';

% blastStructure has been exported manually into an Excel file for closer
% inspection
blastStructure = getBlast(organismID, fastaFile, refModelIDs, refFastaFiles);

% generate a draft model based on homology with default parameters 
% (E value, alignment length and identity)

draftSalb = getModelFromHomology(refModels, blastStructure, organismID);

%% Save variables for later steps

mkdir('scrap');
save('scrap/templateSco.mat', 'modelSco');
save('scrap/draft_ScoHomology.mat', 'draftSalb');