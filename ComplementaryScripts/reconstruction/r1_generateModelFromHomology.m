% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R1: Draft model generation from homology with template model

% Updated 2019-10-21
% Cheewin Kittikunapong

%% Importing the template model and verifying grRules are standardized

% We will build S. albus off its homology with S. coelicolor referring to
% the consensus genome-scale metabolic model (GEM) referred to as Sco-GEM
% The model under development is based on the latest iteration iKS1317 by
% Kumelj et al. 2018 and also draws from models iAA1259, iMK1208 and Sco4,
% the latter which was generated using our RAVEN 2.0 toolbox by drawing 
%information from MetaCyc and KEGG databases

% In this script, we will use the current version of Sco-GEM 1.2.0 
% (tag: v.1.2.0 on https://github.com/SysBioChalmers/Sco-GEM)
% (hash: 616dffb, viewed using git log or, if set, git graph)

modelSco = importModel...
    ('../../ComplementaryData/reconstruction/templateModels/Sco-GEM.xml', false);

modelSco.id = 'Sco-GEM';

% We have to fix some grRules that are not written suitably for use with
% the RAVEN function getModelFromHomology.
% If not fixed, the function will return a series of warnings for reactions
% that have to be adjusted for their grRules

% Import list of reactions with grRules to fix for Sco-GEM 1.2.0
newGrRules = readtable...
    ('../../ComplementaryData/reconstruction/templateModels/Sco-GEM_newGrRules.csv');
newGrRules = table2cell(newGrRules);

% Used my own script to correct the grRules. Note that RAVEN has function
% changeGrRules, which utilizes different inputs
for i = 1:length(newGrRules)
   modelSco = modifyGrRules(modelSco, newGrRules{i,2}, newGrRules{i,3});
end

% You can verify the correct reactions have been modified for grRules by
% exporting the model to an Excel file
% exportToExcelFormat(modelSco, 'scrap/Sco-GEM_newGrRules.xlsx')

%% %% Specify parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb
organismID = 'salb';

% specify FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = '../../ComplementaryData/genome/J1074_protein_XNR_sortLoci.fasta';

% cell array for reference organisms structures and IDs

refModels = modelSco;
refModelIDs = modelSco.id;

% cell array for reference FASTA files may be appended with more templates
refFastaFiles = '../../ComplementaryData/genome/Sco_all_protein.faa';

%% Build Salb draft model from Sco homology

% generate the blastStructure from bidirectional BLASTP between S. albus
% and S. coelicolor FASTA files

blastStructure = getBlast(organismID, fastaFile, refModelIDs, refFastaFiles);

% You can also generate the blastStructure  with FASTA files from different
% assemblies; however, we will use XNR notated assembly, which is based on
% Zaburannyi et al., 2014 and is also the notation in KEGG database

% We will generate the model using non-default parameters in the function:
% E-value threshold of 10^-50, minimum alignment length of 90 and
% identity of 40

% These parameters were chosen to avoid excessively non-selective addition
% of reactions in later steps to restore functional fructose metabolism (r4)

modelSalb = getModelFromHomology(refModels, blastStructure, organismID, {},...
    1, 0, 10^-50, 90, 40, true);
modelSalb.id = 'Salb-GEM';

%% Removing an arbritary metabolite term from Sco-GEM template model

% Before proceeding, we must remove a pseudometabolite that was used
% specifically for Sco-GEM that was added to Salb-GEM

% Switching back acetyl-CoA metabolite terms
find(contains(modelSalb.mets, 'accoa_res_c'))
modelSalb.mets(ans)

% reactions that use this pseudometabolite
% IPPS & MMSAD3

modelSalb = changeRxns(modelSalb, 'IPPS', ...
    '3mob_c + accoa_res_c + h2o_c <=> 3c3hmp_c + coa_c + h_c', 1);

modelSalb = changeRxns(modelSalb, 'MMSAD3', ...
    'coa_c + msa_c + nad_c => accoa_res_c + co2_c + nadh_c', 1);

modelSalb = removeMets(modelSalb, 'accoa_res_c');


%% Export model to .xml and .mat

% make a 'scrap' directory that is ignored by git
mkdir('scrap');

% save the model structure for Sco-GEM and Salb-GEM for later steps
save('scrap/r1_scoGEM_newGrRules.mat', 'modelSco')
save('scrap/r1_draftSalb_scoHomology.mat', 'modelSalb');

% You can save the blastStructure as well to avoid spending time
% regenerating the data
% save('scrap/blastStructure.mat');

