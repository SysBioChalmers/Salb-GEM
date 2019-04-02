% Cheewin Kittikunapong
% Genome-scale metabolic model reconstruction of Streptomyces 'albus' J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% KEGG IDs
% sco = S. coelicolor
% salb = S. albus J1074

% Project started 2019-03-11

%% Loading Sco (S. coelicolor) as template model
% We will build S. albus off its homology with S. coelicolor

% checked with directly opening .mat file: 
% S matrix has same values, same dimensions for all model class params

% choose to import new template from consensus scoGEM
% modelSco = importModel('template_models/ScoGEM.xml', false);

% OR

% choose to load current workspace / model draft
% load from folder 'temp'

load('template_models/temp/modelSco_NADH17b8.mat');
modelSco

% Simplify model for simulations
% model=simplifyModel(model,true,false,true,true);

% Search for genes of interest and their associated GR rules
% str = {};
% findAssociatedGrRules(modelSco, str);

% Search for  reactions with incorrectly formatted grRules and substitute
% with manually re-formatted new strings
% old_str = {}; new_str = {};
% [modelSco, indicesChanged] = changeGrRules(modelSco, old_str, new_str);

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% MATLAB/SalbGEM/seq/KEGG/...

% fastaFile = 'seq/KEGG/GCA_000359525.1_ASM35952v1_protein.faa';
% fastaFile = 'seq/GenBank/S.albus_J1074_AA.fasta';

% identical model with KEGG but with Genes labeled by loci of CDS
fastaFile = 'seq/GenBank/J1074_protein_fasta_sortLoci.txt';

%% Build Salb draft model from Sco homology

% draftModel = getModelFromHomology(models,blastStructure,getModelFor,...
%  preferredOrder,strictness,onlyGenesInModels,maxE,minLen,minIde,mapNewGenesToOld);

% cell array for reference organisms structures and IDs

refModels = modelSco;
refModelIDs = modelSco.id;

% cell array for reference FASTA files may be appended with more templates
refFastaFiles = 'template_models/Sco_all_protein.faa';

blastStructure = getBlast(organismID, fastaFile, refModelIDs, refFastaFiles);

% encountered problem standardizing grRules of template model

%draftSalb = getModelFromHomology(refModels, blastStructure, organismID);

%% Curating draftSalb

% load the current draft model built using ScoGEM homology
load('draft models/ScoHomology.mat');
disp(draftSalb);

%% Build Salb draft model from KEGG Orthologies (KO)

% KEGG Orthology-trained HMMs 
% (prok100_kegg87, prok90_kegg87,prok50_kegg87)

dataDir = 'prok100_kegg87';

% Obtain model from KEGG, generate Hidden Markov Models from dataDir
KEGGAnnotation=getKEGGModelForOrganism(organismID,'','','',false,false);

KEGGmodel = getKEGGModelForOrganism(organismID,fastaFile,dataDir, 'output',...
    false, false, false, false, 10^-30, 0.8, 0.3, -1);

% Remove bad reactions
% None removed
[KEGGModel, removedRxns]=removeBadRxns(KEGGmodel);

[fluxes, metabolite] = makeSomething(KEGGmodel, {'H+'}, true);

% Can make CO2
KEGGmodel.metNames(metabolite)
printFluxes(KEGGmodel, fluxes, false, [], [],...
    '%rxnID (%rxnName):\n\t%eqn: %flux\n');

% TMI. How many reactions? There are 98 active reactions.
countNonZero(fluxes);

balanceStructure=getElementalBalance(KEGGmodel);

%% Build Salb draft model from MetaCyc database

% Summary (referred to as S. albus J1074):
% https://biocyc.org/organism-summary?object=SALB457425

% no exchange reactions (S 1799x1248)
METAmodel = getMetaCycModelForOrganism(organismID, fastaFile);
% open exchange reactions (S 1860x1385)
% METAmodel = getMetaCycModelForOrganism(organismID, fastaFile, 1); 

[METAModel, removedRxns]=removeBadRxns(METAmodel);
