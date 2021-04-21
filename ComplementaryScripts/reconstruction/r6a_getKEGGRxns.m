% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R6-A: Reconstruction of KEGG Model to for new gene-associated reactions

% Updated 2019-10-21
% Cheewin Kittikunapong

% KEGG IDs
% sco = S. coelicolor
% salb = S. albus J1074

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference
clear

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% load FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = '../../ComplementaryData/genome/J1074_protein_XNR_sortLoci.fasta';

% load draft model for reference in later step
load('scrap/r5_draftSalb_curateAnnotations');

%% Build Salb draft model from KEGG Orthologies (KO)

% As of this time, the current KEGG release supported by RAVEN is 91.0.
% KEGG Orthology-trained HMMs include:
% [prok100_kegg94, prok90_kegg94, prok50_kegg94]

% We will use the Hidden Markov Models (HMMs) trained for prokaryotic
% species with a protein redundancy cut-off of 90%
dataDir = 'prok90_kegg94';

% If this is your first time running getKEGGModelForOrganism on your
% device, it is normal for the run time to be longer as RAVEN is generating
% anew the MATLAB .mat files used in reconstruction

% obtain model from KEGG based on annotations alone 
KEGGAnnotation=getKEGGModelForOrganism(organismID);

% generate model from Hidden Markov Models
KEGGHomology = getKEGGModelForOrganism(organismID,fastaFile,dataDir, 'output',...
    false, false, false, false, 10^-50, 0.3, 1, 0);

% merge the models
KEGGDraftModel=mergeModels({KEGGAnnotation KEGGHomology});
KEGGDraftModel=expandModel(KEGGDraftModel);
KEGGDraftModel=contractModel(KEGGDraftModel);

%% Removing Salb-GEM genes and associated rxns from KEGG model

% At this step, we will remove any reactions from the KEGG model that have
% gene associations with genes already existent in Salb-GEM. This is
% because we are looking for any new reactions that are introduced with new
% genes that were not added based on our previous homology reconstruction.

% NOTE: This step used to be step R5 preceding curation of gene
% associations and relevant annotations. However, it was acknowledged that
% curation of grRules for rxns added in steps R3 and R4 will affect the
% following steps' output.

% reset model
KEGGModel = KEGGDraftModel;

% You can choose to export the current model to see all the reactions that
% are included from KEGG for S. albus
% exportToExcelFormat(KEGGModel, 'scrap/KEGGModel_Salb_all.xlsx');

KEGGModel = removeGenes(KEGGModel, modelSalb.genes, true, true);

% Exporting at this stage, you can view all remaining reactions that are
% associated with genes not yet included in the draft model generated only
% from homology with the template model
% exportToExcelFormat(KEGGModel, 'scrap/KEGGModel_RxnsToAdd.xlsx');

%% Before proceeding to the mapping section, check for and/or obtain MetaNetX

% Using the MetaNetX database, we will convert the KEGG metabolite and
% reaction IDs in the current model to their corresponding BiGG IDs that
% are used in common RAVEN and COBRA GEMs

% Depending on the status of the RAVEN toolbox at time of running this
% script, you may or may not have to change to the 'feat/addMetaNetX'
% branch to call the following function:

% Note that by switching the branch you will lose the 'KEGGModel.mat' file
% generated from previous steps. If you seek to run this script again, it
% is recommended to manually copy the file to a separate directory before
% switching to expedite the KEGG model generation in later runs.

% UPDATE: A recent temporary solution implemented here is to automatically 
% download and add the relevant MetaNetX dependencies you will need to run 
% the next set of scripts without having to switch branches manually.

dir = what('RAVEN/external/');
dirBack = what('Salb-GEM/ComplementaryScripts/reconstruction');

if ~exist('RAVEN/external/MetaNetX', 'dir')
    fprintf('Downloading MetaNetX folder from RAVEN repository...\t');
    cd(dir.path);
    websave('MetaNetX.zip', 'https://chalmersuniversity.box.com/shared/static/90mjq7ywbpswuo3slyfyiuhzk3xfj0tm');
    unzip('MetaNetX.zip'); delete MetaNetX.zip; addpath('MetaNetX');
    cd(dirBack.path);
    fprintf('Done!\n');
end

%% Mapping as many reactions and metabolite IDs as possible

KEGG2BiGG_Rxns = mapIDsViaMNXref('rxns',KEGGModel.rxns,'KEGG','BiGG');
KEGG2BiGG_Mets = mapIDsViaMNXref('mets',KEGGModel.mets,'KEGG','BiGG');

% Some metabolites and metabolites do not have corresponding matches between 
% KEGG/BiGG. As a result, the resulting cells will be marked for removal in
% the following steps:

metIndex = find(cellfun(@isempty,KEGG2BiGG_Mets));
rxnIndex = find(cellfun(@isempty,KEGG2BiGG_Rxns));
% Note that there may be several matches for metabolites and reactions 
% (e.g. water either as H2O or OH)

KEGGModel = removeReactions(KEGGModel, KEGGModel.rxns(rxnIndex), true, true);

%% Mapping continued
% This step is repeated from above to update indices to our reduced model:

KEGG2BiGG_Rxns = mapIDsViaMNXref('rxns',KEGGModel.rxns,'KEGG','BiGG');
KEGG2BiGG_Mets = mapIDsViaMNXref('mets',KEGGModel.mets,'KEGG','BiGG');

metIndex = find(cellfun(@isempty,KEGG2BiGG_Mets));
KEGGModel = removeMets(KEGGModel, KEGGModel.mets(metIndex), true, true);

% Consequently, we will export any reactions without matches with 'X' 
rxnIndex = find(cellfun(@isempty,KEGG2BiGG_Rxns));
%KEGGModel.rxns(ans) = strrep(KEGGModel.rxns(ans), 'R', 'X');

% If you would like to view in Excel with new reaction equations based on
% BiGG ideas, execute the following:
% KEGGModel.metNames = KEGG2BiGG_Mets;
% exportToExcelFormat(KEGGModel, 'KEGGModel.xlsx');

%% save workspace

% save KEGG Models to not have to generate again
save('scrap/r6a_Salb_KEGGModels.mat','KEGGAnnotation','KEGGDraftModel', 'KEGGHomology', 'KEGGModel');

% save MNX mapping
save('scrap/r6a_Salb_KEGG2BiGG.mat','KEGG2BiGG_Mets','KEGG2BiGG_Rxns');
