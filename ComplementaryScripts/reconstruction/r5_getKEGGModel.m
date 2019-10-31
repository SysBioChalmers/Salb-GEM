% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R5: Reconstruction of KEGG Model and curation of reactions to add
% to current draft of Salb-GEM

% Updated 2019-10-21
% Cheewin Kittikunapong

% KEGG IDs
% sco = S. coelicolor
% salb = S. albus J1074

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% load FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = '../../ComplementaryData/genome/J1074_protein_XNR_sortLoci.fasta';

% load draft model for reference in later step
load('scrap/r4_draftSalb_checkGrowth');

%% Build Salb draft model from KEGG Orthologies (KO)

% As of this time, the current KEGG release supported by RAVEN is 91.0.
% KEGG Orthology-trained HMMs 
% (prok100_kegg91, prok90_kegg91,prok50_kegg91)

% We will use the Hidden Markov Models (HMMs) trained for prokaryotic
% species with a protein redunancy cut-off of 90%
dataDir = 'prok90_kegg91';

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

% reset model
KEGGModel = KEGGDraftModel;

% You can choose to export the current model to see all the reactions that
% are included from KEGG for S. albus
% exportToExcelFormat(KEGGModel, 'KEGGModel_Salb_all.xlsx');

KEGGModel = removeGenes(KEGGModel, modelSalb.genes);
find(~contains(KEGGModel.grRules, KEGGModel.genes));
KEGGModel2 = removeReactions(KEGGModel, KEGGModel.rxns(ans));

% Exporting at this stage, you can view all remaining reactions that are
% associated with genes not yet included in the draft model generated only
% from homology with the template model
% exportToExcelFormat(KEGGModel, 'KEGGModel_RxnsToAdd.xlsx');
