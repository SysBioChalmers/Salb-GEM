% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R2: Addition of exchange reactions and pseudoreactions from
% template model to draft model

% Updated 2019-10-21
% Cheewin Kittikunapong

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference
clear

%% Load the model structures generated from previous step, r1

% load modified Sco-GEM model
load('scrap/r1_scoGEM_newGrRules.mat');
modelSco = rmfield(modelSco, 'unconstrained');

% load draft Salb-GEM model
load('scrap/r1_draftSalb_scoHomology.mat');

%% find and add the exchange reactions from the template to draft model

% identify the exchange reactions and their indices
[ex_rxns, ex_i] = getExchangeRxns(modelSco);

% display the reaction names
% modelSco.rxnNames(ex_i)

% add the reactions
modelSalb = addRxnsGenesMets(modelSalb, modelSco, ex_rxns, 0, ...
    'Exchange reactions to test for growth in draft model', 1);

%% find and add the pseudoreactions and biomass reaction

% Incidentally, the first 15 reactions in the model.rxns field specify
% these reactions. Thus the following script is appropriate, but care
% should be taken in other cases
pseudoRxns = modelSco.rxns(1:15);

% add pseudoreactions and biomass reaction to the draft model
modelSalb = addRxnsGenesMets(modelSalb, modelSco, pseudoRxns, ...
    0, 'Pseudoreactions added to test for growth in draft model', 1);

%% adding transport reactions that do not have gene associations

% Find which rxns have no designated gene associations
noGene = find(cellfun(@isempty, modelSco.grRules));
Transport = find(getTransportRxns(modelSco));
Transport = intersect(Transport, noGene);
transportRxns = modelSco.rxns(Transport);

% add the reactions to the draft model
modelSalb = addRxnsGenesMets(modelSalb, modelSco, transportRxns, ...
    0, 'Transport reactions without gene associations to test for growth in draft model', 1);

%% Export model to .mat

save('scrap/r2_draftSalb_addScoRxns.mat', 'modelSalb');
