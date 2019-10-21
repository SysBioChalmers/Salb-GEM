% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R4: Checking growth of draft model on known substrates

% Updated 2019-10-21
% Cheewin Kittikunapong

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference
clear

%% Checking production of biomass/growth phenotype with known substrates

% In this step, we will verify S. albus growth on known substrates

% What nutrients can S. albus grow on: 
% carbon sources: glucose, xylose, mannitol, fructose
% does not grow in: isoinositol, rhamnose, sucrose

% (Source: Felipe Lombo group, UniOvi, SynBio4Flav, email correspondence)

% Before running this step, please refer to supplementary script S4
% (s4_checkGrowth_scoGEM.m) to see how certain reactions have been identified 
% to support growth in various conditions in the draft model

%% load models from previous steps

% load model from previous gap-filling step
load('scrap/r3_draftSalb_gapFill.mat');
load('scrap/r1_scoGEM_newGrRules.mat');
modelSco = rmfield(modelSco, 'unconstrained');
ex_rxns = getExchangeRxns(modelSco);
modelSco = removeReactions(modelSco, 'ATPM');
modelSco = removeReactions(modelSco, ex_rxns);

%% Change exchange reaction bounds

% get exchange reaction indices (make sure model structure does not have
% field 'unconstrained', which will cause function to return empty array
modelSalb = rmfield(modelSalb, 'unconstrained');
EX_Rxns = getExchangeRxns(modelSalb);

%% Switching carbon sources

% We will be testing whether the draft model as it is can grow on sources
% besides glucose

% First, set uptake of glucose to null.
modelSalb = setParam(modelSalb, 'lb', 'EX_glc__D_e', 0);
reset = modelSalb;

%% Testing growth on fructose as sole carbon source

modelSalb = setParam(modelSalb, 'lb', 'EX_fru_e', -0.8);
[solution, hsSolOut] = solveLP(modelSalb, 0);

% As shown, the model does not grow on fructose only

% Also, for an unknown reason a subsequent gap-filling step does not
% suggest any additional reactions to satisfy our specified conditions
[ConnectedRxns, cannotConnect, addedRxns, newModel, exitFlag] = fillGaps...
    (modelSalb, modelSco, false, true);

%% Adding reactions to support fructose utilization

% reactions through which carried flux in ScoGEM when FBA with fructose only;
load('scrap/s4_missingRxns_fruPTS1');

% This is a quick fix method to ensure the model can grow on fructose until
% we have more information
modelSalb = addRxnsGenesMets(modelSalb, modelSco, missingRxns_fruPTS1, false, ...
    'Included for functional fructose metabolism');

% verify functional model
[solution, hsSolOut] = solveLP(modelSalb, 0); -1*solution.f

% see flux distribution in the reactions we have added
rxnIndexes = getIndexes(modelSalb, missingRxns_fruPTS1,'rxns');
fluxes = solution.x(rxnIndexes); noFlux = find(~fluxes);
rxnIndexes = rxnIndexes(noFlux);

% see added rxns that do not carry flux in the draft model
% modelSalb.rxns(rxnIndexes)

disp('The following added reactions do not carry flux in Salb-GEM');
disp(modelSalb.rxnNames(rxnIndexes));

modelSalb = setParam(modelSalb, 'lb', 'EX_fru_e', 0);

%% Testing growth on mannitol on current model 

% mannitol only
modelSalb = setParam(modelSalb, 'lb', 'EX_mnl_e', -0.8);
[solution, hsSolOut] = solveLP(modelSalb, 0); -1*solution.f

%% Adding reactions to support mannitol utilization

% reactions through which carried flux in ScoGEM when FBA with fructose only;
load('scrap/s4_missingRxns_mnl');

% By doing the same procedure as with fructose previously, we have added
% five more reactions
modelSalb = addRxnsGenesMets(modelSalb, modelSco, missingRxns_mnl, false,...
    'Included for functional mannitol metabolism');

% Once again, verify a functional model
[solution, hsSolOut] = solveLP(modelSalb, 0); -1*solution.f

% see flux distribution in the subset of reactions we have added for
% mannitol metabolism
rxnIndexes = getIndexes(modelSalb, missingRxns_mnl,'rxns');
fluxes = solution.x(rxnIndexes); noFlux = find(~fluxes);
rxnIndexes = rxnIndexes(noFlux);

% see added rxns that do not carry flux in the draft model
% modelSalb.rxns(rxnIndexes)

disp('The following added reactions do not carry flux in Salb-GEM');
disp(modelSalb.rxnNames(rxnIndexes));

modelSalb = setParam(modelSalb, 'lb', 'EX_mnl_e', 0);

%% Adding reactions to support xylose utilization

modelSalb = setParam(modelSalb, 'lb', 'EX_xyl__D_e', -0.8);
[solution, hsSolOut] = solveLP(modelSalb, 0); -1*solution.f

% Apparently no additional rxns are needed as the model as iscan grow on
% xylose.

% In fact, our starting model that was gap-filled for growth on glucose can
% already support growth on xylose; however upon addition of the reactions
% for fructose and mannitol utilization, the flux through biomass was
% improved

% For this reason, no additional reactions are incorporated.

modelSalb = setParam(modelSalb, 'lb', 'EX_xyl__D_e', 0);

%% Verify no growth for known true negatives

modelSalb= setParam(modelSalb, 'lb', 'EX_inost_e', -0.8);
modelSalb = setParam(modelSalb, 'lb', 'EX_rmn_e', -0.8);
modelSalb = setParam(modelSalb, 'lb', 'EX_sucr_e', -0.8);
[solution, hsSolOut] = solveLP(modelSalb, 0); -1*solution.f

modelSalb= setParam(modelSalb, 'lb', 'EX_inost_e', 0);
modelSalb = setParam(modelSalb, 'lb', 'EX_rmn_e', 0);
modelSalb = setParam(modelSalb, 'lb', 'EX_sucr_e', 0);

modelSalb = setParam(modelSalb, 'lb', 'EX_glc__D_e', -0.8);

% As demonstrated the current model does not grow on all of the above 
% combined even if we have added rxns in for fructose and mannitol
% metabolism. Because we lack further mechanistic information, we will 
% leave the model as is.

% REMEMBER: these reactions have been added to support fructose and
% mannitol utilization BECAUSE fillGaps on RAVEN does not help fill 
% in these rxns. Adding in ALL the rxns carrying flux for them seems
% unrealistic with no way to verify, but the reactions have been marked if
% we do have to come back to fix them.

%% save model

save('scrap/r4_draftSalb_checkGrowth.mat', 'modelSalb');