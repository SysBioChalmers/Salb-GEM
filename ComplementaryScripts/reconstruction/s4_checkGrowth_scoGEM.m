% This scripts examines the current template model's predictions for
% certain carbon sources, which have been observed for S. albus

% Cheewin Kittikunapong
% 2019-09-06

%% Checking growth of Sco-GEM on known S. albus carbon sources
clear;

% first load the Sco-GEM model from our reconstruction steps
load('scrap/r1_scoGEM_newGrRules.mat');
modelSco = rmfield(modelSco, 'unconstrained');
modelSco = setParam(modelSco,'lb','BIOMASS_SCO',0);
modelSco = setParam(modelSco,'ub','BIOMASS_SCO',1000);
modelSco = setParam(modelSco, 'obj', 'BIOMASS_SCO', 1);

% We will be looking for the reactions that have not been added to Salb-GEM
% thus far. To find this information, we must load the Salb-GEM draft we
% have generated after the gap-filling step, R3.
load('scrap/r3_draftSalb_gapFill.mat');

% identify the reactions that do not exist in Salb-GEM so far
% We will see which subset of reactions in this list does carry flux for
% various carbon source conditions
missingRxns = setdiff(modelSco.rxns, modelSalb.rxns);

% First, verify growth on glucose
modelSco = setParam(modelSco, 'lb', 'EX_glc__D_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -sol.f
printFluxes(modelSco, sol.x,true,10^-8)

%% switch carbon sources
modelSco = setParam(modelSco, 'lb', 'EX_glc__D_e', 0);
reset = modelSco;

%% Growth on fructose only

% In case we need to back-track quickly
modelSco = reset;

% growth on fructose
modelSco = setParam(modelSco, 'lb', 'EX_fru_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -1 .* sol.f
printFluxes(modelSco, sol.x,true,10^-8)

% The fluxes through biomass between FBA with glucose and fructose exchanges 
% are virtually identical. This should be expected given how close the
% pathways are with each other.

% Now to identify the reactions carrying flux that we may consider adding
% to our draft model
rxnIndexes = find(sol.x); missingRxns_fru = modelSco.rxns(rxnIndexes);
missingRxns_fruPTS2 = intersect(missingRxns, missingRxns_fru);
save('scrap/s4_missingRxns_fruPTS2.mat', 'missingRxns_fruPTS2');

%% Growth on fructose only (via a different uptake route)
% There currently lacks unequivocal evidence on how exactly fructose is
% taken up into the cell in S. coelicolor, although we do know that it is
% via a PEP:phosphotransferase system common in bacteria.

% Two separate models from which Sco-GEM is derived have specified
% different routes of uptake for fructose, either converting extracellular
% fructose to intracellular fructose 6-phosphate or fructose 1-phosphate

% While both reactions do exist in the model, flux only preferentially
% travels through FRUpts2, which produces F6P

% By removing the fructose PEP:phosphotransferase FRUpts2, the flux must
% go through fructose 1-phosphate instead of fructose 6-phosphate

modelSco = removeReactions(modelSco, 'FRUpts2');

modelSco = setParam(modelSco, 'lb', 'EX_fru_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -1 .* sol.f
printFluxes(modelSco, sol.x,true,10^-8)

% After identifying the reactions, we can see that the subset through
% FRUpts1 requires more steps, which could be why the model only prefers
% flux be carried through FRUpts2

rxnIndexes = find(sol.x); missingRxns_fru = modelSco.rxns(rxnIndexes);
missingRxns_fruPTS1 = intersect(missingRxns, missingRxns_fru);
save('scrap/s4_missingRxns_fruPTS1.mat', 'missingRxns_fruPTS1');

%% Growth on mannitol only

modelSco = reset;
% growth on mannitol
modelSco = setParam(modelSco, 'lb', 'EX_mnl_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -1 .* sol.f
printFluxes(modelSco, sol.x,true,10^-8)

% different fluxes here, makes sense - different pathways altogether

rxnIndexes = find(sol.x); missingRxns_mnl = modelSco.rxns(rxnIndexes);
missingRxns_mnl = intersect(missingRxns, missingRxns_mnl);
save('scrap/s4_missingRxns_mnl.mat', 'missingRxns_mnl');

%% Growth on xylose only

modelSco = reset;
% growth on xylose
modelSco = setParam(modelSco, 'lb', 'EX_xyl__D_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -1 .* sol.f
printFluxes(modelSco, sol.x,true,10^-8)

rxnIndexes = find(sol.x); missingRxns_xyl = modelSco.rxns(rxnIndexes);
missingRxns_xyl = intersect(missingRxns, missingRxns_xyl);
save('scrap/s4_missingRxns_xyl.mat', 'missingRxns_xyl');

%% Growth on inositol only

modelSco = reset;
% growth on inositol
modelSco = setParam(modelSco, 'lb', 'EX_inost_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -1 .* sol.f
printFluxes(modelSco, sol.x,true,10^-8)

%% Growth on rhamnose only

modelSco = reset;
% growth on rhamnose
modelSco = setParam(modelSco, 'lb', 'EX_rmn_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -1 .* sol.f
printFluxes(modelSco, sol.x,true,10^-8)

%% Growth on sucrose only

modelSco = reset;
% growth on sucrose
modelSco = setParam(modelSco, 'lb', 'EX_sucr_e', -0.8);
[sol, ~] = solveLP(modelSco, 0); -1 .* sol.f
%printFluxes(modelSco, sol.x,true,10^-8)

% infeasible,does not grow on sucrose
% sucrose transporter (SUCRabc) is bounded to 0 in scoGEM