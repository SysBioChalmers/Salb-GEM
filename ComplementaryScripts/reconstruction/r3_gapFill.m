% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R3: Gap-filling of draft model using template model reactions

% Updated 2019-10-21
% Cheewin Kittikunapong

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference
clear

% For this step, you will need to change the parameters for mergeModels
% in RAVEN fillGaps to sort by metabolite IDs rather than by metabolite
% names. This is due to how Sco-GEM is structured by having non-unique
% metabolite names. If not accounted for, mergeModels will incorrectly omit
% some metabolites while parsing.

% Go to path (~/RAVEN/core/fillGaps.m) and correct the following line 144:
% From:
% allModels=mergeModels([{model};models],'metNames',true);
% To:
% allModels=mergeModels([{model};models],'mets',true);


%% load models from previous steps

% load current draft of modelSalb from step R2
load('scrap/r2_draftSalb_addScoRxns.mat');

% load template with grRules from step R1
load('scrap/r1_scoGEM_newGrRules.mat');
modelSco = rmfield(modelSco, 'unconstrained');

% We will remove ATP maintenance and exchange reactions from modelSco to not
% interfere with gapFill. Otherwise, these reactions may be suggested
% erroneously.

ex_rxns = getExchangeRxns(modelSco);
modelSco = removeReactions(modelSco, 'ATPM');
modelSco = removeReactions(modelSco, ex_rxns);

%% Removing Sco-GEM tailored pseudoreactions

% As some pseudoreactions were implemented specifically in the context of
% Sco-GEM, we will remove these reactions prior to the gap-filling step

modelSalb = removeReactions(modelSalb, 'BIOMASS_SCO_tRNA');
modelSalb = removeReactions(modelSalb, 'PROTEIN_PSEUDO_tRNA');

% For ease of later reference, we will also changed how the biomass
% reaction is structured by adding an actual pseudometabolite for biomass
modelSalb = changeRxns(modelSalb, 'BIOMASS_SCO', '75.79 ATP[c] + 75.79 H2O[c] + carbohydrate pseudometabolite[c] + cell wall pseudometabolite[c] + DNA pseudometabolite[c] + lipid pseudometabolite[c] + misc pseudometabolite[c] + protein pseudometabolite[c] + RNA pseudometabolite[c] => 75.79 ADP[c] + 75.79 H+[c] + 75.79 Phosphate[c] + biomass[c]', 3, '', true);
modelSalb.mets{end} = 'biomass_c';
i = getIndexes(modelSalb, 'BIOMASS_SCO', 'rxns');
modelSalb.rxns{i} = 'BIOMASS_SALB'; clear('i');

% add an exchange reaction for biomass labelled as growth
% This is done so that in later steps of analysis, we can easily view the
% fluxes of the biomass function and exchange reactions using printFluxes
modelSalb = addExchangeRxns(modelSalb, 'out','biomass_c'); 
modelSalb.rxns{end} = 'growth'; modelSalb.rxnNames{end} = 'biomass equation';

%% Updating biomass pseudoreactions to S. albus profile

% We will also be updating some pseudoreactions for DNA, RNA and amino
% acids with available information using the S. albus genome

% Load biomass information
fid             = fopen([ '../../ComplementaryData/biomass/biomassCuration.csv']);
loadedData      = textscan(fid, '%q %q %q %f','delimiter', ',', 'HeaderLines', 1);
fclose(fid);

biomass.name         = loadedData{1};    
biomass.mets         = loadedData{2};
biomass.pseudorxn    = loadedData{3};    
biomass.coeff        = loadedData{4};

% DNA pseudoreaction
% Obtain the relevant rows for information
indexes = find(contains(biomass.pseudorxn, 'DNA'));
% Update stoichiometries
equations.mets          = biomass.mets(indexes);
equations.stoichCoeffs  = biomass.coeff(indexes);
% Update changes in the model structure
modelSalb = changeRxns(modelSalb, 'DNA_PSEUDO', equations, 1);

% RNA pseudoreaction
indexes = find(contains(biomass.pseudorxn, 'RNA'));
equations.mets          = biomass.mets(indexes);
equations.stoichCoeffs  = biomass.coeff(indexes);
modelSalb = changeRxns(modelSalb, 'RNA_PSEUDO', equations, 1);

% Protein pseudoreaction
indexes = find(contains(biomass.pseudorxn, 'protein'));
equations.mets          = biomass.mets(indexes);
equations.stoichCoeffs  = biomass.coeff(indexes);
modelSalb = changeRxns(modelSalb, 'PROTEIN_PSEUDO', equations, 1);

% You can verify proper modification using below function
% constructEquations(modelSalb, '[macromolecule]_PSEUDO')

%% Use RAVEN fillGaps

% use biomass production as objective function for gap-filling
modelSalb = setParam(modelSalb,'obj','BIOMASS_SALB',0);
modelSalb = setParam(modelSalb,'obj','growth',1);

% set biomass production to arbitrary low flux, to force gap-filling to
% produce biomass in the draft model
modelSalb = setParam(modelSalb,'lb','growth',0.01);
modelSalb = setParam(modelSalb,'ub','growth',1000);

% Note from fillGaps documentation:
% define a biomass equation and set the lower bound to >0 with 
% useModelConstraints=true would then give the smallest set of reactions
% that have to be included in order for the model to produce biomass

[ConnectedRxns, cannotConnect, addedRxns, modelSalb, exitFlag] = fillGaps(modelSalb, modelSco, false, true);
modelSalb.id = 'Salb-GEM'
rxnIndexes = getIndexes(modelSalb, addedRxns, 'rxns');

modelSalb.rxnNotes(rxnIndexes) = {'Included from gap-filling with RAVEN'};

%% Verify growth by FBA and examine fluxes

modelSalb = setParam(modelSalb,'lb','growth',0);

[solution, hsSolOut] = solveLP(modelSalb, 0); -solution.f

% see what the exchange fluxes are (glucose uptake, etc.)

% NOTE:
% If you did not change the parameters of mergeModel in RAVEN fillGaps as
% specified in the SUGGESTION section at the top, your output from
% printFluxes will return fluxes for reactions that are not supposed to be
% exchange reactions. This is because RAVEN incorrectly parsed the
% reactions during the gap-filling stage and omitted certain metabolites
% from being incorporated into the model. Specifically for the reactions
% you see in printFluxes, their equations are missing one complete side
% leading it to be confused for an exchange reaction.

printFluxes(modelSalb, solution.x, true, 10^-8, '', '%rxnID\t %flux\t lb=%lower \t ub=%upper\n');

%% Check and validate added reactions in gap-filling

% some added reactions are (sometimes erroneous) duplicates of rxns
% existing already in the model

% (This is an artifact introduced by how the current mergeModels function
% operates; can address using a modified script found on fix/mergeModels)

%dup = find(contains(modelSalb.rxns, '_scoGEM'));
%modelSalb = removeReactions(modelSalb, modelSalb.rxns(dup), true, true, true);

% see flux distribution along added reactions from gapFill
% remove duplication reactions from added reactions list
%dup = find(contains(addedRxns, '_scoGEM'));
%addedRxns(dup) = [];

index = getIndexes(modelSalb, addedRxns, 'rxns');
q_addedRxns = solution.x(index);
names_addedRxns = modelSalb.rxnNames(index);

%% save model

%save('scrap/r3_draftSalb_gapFill.mat', 'modelSalb');
