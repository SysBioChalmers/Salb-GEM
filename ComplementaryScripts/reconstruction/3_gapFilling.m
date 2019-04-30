% Cheewin Kittikunapong
% Genome-scale metabolic model reconstruction of Streptomyces 'albus' J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11
% Script written 2019-04-16

%% Find all blocked metabolites

load('scrap/draft_ScoHomology_gapFill.mat');
testModel = setParam(draftModel_new, 'obj', 'BIOMASS_SCO', 1);

% identify deadend metabolities

outputMets = detectDeadEnds(testModel);
DeadEnds = testModel.mets(outputMets)

% identify corresponding reactions / find blocked reactions

[rxnList, rxnFormulaList] = findRxnsFromMets(testModel, DeadEnds)
%[allGaps, rootGaps, downstreamGaps] = gapFind(testModel);

%% Using COBRA fastGapFill

%[consistModel, consistMatricesSUX, BlockedRxns] = prepareFastGapFill...
%    (model, listCompartments, epsilon, filename, dictionary_file, blackList)

[consistModel. consistMatricesSUX, BlockedRxns] = prepareFastGapFill...
    (testModel);

% NOTE 19-04-16
% INFEASIBLE solution
% similar case to RAVEN fillGaps

%% Load template model

load('scrap/templateSco.mat');
template = modelSco;

%% Using RAVEN fillGaps function
%[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps...
%    (model,models,allowNetProduction,useModelConstraints,supressWarnings,rxnScores,params)

[ConnectedRxns, cannotConnect, addedRxns, newModel, exitFlag] = fillGaps...
    (testModel, template);

% NOTE 19-04-16
% error in fillGaps (line 188)
% originalFlux=haveFlux...

% error in haveFlux (line 55)
% J(abs(sol.x(minIndexes))>cutOff) = true
% sol is empty due to infeasible problem.. why?

% NOTE 19-04-17
% Opened up all sco ex_Rxns
% returned:
%   - 483 Connected Rxns
%   - 897 cannotConnect
%   - 0 added Rxns

% NOTE 19-04-26
% Added transport reactions
% returned:
%   - 433 Connected Rxns
%   - 923 cannotConnect
%   - 0 added Rxns

% NOTE 19-04-29
% Generated merged model with essentialModel (SM function)
% returned:
%   - 243 Connected Rxns
%   - 880 cannotConnect

%% Get Essential Model for draftSalb

% essentialModel = getEssentialModel(model, template, biomassRxn)
addpath('../../ComplementaryScripts');
essentialModel = getEssentialModel(testModel, template,'BIOMASS_SCO');

newModel = mergeModels({draftModel_new, essentialModel});

[solution, hsSolOut] = solveLP(newModel);

% 425 deadEnd Metabolites
% 358 blockedRxns

% no growth