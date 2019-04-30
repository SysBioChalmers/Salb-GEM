% Cheewin Kittikunapong
% Genome-scale metabolic model reconstruction of Streptomyces 'albus' J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11
% Script written 2019-04-12

%% Addition of exchange and pseudoreactions from scoGEM
%initCobraToolbox(0); % don't check for updates

% load the template model scoGEM
ScoGEM = importModel...
    ('../../ComplementaryData/reconstruction/templateModels/scoGEM.xml');

% load ScoHomology.mat draftSalb model
load('scrap/draft_ScoHomology.mat');

% find the exchange reactions and their indices in the template
[ex_rxns, ex_i] = getExchangeRxns(ScoGEM, 'both');

% identify the reaction names
ex_rxnNames = ScoGEM.rxnNames(ex_i);

% add exchange reactions from scoGEM to the Salb draft model
draftModel_ex = addRxnsGenesMets(draftSalb,ScoGEM, ex_rxns, 0, ...
    'Added to test for growth in draft model', 1);

% open the gates!
draftModel_ex.lb(1525:end) = -1000;

% find the pseudoreactions and biomass reaction
pseudoRxns = ScoGEM.rxns(1:15);

% add pseudoreactions and biomass reaction to Salb draft model
draftModel_biomass = addRxnsGenesMets(draftModel_ex, ScoGEM, pseudoRxns, ...
    0, 'Added to test for growth in draft model', 1);

%% Add non-GeneAssociated Rxns from Sco

% Find which rxns have no designated gene associations
noGene = find(cellfun(@isempty, ScoGEM.grRules));
Transport = find(getTransportRxns(ScoGEM));
Transport = intersect(Transport, noGene);
transportRxns = ScoGEM.rxns(Transport);
draftModel_new = addRxnsGenesMets(draftModel_biomass, ScoGEM, transportRxns, ...
    0, 'Added to test for growth in draft model', 1);


%% 26 apr 2019

LAminoAcids = {'L-Alanine' 'L-Arginine' 'L-Asparagine' 'L-Aspartate' 'L-Cysteine'...
    'L-Glutamine' 'L-Glutamate' 'Glycine' 'L-Histidine' 'L-Isoleucine' 'L-Leucine'...
    'L-Lysine' 'L-Methionine' 'L-Phenylalanine' 'L-Proline' 'L-Serine' ...
    'L-Threonine' 'L-Tryptophan' 'L-Tyrosine' 'L-Valine'};

% Mark indices of amino acids
amino = find(contains(draftModel_new.metNames,LAminoAcids));

% Find L-amino acids
draftModel_new.metNames(amino)

% Filter out duplicates, still containing AA derivatives containing string
% in MetName
unique(draftModel_new.metNames(amino))

%% Define objective function

% define the biomass reaction as objective function
testModel = setParam(draftModel_new, 'obj', 'BIOMASS_SCO', 1);
%testModel = changeObjective(draftModel_new, {'BIOMASS_SCO'});

% test FBA
%solutionFBA = optimizeCbModel(testModel);

[solution, hsSolOut] = solveLP(testModel, 0);

%% Save variables for later steps

save('scrap/draft_ScoHomology_gapFill', 'draftModel_new');