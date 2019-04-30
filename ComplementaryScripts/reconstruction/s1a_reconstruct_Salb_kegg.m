% Cheewin Kittikunapong
% Genome-scale metabolic model reconstruction of Streptomyces 'albus' J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Reconstruction of draft model from KEGG

% KEGG IDs
% sco = S. coelicolor
% salb = S. albus J1074

% Script written 2019-04-15

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% load FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = 'seq/GenBank/J1074_protein_fasta_sortLoci.fasta';

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