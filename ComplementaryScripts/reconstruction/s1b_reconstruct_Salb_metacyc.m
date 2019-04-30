% Cheewin Kittikunapong
% Genome-scale metabolic model reconstruction of Streptomyces 'albus' J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Reconstruction of draft model from MetaCyc

% Script written 2019-04-15

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% load FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = 'seq/GenBank/J1074_protein_fasta_sortLoci.fasta';

%% Build Salb draft model from MetaCyc database

% Summary (referred to as S. albus J1074):
% https://biocyc.org/organism-summary?object=SALB457425

% no exchange reactions (S 1799x1248)
METAmodel = getMetaCycModelForOrganism(organismID, fastaFile);
% METAmodel = getMetaCycModelForOrganism(organismID, fastaFile, 1); 

[METAModel, removedRxns]=removeBadRxns(METAmodel);
