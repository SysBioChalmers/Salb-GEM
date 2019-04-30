% Cheewin Kittikunapong

% Script written 2019-04-25

%initCobraToolbox(0);

%% Loading Sco and Sclav as templates

%modelSco = importModel('template_models/ScoGEM.xml', false);

% load standardized grRules for modelSco
modelSco = importModel...
    ('../../ComplementaryData/reconstruction/templateModels/scoGEM_stdGrRules.xml',false);

% modelSclav = importModel('template_models/SclavGEM.xml', false, true);
% NOTE 2019-04-25
% remove Lines 598-600 to import modelSclav via RAVEN
% SBML2 format, missing metabolite IDs in S matrix

% ad hoc solution: read and rewrite in CobraToolbox then import in RAVEN
%model = readCbModel('../../ComplementaryData/reconstruction/templateModels/SclavGEM.xml');
%writeCbModel(model, 'fileName','SclavGEMCb.xml','format','sbml');

modelSclavCb = importModel('../../ComplementaryData/reconstruction/templateModels/SclavGEMCb.xml', false);

% verify field elements not missing
% exportToExcelFormat(modelSclavCb, 'SclavCB');

% load rewritten template model for modelSclav
modelSclav = modelSclavCb; clear('modelSclavCb');
modelSclav.id = 'SclavGEM';

%% Designate parameters and information

% Source: https://www.genome.jp/kegg-bin/show_organism?org=salb

organismID = 'salb';

% load FASTA file with headers sorted by gene loci (XNR_xxxx)
fastaFile = '../../ComplementaryData/genome/J1074_protein_fasta_sortLoci.txt';

%% Build Salb draft model from Sco and Sclav homology

% Error using getModelFromHomology (line 129)
%Less than 5% of the genes in the template model with model.id "Streptomyces_clavuligerus_model_Sclav_iLT1021_xml"
%can be found in the blastStructure. Ensure that the protein FASTA
%used in getBlast and the template model used in getModelFromHomology
%use the same style of gene identifiers

% cell array for reference organisms structures and IDs

refModels = {modelSco, modelSclav};
refModelIDs = {modelSco.id, modelSclav.id};

% cell array for reference FASTA files may be appended with more templates
refFastaFiles = {'../../ComplementaryData/genome/Sco_all_protein.faa',...
    '../../ComplementaryData/genome/Sclav_all_protein.faa'};

% blastStructure has been exported manually into an Excel file for closer
% inspection
blastStructure = getBlast(organismID, fastaFile, refModelIDs, refFastaFiles);

%% generate a draft model based on homology with default parameters 
% (E value, alignment length and identity)

draftSalb = getModelFromHomology(refModels, blastStructure, organismID);

% NOTE: 2019-04-30
% different metabolite IDs (metaid, id, name)
% different reactions IDs
