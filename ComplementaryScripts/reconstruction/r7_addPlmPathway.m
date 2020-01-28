% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R7: Manual addition of paulomycin biosynthetic pathway

% Updated 2020-01-21
% Cheewin Kittikunapong

% Information regarding the pathway has been taken from Gonzalez et al.,
% 2016 (https://doi.org/10.1186/s12934-016-0452-4)

% The pathway has been further defined by inactivation of genes encoding 
% various enzymes involved in transfer of glycosyl and acyl moieties and 
% paulic acid and paulomycose biosynthesis.The study has enabled the 
% assignment of the genes to the individual steps along the biosynthetic 
% pathway based on identification of products and intermediates accumulated 
% in mutant strains. Using this information on the proposed pathway steps 
% and annotations that have been made for certain reactions and genes on 
% KEGG, the corresponding reactions of the paulomycin pathway have been 
% manually introduced into the metabolic network reconstruction.

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference

clear;
load('scrap/r6b_draftSalb_addKEGGRxns');
%% Addition of new metabolites in the pathway

% Metabolites information will be loaded from a text file. 

fid         = fopen('../../ComplementaryData/reconstruction/paulomycinMets.txt');
loadedData  = textscan(fid,'%q %q %q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_', 'HeaderLines', 1);
fclose(fid);

clear metsToAdd

% define the metsToAdd structure with metabolite names, compartment and
% chemical formulae
metsToAdd.metNames      = loadedData{1}; numel(metsToAdd.metNames)
metsToAdd.compartments  = loadedData{2};
metsToAdd.metFormulas   = loadedData{3};
metsToAdd.mets          = generateNewIds(modelSalb,'mets','s_',length(metsToAdd.metNames));
% there is a bug with generation of Ids for metabolites
metsToAdd.mets          = regexprep(metsToAdd.mets,'s','m');
model                   = addMets(modelSalb,metsToAdd); clear metsToAdd;

%% Addition of reactions in the paulomycin biosynthetic pathway

fid         = fopen('../../ComplementaryData/reconstruction/paulomycinRxns.txt');
loadedData  = textscan(fid,'%q %q %q %q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_', 'HeaderLines', 1);
fclose(fid);

%%

clear rxnsToAdd
rxnsToAdd.rxnNames      = regexprep(loadedData{1},'***','');
rxnsToAdd.equations     = regexprep(loadedData{2},'***','');
rxnsToAdd.grRules       = regexprep(loadedData{3},'***','');
rxnsToAdd.rxns          = generateNewIds(modelSalb,'rxns','r_',length(rxnsToAdd.rxnNames));
model                   = addRxns(model,rxnsToAdd,3,'',true,true); clear rxnsToAdd

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

modelSalb = model;

%% Save the model structure for later steps

save('scrap/r7_draftSalb_addPlmPathway', 'modelSalb');