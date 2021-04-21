% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R5: Curation of annotations for reactions added for function

% Updated 2019-10-21
% Cheewin Kittikunapong

% KEGG IDs
% sco = S. coelicolor
% salb = S. albus J1074

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference
clear;
load('scrap/r4_draftSalb_checkGrowth');
load('scrap/r1_scoGEM_newGrRules');

modelSco_new = importModel('../../ComplementaryData/reconstruction/templateModels/Sco-GEM_1.3.0.xml');

%% Reactions with grRules from template

% With RAVEN fillGaps, we added reactions with grRules in SCO terms.
% Additionally, we had also imported reactions to support fructose and 
% mannitol metabolism. Before proceeding further, we must match the 
% corresppmdomg SCO genes to their XNR counterparts.

rxnIdx      = strfind(modelSalb.grRules,'SCO');
rxnIdx      = ~cellfun('isempty',rxnIdx); % Get reaction indexes

fprintf('These reactions from gap-filling have to be checked for GPR:\n');
disp(modelSalb.rxns(rxnIdx));

%% Correcting the grRules

% Importing the list of updated grRules
fid         = fopen(['../../ComplementaryData/reconstruction/updatedGrRules.csv']);
loadedData  = textscan(fid,'%s %s %s %s %s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
rxns        = loadedData{1};
newGrRules  = loadedData{4};
for i = 1:length(rxns)
    if find(ismember(modelSalb.rxns, rxns(i)))
        modelSalb = changeGrRules(modelSalb,rxns(i),newGrRules(i),true);
    end
end
%%
% Removing remaining SCO genes from model structure
geneIndex = find(contains(modelSalb.genes, 'SCO'));
modelSalb = removeGenes(modelSalb, modelSalb.genes(geneIndex));

% Removing any unused metabolite
% modelSalb = removeMets(modelSalb,all(modelSalb.S == 0,2),false,true,true,true);

%% Remove unconnected non-gene associated reactions
subGraphs = getAllSubGraphs(modelSalb);

% Find which reactions have no gene associated
rxnToRemove     = [];
for i=2:size(subGraphs,2)
    metIdx      = subGraphs(:,i);
    rxnIdx      = modelSalb.S(metIdx,:);
    [~,col,~]   = find(rxnIdx);
    col         = unique(col);
    grRules     = modelSalb.grRules(col);
    if isempty(grRules)
         rxnToRemove = [rxnToRemove; col];
    end
end

modelSalb = removeReactions(modelSalb,rxnToRemove,true,true,true);

%% Correct faulty grRules where the same complex is representated multiple times

for n = 1:length(modelSalb.grRules)
    if any(modelSalb.grRules{n})
        noAnd = strfind(modelSalb.grRules(n),'and');
        noAnd = any(vertcat(noAnd{:})); % Give 0 if no 'and' is present.
        if noAnd == 0
            geneList = transpose(cell(unique(regexp(modelSalb.grRules{n},'[)(]*|( and )*|( or )*','split'))));
            geneList = regexprep(geneList,'[(*)*]','');
            if length(geneList) == 1
                newgrRule = geneList;
            else
                newgrRule = geneList{1};
                for k = 2:length(geneList)
                    newgrRule = [newgrRule ' or ' geneList{k}];
                end
            end
            modelSalb.grRules(n) = cellstr(newgrRule);
        end
    end
end

%% Correcting the grRules for complexes
% In later analyses, it was discovered that the grRules for complexes
% involving parantheses and ANDs were not accounted for by the previous
% step, thus redundant GPRs were assigned. More importantly, GPRs may also
% be erroneous if isozymes existed for a particular homologue pair
% resulting in GPRs that were written incorrectly upon obtainment of the
% initial draft model from protein homology with Sco-GEM.

% Importing the list of fixed grRules
fid         = fopen(['../../ComplementaryData/reconstruction/fixedGrRules_complexes.csv']);
loadedData  = textscan(fid,'%s %s %s %s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
rxns        = loadedData{1};
newGrRules  = loadedData{5};
for i = 1:length(rxns)
    if find(ismember(modelSalb.rxns, rxns(i)))
        modelSalb = changeGrRules(modelSalb,rxns(i),newGrRules(i),true);
    end
end

%% Update MIRIAM annotations using our template model

% Certain reactions added after step R1 when we generated the draft model
% from homology with the template omitted any annotation information. We
% will update all known annotation using template as basis.

[match, matchIdx]   = ismember(modelSco.rxns,modelSalb.rxns);
modelSalb.rxnMiriams(matchIdx(match)) = modelSco.rxnMiriams(match);

[match, matchIdx]   = ismember(modelSco.mets,modelSalb.mets);
modelSalb.metMiriams(matchIdx(match)) = modelSco.metMiriams(match);

%% Update metabolite species with chemical formulae annotation
mets = intersect(modelSalb.mets, modelSco_new.mets);
idx1 = getIndexes(modelSalb, mets, 'mets');
idx2 = getIndexes(modelSco_new, mets, 'mets');
modelSalb.metFormulas(idx1) = modelSco_new.metFormulas(idx2);

%% Save the model in .mat format

save('scrap/r5_draftSalb_curateAnnotations.mat','modelSalb');
