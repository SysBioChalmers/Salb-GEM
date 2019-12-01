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

modelSalb = changeGrRules(modelSalb,rxns,newGrRules,true);

% Removing remaining SCO genes from model structure
geneIndex = find(contains(modelSalb.genes, 'SCO'));
modelSalb = removeGenes(modelSalb, modelSalb.genes(geneIndex));

% Removing any unused metabolite
modelSalb = removeMets(modelSalb,all(modelSalb.S == 0,2),false,true,true,true);

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
    if isempty(grRules{1})
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

%% Update MIRIAM annotations using our template model

% Certain reactions added after step R1 when we generated the draft model
% from homology with the template omitted any annotation information. We
% will update all known annotation using template as basis.

[match, matchIdx]   = ismember(modelSco.rxns,modelSalb.rxns);
modelSalb.rxnMiriams(matchIdx(match)) = modelSco.rxnMiriams(match);

[match, matchIdx]   = ismember(modelSco.mets,modelSalb.mets);
modelSalb.metMiriams(matchIdx(match)) = modelSco.metMiriams(match);

%% Save the model in .mat format

save('scrap/r5_draftSalb_curateAnnotations.mat','modelSalb');
