% Genome-scale metabolic model reconstruction of Streptomyces albus J1074
% (As of October 2018 referred to as Streptomyces albidoflavus)

% Project started 2019-03-11

% Script R7: Manual addition of paulomycin biosynthetic pathway

% Updated 2020-01-21
% Cheewin Kittikunapong

%% Suggestion:
% You may want to clear your workspace to avoid an excessive number of
% variables for ease of reference

clear;

load('scrap/r7_draftSalb_addPlmPathway');
modelSco = importModel('../../ComplementaryData/reconstruction/templateModels/Sco-GEM.xml');

% Convert to COBRA format for ease of modification of fields
modelCb = ravenCobraWrapper(modelSalb);

% NOTE:
% edits have been made to the wrapper function in order to be able to parse
% fields for pfam, uniprot, ec, KO and ncbi IDs

%% Annotate each gene with KEGG gene IDs ('[organism id]:')

for i = 1:numel(modelCb.genes)
    modelCb.geneiskegg__46__genesID(i,1) = strcat('salb:',modelCb.genes(i));
end

%% Load EC number information from KEGG FTP and annotate corresponding genes

fid         = fopen(['../../ComplementaryData/reconstruction/KEGGFTP/salb_enzyme.list']);
loadedData  = textscan(fid,'%s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
genes    = loadedData{1}; genes = regexprep(genes,'salb:','');
EC           = loadedData{2}; EC = regexprep(EC,'ec:','');

geneIdx = find(ismember(modelCb.genes, genes));

modelCb.geneisecID = cell(numel(modelCb.genes),1);

for i = 1:numel(geneIdx)
    j = find(ismember(genes, modelCb.genes(geneIdx(i))));
        modelCb.geneisecID(geneIdx(i)) = join(EC(j),'; ');
end

%% Load KO information from KEGG FTP and annotate corresponding genes

fid         = fopen(['../../ComplementaryData/reconstruction/KEGGFTP/salb_ko.list']);
loadedData  = textscan(fid,'%s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
genes        = loadedData{1}; genes = regexprep(genes,'salb:','');
KO           = loadedData{2}; KO = regexprep(KO,'ko:','');

geneIdx = find(ismember(modelCb.genes, genes));

modelCb.geneiskoID = cell(numel(modelCb.genes),1);

for i = 1:numel(geneIdx)
    j = find(ismember(genes, modelCb.genes(geneIdx(i))));
        modelCb.geneiskoID{geneIdx(i)} = strjoin(KO(j),'; ');
end

%% Load pfam information from KEGG FTP and annotate corresponding genes

fid         = fopen(['../../ComplementaryData/reconstruction/KEGGFTP/salb_pfam.list']);
loadedData  = textscan(fid,'%s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
genes        = loadedData{1}; genes = regexprep(genes,'salb:','');
pfam           = loadedData{2}; pfam = regexprep(pfam,'pf:','');

geneIdx = find(ismember(modelCb.genes, genes));

modelCb.geneispfamID = cell(numel(modelCb.genes),1);

for i = 1:numel(geneIdx)
    j = find(ismember(genes, modelCb.genes(geneIdx(i))));
        modelCb.geneispfamID{geneIdx(i)} = join(pfam(j),'; ');
end

%% Load UniProt information from KEGG FTP and annotate corresponding genes

fid         = fopen(['../../ComplementaryData/reconstruction/KEGGFTP/salb_uniprot.list']);
loadedData  = textscan(fid,'%s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
genes        = loadedData{1}; genes = regexprep(genes,'salb:','');
uprot           = loadedData{2}; uprot = regexprep(uprot,'up:','');

geneIdx = find(ismember(modelCb.genes, genes));

modelCb.proteinisuniprotID = cell(numel(modelCb.genes),1);

for i = 1:numel(geneIdx)
    j = find(ismember(genes, modelCb.genes(geneIdx(i))));
        modelCb.proteinisuniprotID{geneIdx(i)} = strjoin(uprot(j),'; ');
end

%% Load NCBI ID information from KEGG FTP and annotate corresponding genes

fid         = fopen(['../../ComplementaryData/reconstruction/KEGGFTP/salb_ncbi-proteinid.list']);
loadedData  = textscan(fid,'%s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
genes        = loadedData{1}; genes = regexprep(genes,'salb:','');
ncbi           = loadedData{2}; ncbi = regexprep(ncbi,'ncbi-proteinid:','');

geneIdx = find(ismember(modelCb.genes, genes));

modelCb.geneisncbiID = cell(numel(modelCb.genes),1);

for i = 1:numel(geneIdx)
    j = find(ismember(genes, modelCb.genes(geneIdx(i))));
        modelCb.geneisncbiID{geneIdx(i)} = strjoin(ncbi(j),'; ');
end

%% Load GO number information from KEGG FTP and annotate corresponding genes

fid         = fopen(['../../ComplementaryData/reconstruction/XNR_GO.tsv']);
loadedData  = textscan(fid,'%s %s %s','delimiter','\t', 'HeaderLines',1); fclose(fid);
genes        = loadedData{1};
GO           = loadedData{2};

geneIdx = find(ismember(modelCb.genes, genes));

modelCb.geneisgoID = cell(numel(modelCb.genes),1);

for i = 1:numel(geneIdx)
    j = find(ismember(genes, modelCb.genes(geneIdx(i))));
        modelCb.geneisgoID{geneIdx(i)} = strjoin(GO(j),'; ');
end

%% Addition of annotation fields to the RAVEN model structure

% convert the modified/annotated COBRA model back to RAVEN format
model = ravenCobraWrapper(modelCb);

% We should then add the fields of this new RAVEN model to the original
% model. This is done because the use of ravenCobraWrapper will cause some
% information to be lost across the two conversions
modelSalb.geneMiriams = model.geneMiriams;

%% Save model structure

save('scrap/r8_draftSalb_annotateGenes', 'modelSalb');