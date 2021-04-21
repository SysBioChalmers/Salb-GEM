% Script to analyze contribution of each KEGG rxn towards growth rate

% model = r5
load('reconstruction/scrap/r5_draftSalb_curateAnnotations.mat');
model = modelSalb;
% modelSalb = final
load('reconstruction/scrap/r7_draftSalb_addPlmPathway.mat');
% load list of KEGG rxns added
keggRxns = {'2AMACHYD','34DHOXPEGOX','3HOPCD','3OAR40_1','3HXKYNOXDA','3MOXTYROX','ACALD','41R2A1H12BOOX','42A12BOOX','4OT2','5HOXINOXDA','6TXAN5MPAML','AACTOOR','AGTi','ACAS_2ahbut','ASP4DC','ASPO2y','GLYTA','AMID2','AOBUTDs','APENTAMAH','APRTO2','KAS15','ASNTRAT','LALDO2x','SAGH','TREH','MMM2','BATHORHODOPSINI','BDG2HCGHD','BLACT','BTCOARx','CBMD','COA1819ZD9DS','CYSSADS','DDMCITALOR','DHBD','DMT','ECH_3hivcoa','SPT_syn','FALGTHLs','FDMO2','GLUKA','FFSD','G1PTT','LCGH','GALO','GLNSP2','GLNTRAT','P5CRx','MSAR','SPTc','GMHEPAT','PGMT_B','HEX10','HIBH','HKNDDH','HKNTDH','HMNO','HPROa','GUAD','ICHORT','AHSERL4','ASNN','FEROc','G3PT','LPPGS','LUMIRHODOPSINI','MCITD','MCITL2','MCITS','ASPO2','GLYCLTDy','NBAH','NFORGLUAH','NTRLASE','NTRLASE2','NTRLASE3','OPAH','AHAL','PDE1','PDH','PEAMNO','HSERTA','PTRCAT1','RAFH','RBK_Dr','RBK_L1','RBP4E','RBTDG','RDH1','RE0691C','RE2034C','3HPCOAHYD','LACZ','SELCYSS','SELGTHR','SELNPS','SERD_D','SPHPL','ATHRDHr','FADRx','THMP','SAGH_glycan','TRYPTAOX','TYROXDAc','VNTDM'};

%% Analyze flux through each rxn upon their addition
x = [];

for i = 1:length(keggRxns)
    model1 = addRxnsGenesMets(model, modelSalb, keggRxns(i), true);
    sol = solveLP(model1,1);
    % track growth rate upon addition of each rxn
    x = [x; -sol.f];
    idx = getIndexes(model1, keggRxns(i), 'rxns');
    % track flux through each rxn upon addition
    %x = [x; sol.x(idx)];
end

%% Analyze flux through each rxn in the final model
% Some rxns that carried flux in the previous section do not do so in the
% final model as there are alternate pathways available.
x = [];

for i = 1:length(keggRxns)
    idx = getIndexes(modelSalb, keggRxns(i), 'rxns');
    sol = solveLP(modelSalb,1);
    % track flux through rxn in the final model
    x = [x; sol.x(idx)];
end

%% Brute-force identification of problematic rxns
% From looking through the data from previous sections, we now have to
% identify the set of reactions responsible for the pronounced increase in
% growth and nonsensical production rate of malonyl-CoA.

% combination matrix of 1 and 0 to remove KEGG rxns combinatorially
rxnsToCheck = {'ASPO2y', 'ASPO2', 'FADRx', 'GLYTA', 'AGTi', 'SPT_syn', 'SPTc', 'GLYCLTDy', 'ASNN', 'GUAD', 'PGMT_B', 'P5CRx', 'HPROa'};% generate matrix
m = length(rxnsToCheck);
mat = logical(dec2bin(0:2^m-1,m)-'0'); % Thanks to Matt J: https://www.mathworks.com/matlabcentral/answers/114863-how-do-i-make-a-2-m-x-m-dimensional-matrix-containing-all-possible-combinations-of-1s-and-0s

%% Run the inefficient algorithm to find the reactions

mu = [];
q = [];
for i = 1:length(mat)
    model = removeReactions(modelSalb, rxnsToCheck(mat(i,:)));
    sol = solveLP(model);
    % track growth rates upon deletion
    mu = [mu; -sol.f];
    model = addExchangeRxn(model, 'malcoa_c', 0, 1000);
    model = changeRxnBounds(model, 'BIOMASS_SALB', -sol.f * 0.99, 'l');
    model = changeObjective(model, model.rxns(end));
    sol = solveLP(model);
    % track malonyl-CoA production upon deletion
    q = [q; -sol.f];
end

%% Analyze results

% Find the reaction combinations upon deletion resulted in a decrease in
% growth rate
muDiff = 0.1166 - mu;
idx = find(muDiff > 1e-4);

% find the combination with least rxns
n = [];
for i = 1:length(idx)
    n = [n; length(find(mat(idx(i),:)))];
end
% the greatest common denominator is P5CRx and HPROa
idxRxns = find(n == min(n));

%% Verify results

q = [];
for i = 1:length(idx)
    model = removeReactions(modelSalb, rxnsToCheck(mat(idx(i),:)));
    model = addExchangeRxn(model, 'malcoa_c', 0, 1000);
    model = changeRxnBounds(model, 'BIOMASS_SALB', mu(idx(i)) * 0.99, 'l');
    model = changeObjective(model, model.rxns(end));
    sol = solveLP(model);
    q = [q; -sol.f];
end
