function essentialModel=getEssentialModel(model, templateModel, growthRxn)
% getEssentialModel
%   Obtains the reduced subset of template model, having only the missing
%   reactions, which are missing in draftModel to satisfy the growth
%   objective. The function is useful, when one wants to get a list of
%   reactions from template model, which can be used for gap filling for
%   the reconstructed model 
%
%   model                   a model in reconstruction
%   templateModel           a template model, from which essential
%                           reactions should be obtained
%   growthRxn               growth reaction ID
%
%
%   essentialModel=getEssentialModel(model, templateModel)
%
%   Simonas Marcisauskas, 2015-09-18
% 	simmarc@chalmers.se
%
%   >>> removed media information (2016-04-18)
%   >>> added links to SM functions (2016-04-18)
%   >>> added growthRxn (2016-04-18)

sol=solveLP(templateModel,1);
templateModel=setParam(templateModel,'eq',growthRxn,abs(sol.f));

% Performing random sampling for the template model;
solutions=randomSampling(templateModel, 1000, false);

vector=zeros(numel(templateModel.rxns),1);
solutions=abs(solutions);
for i=1:numel(vector)
    if(sum(solutions(i,:))==0)
        vector(i)=1;
    end
end

reducedTemplateModel=removeReactions(templateModel, vector, true, true, true);
[essRxns]=ismember(reducedTemplateModel.rxns,model.rxns);
essentialModel=removeReactions(reducedTemplateModel,essRxns,true,true,true);