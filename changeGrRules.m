function [model_new, Rules2Change] = changeGrRules(model, OldRule, NewRule)

% Cheewin Kittikunapong 2019-03-22

% define string (str) globally in Command Window

% modelSco = importModel('ScoGEM.xml', false);

% OldRule = ;

% NewRule = ;

% find all occurrences of the grRule string
Rules2Change = find(contains(model.grRules, OldRule));

% This code was kept from before converting to a function
% check rxns match with indices of grRule
% model.rxns{Rules2Change};

%% replace grRules with appropriately formatted rules

for i = 1:length(Rules2Change)
    

    model.rxns(Rules2Change(i))
    model.grRules{Rules2Change(i)} = NewRule;
    
end

% verify correct substitution
disp(model.rxns{Rules2Change});
disp(model.grRules{Rules2Change});

model_new = model;

end
