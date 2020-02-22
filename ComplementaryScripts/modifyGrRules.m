function model_new = modifyGrRules(model, OldRule, NewRule)
% Cheewin Kittikunapong 2019-03-22

% find all occurrences of the grRule string
Rules2Change = find(contains(model.grRules, OldRule));


%% replace grRules with appropriately formatted rules

for i = 1:length(Rules2Change)
    
    model.grRules{Rules2Change(i)} = NewRule;
    
end

model_new = model;

end
