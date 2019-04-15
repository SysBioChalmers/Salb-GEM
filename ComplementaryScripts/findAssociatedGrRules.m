function [where, total] = findAssociatedGrRules(model, grRule)

% Cheewin Kittikunapong 2019-03-18

% load grRule by default or insert set of genes of interest as cell array
% function to see all occurrences of a gene in reaction grRule

%(SCO5281 and (SCO2181 or SCO7123 or SCO1268) and (SCO0884 or SCO2180 or SCO4919))
toRemove = ["(", ")"];
toSplit = [" and ", " or "];
string = erase(grRule, toRemove);
string = split(string, toSplit);

length(string);

for i = 1:length(string)
    a = contains(model.grRules,string(i));
    matrix(:,i) = a;
    where{i,:} = find(a);
end

total = sum(matrix);

for i = 1:length(string)
    fprintf('Gene %s is involved in %d grRule(s) at indices: \n', string{i}, total(i));

    disp(where{i})
end

