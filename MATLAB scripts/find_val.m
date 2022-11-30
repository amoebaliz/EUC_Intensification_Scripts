for j = 1:length(C)
    k(j) = length(C{j});
end

I = find(k<length(C{1}));
[length(C{1}) I(1)]
