% Function [M] = mrepExp(Mlist)
%
%   mrepExp(Mlist) computes the expmap of mreps in Mlist.
%
% INPUT : Mlist - cell array of combined(tube, quad) atoms
% OUTPUT : M - cell array of combined(tube, quad) atoms

function M = mrepExp(Mlist)

[row, col] = size(Mlist);

for r = 1:row
    for c = 1:col
    M{r, c} = Exp(Mlist{r, c}); 
    end
end

return;