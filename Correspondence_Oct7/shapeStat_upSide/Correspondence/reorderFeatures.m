% this function group the features as [x1, x2, ....xn, y1, y2, ...,yn, z1,
% z2...,zn]'
function f_reordered = reorderFeatures(f_orig, dim)
[num_features, num_samples] = size(f_orig);
num_groups = num_features / dim;
f_reordered = [];
for i = 1: dim
    for j = 1: num_groups
        f_vector = f_orig((j-1) * dim + i, :);
        f_reordered = [f_reordered; f_vector];
    end
end
end