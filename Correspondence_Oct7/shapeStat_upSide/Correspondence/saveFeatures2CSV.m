function saveFeatures2CSV(spokes, links)
feature_mat = [spokes; links];
[d, n] = size(feature_mat);
y_mat = ones(1, n);
feature_mat = [feature_mat; y_mat];

csvwrite('spikeslab.csv', feature_mat);
end