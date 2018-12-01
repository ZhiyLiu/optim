close all;
mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{'CPNS scores'}};
legendstr = {{''}};
[U S V] = svd(ZComp);
drawScatPlot(U, mcolor, legendstr, markerstr, titlestr);