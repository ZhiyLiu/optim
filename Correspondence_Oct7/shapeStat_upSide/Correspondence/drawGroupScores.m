close all;
%%
figure; hold on;
draw_data = outstruct.MatrixJoint{1,1};
legendstr = {{'orginal 10 cases' 'Shift group A'}};

[d n] = size(draw_data);
mcolor = [];
markerstr = [];
for i = 1:n
    if i < n/2 + 1
        markerstr = strvcat(markerstr, 'o');
        mcolor = [mcolor; 1 0 0];
    else
        markerstr = strvcat(markerstr, '+');
        mcolor = [mcolor; 0 0 1];
    end
end

titlestr = {{'Joint matrix for data block 1'}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);

%%
figure; hold on;
draw_data = outstruct.MatrixJoint{1,2};
legendstr = {{'orginal 10 cases' 'Shift group A'}};

[d n] = size(draw_data);
mcolor = [];
markerstr = [];
for i = 1:n
    if i < n/2 + 1
        markerstr = strvcat(markerstr, 'o');
        mcolor = [mcolor; 1 0 0];
    else
        markerstr = strvcat(markerstr, '+');
        mcolor = [mcolor; 0 0 1];
    end
end

titlestr = {{'Joint matrix for data block 2'}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);

%%
figure; hold on;
draw_data = outstruct.MatrixIndiv{1,1};
legendstr = {{'orginal 10 cases' 'Shift group A'}};

[d n] = size(draw_data);
mcolor = [];
markerstr = [];
for i = 1:n
    if i < n/2 + 1
        markerstr = strvcat(markerstr, 'o');
        mcolor = [mcolor; 1 0 0];
    else
        markerstr = strvcat(markerstr, '+');
        mcolor = [mcolor; 0 0 1];
    end
end

titlestr = {{'Individual matrix for data block 1'}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);

%%
figure; hold on;
draw_data = outstruct.MatrixIndiv{1,2};
legendstr = {{'orginal 10 cases' 'Shift group A'}};

[d n] = size(draw_data);
mcolor = [];
markerstr = [];
for i = 1:n
    if i < n/2 + 1
        markerstr = strvcat(markerstr, 'o');
        mcolor = [mcolor; 1 0 0];
    else
        markerstr = strvcat(markerstr, '+');
        mcolor = [mcolor; 0 0 1];
    end
end

titlestr = {{'Individual matrix for data block 2'}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);
