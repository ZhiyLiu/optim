close all;
%%
figure; hold on;
draw_data = ZComp;
legendstr = {{''}};

[d n] = size(draw_data);
mcolor = [];
markerstr = [];
group = 3;
num_per_group = n/group;
iter_group = 1;
for i = 1:n
    % 2 groups
    if i < n/2 + 1
        markerstr = strvcat(markerstr, 'o');
        mcolor = [mcolor; 1 0 0];
    else
        markerstr = strvcat(markerstr, 'x');
        mcolor = [mcolor; 0 0 0];
    end
    
    % 3 groups
%     if i < n/3 + 1
%         markerstr = strvcat(markerstr, 'o');
%         mcolor = [mcolor; 1 0 0];
%     elseif i > n/3 && i < 2*n/3 + 1
%         markerstr = strvcat(markerstr, '+');
%         mcolor = [mcolor; 0 0 1];
%     else
%         markerstr = strvcat(markerstr, 'x');
%         mcolor = [mcolor; 0 0 0];
%     end
end

titlestr = {{'Scatter plot of ZComp'}}
%[U S V]=svd(draw_data);
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);
