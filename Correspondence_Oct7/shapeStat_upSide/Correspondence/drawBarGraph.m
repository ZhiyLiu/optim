% usually loadings are much more important in homogeneous population, since
% it makes more sense to measure for each single feature how it varies among
% objects. In loadings, each row is one primitive in graph
function drawBarGraph(data, group_size, feature_name, barcolor, titlestr)

if(length(data) > 0)
    
    for i = 1:length(feature_name)
        xdata = (i-1)*group_size + 1 : i * group_size;
        ydata = data((i-1)*group_size + 1 : i * group_size)';
        
        b = bar(xdata, ydata, barcolor(i)); 
        b.EdgeColor = barcolor(i);
        hold on;
    end 
    title(titlestr);
    legend(feature_name);
end


end