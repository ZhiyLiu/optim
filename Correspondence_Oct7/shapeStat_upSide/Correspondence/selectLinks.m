%%% This function select links that link to the same object among
%%% population.
%%% Input: feature_mat: construct from all linking structures
%%% Input: dimPerLink: dimension of one link
%%% Input: index of property-linkTo in link definition
%%% Output: idList- the ids of selected links
%%% Output: links- selected links matrix, dimension: (nLinks*dimPerLink, nSamples)
function [idList, links] = selectLinks(feature_mat, dimPerLink, idOfLinkTo)
    [numRows, numCols] = size(feature_mat);
    numLinks = numRows / dimPerLink;
    
    idList = [];
    links = [];
    
    for i = 1:numLinks
        index_link_to = (i-1) * dimPerLink + idOfLinkTo;
        unique_val = unique(feature_mat(index_link_to, :));
        if size(unique_val,2) >1 || unique_val == -1
            % the link is variable among population, then dismiss this link
            % OR the link is always link to nothing, dismiss it
            continue;
        end
        
        links = [links; feature_mat((i-1) * dimPerLink + 1: i * dimPerLink, :)];
        idList = [idList; i-1];
    end
end