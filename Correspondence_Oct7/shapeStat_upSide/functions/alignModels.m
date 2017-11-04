% ===================================================
% Align models
%
%   alignMode: 1 - translation only 
%              2 - translation and rotation
%              3 - translation, rotation and scaling
% ---------------------------------------------------

function [alignedDataStruct] = alignModels(alignMode, header, rawDataStruct)
   
if (~ismember(alignMode, 1:3))
    disp('  Not a valid alignment mode');
    disp('  ');
    return;
end

% ***************************
% Get parameters            *
% ***************************
nAtoms = header.nAtoms;
nSamps = header.nSamps;
atomMatrix = rawDataStruct.atoms;

aModes = {'translation', 'trans+rot', 'trans+rot/scaling'};

% IDs of sample atoms for visualization
visIDs = [1 floor(nAtoms/2) nAtoms];
nVis = length(visIDs);
colors = 'brg';

%-----------------
% Before alignment
%-----------------
figure(1), clf
if (alignMode == 1)
    subplot(121), hold on
else
    subplot(221), hold on
end

for i = 1:nVis        
    selectedAtoms = cVec( [atomMatrix{:, visIDs(i)}] );
    % (x, y, z)
    plot3(selectedAtoms(:,1), selectedAtoms(:,2), selectedAtoms(:,3), ['.' colors(i)] );
end
title('Before alignment');
axis([0 1 0 1]); view(-5, 12);
%-----------------

% ****************************************
% Align models                           *
% ****************************************

% [1] Translation alignment: center each model                                

% JJ => This is global alignment. With Tube, how? do it in the same way
Av = zeros(3, nSamps);

% 1) find the mean of all positions of medial atoms    
for i = 1:nSamps
    pos = zeros(nAtoms, 3);
    for j = 1:nAtoms 
        atomVec = cVec(atomMatrix{i,j});
        pos(j, :) = atomVec(1:3);
    end
    Av(:, i) = -mean(pos)';  % find the mean of all positions of medial atoms    
end

% 2) center each model and update atomMatrix and atomsInCell 
for i = 1:nAtoms
    for j = 1:nSamps
        atom = atomMatrix{j,i};
        atomMatrix{j,i} = set(atom, 'pos', get(atom, 'pos') + Av(:,j));
    end
end

%--------------------------------
% After aligning on translation
%--------------------------------
if (alignMode == 1)
    subplot(122), hold on
else
    subplot(222), hold on
end 
   
for i = 1:nVis        
    selectedAtoms = cVec( [atomMatrix{:, visIDs(i)}] );
    % (x, y, z)
    plot3(selectedAtoms(:,1), selectedAtoms(:,2), selectedAtoms(:,3) , ['.' colors(i)]  );
end
title('After translation');
axis([-0.5 .5 -.5 .5]); view(-5, 12);
%-----------------

if (alignMode == 2 || alignMode == 3)           % rotation or rotation/scaling
    Aqr = repmat([1;0;0;0;1], 1, nSamps);
    
    Fnow = alignmentMetric(atomMatrix);

    err = 1e3;
    epsilon = 3e-3; %1e-3
    cnt = 0;
    cntMAX = 1e3;

    h = waitbar(0, 'Aligning models...');
    while(err > epsilon && cnt < cntMAX)
        
        % Align each model to the mean of the rest  
        for i = 1:nSamps
            
            Ai = atomMatrix(i, :) ;  % cell array
            Arest = mrepMean( atomMatrix(setdiff(1:nSamps, i), :) , header); % cell array           
            opti = alignMreps(Ai, Arest, alignMode);
            Aqr(:,i) = [QuatProd(opti(1:4), Aqr(1:4,i)); opti(5)*Aqr(5)];
            atomMatrix(i, :) = applyMrepTransf(opti, Ai);
            
        end
    
        Fnext = alignmentMetric(atomMatrix);
        
        err = abs(Fnext - Fnow);
        Fnow = Fnext;
        cnt = cnt+1;
        waitbar(max(min(1, epsilon/err), cnt/cntMAX), h, ...
            ['Aligning models.  Current error: ' num2str(err)]);
    end
    close(h); drawnow;

    if (cnt == cntMAX & err > epsilon)  % results might not be correct
        uiwait(warndlg(['Max number of iterations attained.  Alignment result might' ...
            ' not be correct!'], 'Taking too long', 'modal'));
    end
    
    disp(['    Alignment metric is ' num2str(Fnow)]);

    %-----------------------------------------------------
    % After aligning on rotation or rotation with scaling
    %-----------------------------------------------------
    subplot(223), hold on
    for i = 1:nVis

        selectedAtoms = cVec( [atomMatrix{:, visIDs(i)}] );
          % (x, y, z)
        plot3(selectedAtoms(:,1), selectedAtoms(:,2), selectedAtoms(:,3), ['.' colors(i)]  );
    end
    
    title(['After ' aModes{alignMode}]);
    axis([-0.5 .5 -.5 .5]); view(-5, 12);
           
    subplot(224), hold on
    for i = 1:nAtoms   
        selectedAtoms = cVec( [atomMatrix{:, i}] );
        plot3(selectedAtoms(:,1), selectedAtoms(:,2), selectedAtoms(:,3), '.' );
    end
    title(['All atoms: after ' aModes{alignMode}]);
    axis([-0.5 .5 -.5 .5]); view(-5, 12);
    %-----------------------------------------------------
else
    Aqr = repmat([1;0;0;0;1], 1, nSamps);
end

% *****************
% Save results    *
% *****************
alignedDataStruct = rawDataStruct;
alignedDataStruct.atoms = atomMatrix;
alignedDataStruct.alignTransf = [Aqr; Av];                     
alignedDataStruct.alignMode = aModes{alignMode};

