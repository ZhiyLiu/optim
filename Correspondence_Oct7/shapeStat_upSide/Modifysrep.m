%setup;

% data_dir = fullfile( pwd );
% 
% [info, m3dData] = readM3d( data_dir, '*.m3d' );
% nRows = info.nRows ;
% nCols = info.nCols ;
% 
% 
%  for i = 1 : nRows,            
%      for j = 1 : nCols,
%         atomIndex = nCols * (i-1) + (j-1) + 1 ;
%         
%         if( ( i == 1 ) || ( i == nRows ) || ( j == 1 ) || ( j == nCols ) )
%             atomType = 1 ;
%         else
%             atomType = 0 ;
%         end
%         
%         if atomType == 1
%         %get position of end spoke
%         rVals = get( m3dData.atoms{1, atomIndex}, 'r' ) ;
%         posVals = get( m3dData.atoms{1,atomIndex}, 'pos' );
%         Uvals = get( m3dData.atoms{1, atomIndex}, 'U' );
%  
%         delta = 1.25*(rVals(1)+rVals(2))/2-rVals(3)
%         rVals(3)=rVals(3)+delta
%         posVals=posVals-delta*Uvals(3,:)'
%         
%         m3dData.atoms{1, atomIndex}=set(m3dData.atoms{1, atomIndex},'pos',posVals);
%         
%         m3dData.atoms{1, atomIndex}=set(m3dData.atoms{1, atomIndex},'r',rVals); 
%         end
%      end
%  end
[info, m3dData] = readM3d( '/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/', 'ellipsoid.m3d' );
UVals = get( m3dData.atoms{1, 26}, 'U' ) ;
theta = 0.3489;
rotMatrix = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
temp = UVals(:,1);
UVals(:, 1) =  rotMatrix*UVals(:,1);
m3dData.atoms{1, 26}=set(m3dData.atoms{1, 26},'U',UVals); 

% nRows = info.nRows; nCols = info.nCols;
% rVals = get( m3dData.atoms{1, 13}, 'r' ) ;
% rVals(2) = 2*rVals(2);
% m3dData.atoms{1, 13}=set(m3dData.atoms{1, 13},'r',rVals); 
%m3dData.atoms{1, 22}=set(m3dData.atoms{1, 22},'r',rVals); 
writeM3d('/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/modify_angle.m3d', info, -1, m3dData.atoms);
%[ neg_CPNS_data, ~ ] = getCPNSDataFromDirectory( data_neg_dir ,side) ;