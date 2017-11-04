matlabScriptFolder='/work/ltu/Correspondence_Oct7/';


addpath(strcat(matlabScriptFolder,'/shapeStat_upSide/Correspondence/'));

% Add these path to make sure the calculateCPNS_eigenvalue canworks
addpath(strcat(matlabScriptFolder,'/shapeStat_upSide/functions/'), strcat(matlabScriptFolder,'/shapeStat_upSide/ios/'));
addpath(strcat(matlabScriptFolder,'/shapeStat_upSide/'), strcat(matlabScriptFolder,'/shapeStat_upSide/CPNS/'));
addpath(strcat(matlabScriptFolder,'/shapeStat_upSide/CPNS/Sungkyu_additional_functions/subfunctions'));
  

geoEntropy(dataDir); 

quit;


