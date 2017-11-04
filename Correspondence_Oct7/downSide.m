matlabScriptFolder='/work/ltu/Correspondence_Oct7/';


addpath(strcat(matlabScriptFolder,'/shapeStat_downSide/Correspondence/'));

% Add these path to make sure the calculateCPNS_eigenvalue canworks
addpath(strcat(matlabScriptFolder,'/shapeStat_downSide/functions/'), strcat(matlabScriptFolder,'/shapeStat_downSide/ios/'));
addpath(strcat(matlabScriptFolder,'/shapeStat_downSide/'), strcat(matlabScriptFolder,'/shapeStat_downSide/CPNS/'));
addpath(strcat(matlabScriptFolder,'/shapeStat_downSide/CPNS/Sungkyu_additional_functions/subfunctions'));
  

geoEntropy(dataDir); 

quit;


