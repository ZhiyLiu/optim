matlabScriptFolder='/work/ltu/Correspondence_Oct7/';


addpath(strcat(matlabScriptFolder,'/shapeStat_crestSide/Correspondence/'));

% Add these path to make sure the calculateCPNS_eigenvalue canworks
addpath(strcat(matlabScriptFolder,'/shapeStat_crestSide/functions/'), strcat(matlabScriptFolder,'/shapeStat_crestSide/ios/'));
addpath(strcat(matlabScriptFolder,'/shapeStat_crestSide/'), strcat(matlabScriptFolder,'/shapeStat_crestSide/CPNS/'));
addpath(strcat(matlabScriptFolder,'/shapeStat_crestSide/CPNS/Sungkyu_additional_functions/subfunctions'));
  

geoEntropy(dataDir); 

quit;


