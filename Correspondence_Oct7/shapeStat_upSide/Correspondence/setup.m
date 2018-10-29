% first set up path for references
baseDir = '/playpen/software/shapestat/shapeStat';

addpath(    baseDir, ...
            [baseDir filesep 'ios'], ...
            [baseDir filesep 'functions'],...
            [baseDir filesep 'PGA' filesep 'output-plot'],...
            [baseDir filesep 'PGA'], ...
            [baseDir filesep 'CPNS'], ...
            [baseDir filesep 'CPNS' filesep 'Sungkyu_additional_functions'],...
            [baseDir filesep 'CPNS' filesep 'Sungkyu_additional_functions' filesep 'subfunctions']...        
        );
addpath( genpath( fullfile( pwd, 'BatchAdjust' ) ) );
addpath( genpath( fullfile( pwd, 'General' ) ) );
addpath( genpath( fullfile( pwd, 'Smoothing' ) ) );