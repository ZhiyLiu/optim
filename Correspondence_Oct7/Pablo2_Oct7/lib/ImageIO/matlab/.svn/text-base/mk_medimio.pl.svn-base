$pf = "C:\"/Program Files\"";

$mex = "$pf/MATLAB/R2008b/bin/mex.bat -v";

$includeDirs = " 
-I../include 
-I../../include 
-I../Templates 
-I../utilities/include 
-I../Analyze/include 
-I../Dicom/include 
-I../Meta/include 
-I../Gipl/include 
-I../PlanIm/include 
-I../planio/include 
-I../PlanIM/include 
-I../../m3d/include 
-I../../paul_code/include
";

$linkDirs = "
-L../../register/lib
-L../../seurat/lib
-L../../planes/lib
-L../../match/lib
-L../../m3d/lib
-L../../paul_code/lib
-L../../zlib/lib
-L../../ImageIO/lib
-L$pf/CLAPACK/lib
-LR:/libraries/flvw-1.0/lib
-LR:/libraries/fltk/lib
-L$pf/\"Microsoft Visual Studio\"/VC98/lib
";

$libs = "
-lcomctl32
-lWs2_32
-lregister
-lmatch
-lplanes
-lseurat
-lPlanIm
-lmixedImageIO
-lAnalyze
-lGipl
-lMeta
-lutilities
-lm3d
-lpaul_code
-lclapackMD
-lflvwd
-lfltkgld
-lfltkd
-lopengl32
-lglu32
-lzlib
";

$cmd = "$mex $includeDirs $linkDirs $libs readmedim.cpp";
$cmd =~ s/ *\n+ */ /g;

print $cmd, "\n";

system($cmd);
