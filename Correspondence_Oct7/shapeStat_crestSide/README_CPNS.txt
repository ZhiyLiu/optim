% calculateCPNS.m:  Function to calculate the CPNS mean and modes for a
%                   collection of s-rep models in a folder. This program
%                   writes out the mean model and the eigenmode information 
%                   in a "*.m3d" file which can be read by Binary Pablo. 
%                   The program also outputs a "*.jpeg" image containing 
%                   the plot of the eigenvalues (first 20)

* RUNNING THE SCRIPT IN MATLAB

1.Add to matlab search paths the paths to '/functions', '/ios',
  '/output-plot', 'CPNS' and the directory where the five class directories
  (@QuadPrimitive, @QuadSymPrimitive, @SrepQuadPrimitive, @TubePrimitive,
  @TubeSymPrimitive) are.


2. Specifiy input parameters to the script
    
    dataDir     - input directory where training m3d files are stored
	
	dataDir => The folder where the s-reps are stored. This is also the
               folder where the output mean file and the eigenvalues plot
               is written

	outputFilename => 	The name of the output mean file. If this includes the path to a 
						directory, then the output mean file is written to that directory
						Otherwise it is written to the dataDir.  
 
	optionPNS => specify how CPNS should operate

	optionPNS 	= 0 => use SMALL circles or BIG circles depending upon analysis
				= 1 => always use SMALL circles
				= 2 => always use BIG circles
				
3. Check the values of control flags

	Ensure proper values are set for control flags inside the calculateCPNS() matlab program. 
	They are listed and described in Section 2 of the program itself.


4. Run the script. 

	calculateCPNS( dataDir, outputFilename, optionPNS )   

* OUTPUT FILES

<outputFilename>.m3d: 				mean m3d model with CPNS statistics (that can be read by Pablo)
<outputFilename>_CPNS_evals.jpg:	jpeg image containing a plot of Eigenvalues and cumulative Eigenvalue proportions

* CPNS THEORY AND DOCUMENTATION 

The PDF files inside the '/CPNS/documentation' contain files which describe the CPNS theory as well as how CPNS programs are designed in both Matlab
and C++ (Pablo).
'CPNS for s-reps.pdf': description of how CPNS can be applied to an s-rep (written by Sungkyu Jung)
'CPNS Implementation *.pdf': A series of notes on how CPNS has been implemented in Code (by Dibyendusekhar Goswami)
'CPNS Implementation 4: CPNS Software Design.pdf': Notes on how CPNS is implemented in C++ and Matlab (by Dibyendusekhar Goswami)

* SAMPLE DATA	

Some sample data (s-reps) to test the CPNS program is located in the folder "pablo_matlab/shapeStat/CPNS/sample_data"

* NOTE

For UNIX users, if the program has been copied from a windows computer, change the permission of the file "pablo_matlab/shapeStat/ios/linflat" to executable (x)



 
