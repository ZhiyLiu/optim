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


3. Run the script. 

	calculateCPNS( dataDir, outputFilename, optionPNS )   

*OUTPUT FILES

<outputFilename>.m3d: 				mean m3d model with CPNS statistics (that can be read by Pablo)
<outputFilename>_CPNS_evals.jpg:	jpeg image containing a plot of Eigenvalues and cumulative Eigenvalue proportions






 
