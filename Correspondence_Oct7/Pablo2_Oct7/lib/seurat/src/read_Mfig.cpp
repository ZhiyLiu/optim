/********************************************************************************/
/*																				*/
/*  	File	:  read_Mfig.cpp													*/
/*																				*/
/*	Description:  class function for reading an Mfig (Slicefig or Quadfig)		*/
/*		into a Subdivsurf object												*/
/*																				*/
/*																				*/
/*	Project :  Seurat															*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	: 10. April 2000													*/
/*																				*/
/*	Modifications:																*/
/*		30. May -- added support for Parentset field in new Diatomgrid.			*/
/*			Allocate parallel array to dlist and initialize with list-index.	*/
/*		24. June 02 -- removed all .plist Parentset code						*/
/********************************************************************************/

#define D_SUBDIVSURF
#define D_XFERLIST
#define D_PSEUDOSET
#include "Shapedepend.h"
#include <stdio.h>

using namespace ThallCode;

/********************************************************************************/
/* skip_comments() -- standard file i/o function -- skip any line beginning	*/
/*	with a # sign---also skips any preliminary white-space, including	*/
/*	carriage-returns.							*/
/********************************************************************************/
#include <ctype.h>

static bool skip_comments(FILE *infile)
{
    char mychar;

    // Read past any line beginning with a # sign
    while ((mychar = getc(infile)) == '#' || isspace(mychar)) 
        while ((mychar = getc(infile)) != '\n');

    // If line doesn't begin with a # sign, push its first char back into stream
    ungetc(mychar, infile);

    return true;
}

/********************************************************************************/
/* read_mfig() -- load a figure from the given infile into the Slicefig		*/
/********************************************************************************/
bool Subdivsurf::read_mfig(FILE *infile)
{
/*
    char filetype[25];
    int numrows, numslices;
    DbVector3 invec;
    Quat inquat;
    double inphi, inr, inrho;
    int ftype;

    skip_comments(infile);
    fscanf(infile, "%s", filetype);
    if (strcmp(filetype, "QUADFIG") == 0)
        ftype = 0;
    else if (strcmp(filetype, "SLICEFIG") == 0)
        ftype = 1;
    else if (strcmp(filetype, "SLICEFIG2") == 0)
        ftype = 2;
    else {
        fprintf(stderr, "ERROR: input file not of allowable\n");
 	fprintf(stderr, "      but instead is %s\n", filetype);
 	return false;
    }
    // else continue with read

    skip_comments(infile);
    fscanf(infile, "%d %d", &numrows, &numslices);

    DbVector3 offset_fig;
    Quat q_fig;
    double scale_fig;

    // if new Slicefig2, input global spatial state, but ignore for mesh
    if (ftype == 2) {

        skip_comments(infile);
        fscanf(infile, "%lf %lf %lf", &offset_fig.X(), &offset_fig.Y(),
	                              &offset_fig.Z());
        fscanf(infile, "%lf %lf %lf %lf", &q_fig(X), &q_fig(Y), &q_fig(Z),
	                                  &q_fig(W));
        fscanf(infile, "%lf", &scale_fig);
    }

    // input figure center, magnification and frame, but ignore for mesh
    skip_comments(infile);
    fscanf(infile, "%lf %lf %lf", &offset_fig.X(), &offset_fig.Y(),
	                          &offset_fig.Z());
    fscanf(infile, "%lf %lf %lf %lf", &q_fig(X), &q_fig(Y), &q_fig(Z),
	                              &q_fig(W));
    fscanf(infile, "%lf", &scale_fig);

	// Assign Subdivgrid's fig_grid pointer to an array with the given Diatoms
    fig_grid->dlist = new Diatom[numrows*numslices];

    skip_comments(infile);
    // loop through dlist and enter Diatom values from file
    for (int dindex = 0; dindex < numrows*numslices; dindex++) {
		fscanf(infile, "%lf %lf %lf", &invec.X(), &invec.Y(), &invec.Z());
		fscanf(infile, "%lf %lf %lf %lf",
		       &inquat(X), &inquat(Y), &inquat(Z), &inquat(W));
		fscanf(infile, "%lf %lf %lf", &inphi, &inr, &inrho);

		fig_grid->dlist[dindex].p_val() = invec;
		fig_grid->dlist[dindex].q_val() = inquat;
		fig_grid->dlist[dindex].phi_val() = inphi;
		fig_grid->dlist[dindex].r_val() = inr;
		fig_grid->dlist[dindex].rho_val() = inrho;
		fig_grid->dlist[dindex].eta_val() = DEFAULT_ETA;
		fig_grid->dlist[dindex].clear_vflag();
    }
    fig_grid->rows = numrows;
    fig_grid->cols = numslices;

    for (int row = 0; row < numrows; row++)
        for (int slice = 0; slice < numslices; slice++) {
	    int ndex = slice*numrows + row;
	    if (row == 0 || row == numrows - 1 || slice == 0 || slice == numslices - 1) {
	        if ((row == 0 && slice == 0) || (row == 0 && slice == numslices - 1) 
		    		     	   || (row == numrows - 1 && slice == 0)
		    		|| (row == numrows - 1 && slice == numslices - 1)) 
		    fig_grid->dlist[ndex].Diatomtype(CORNER_M);
		else 
	            fig_grid->dlist[ndex].Diatomtype(EDGE_M);
            }
	    else 
	        fig_grid->dlist[ndex].Diatomtype(INTERNAL_M);
	}


	// True bad mesh
	fig_grid->true_badmesh();

    fclose(infile);
*/
    return true;   
}

/********************************************************************************/
/* modify_Diatom() -- use current state to modify the given Diatom		*/
/********************************************************************************/
void Subdivsurf::modify_Diatom(Diatom& atom)
{
/*
    // First, do local coordinate transformations
    DbVector3 vec(0, 0, 0);
    vec -= modelcoords.Objpos;

    atom.selfxlate_p(vec);

    modelcoords.Objorient.rotate_vec(atom.p_val());

    atom.selfrot_Diatom(modelcoords.Objorient);
    atom.selfzoom_r(modelcoords.Objscale);
    atom.p_val() *= modelcoords.Objscale;

    // Now, do global rotation, scaling, and translation
    globalloc.Objorient.rotate_vec(atom.p_val());
    atom.selfrot_Diatom(globalloc.Objorient);
    atom.selfzoom_r(globalloc.Objscale);
    atom.p_val() *= globalloc.Objscale;

    atom.selfxlate_p(globalloc.Objpos);
*/
}

/********************************************************************************/
/* read_mfig_output_m3d() -- load a figure from the given infile into the Slicefig		*/
/********************************************************************************/
bool Subdivsurf::read_mfig_output_m3d(FILE *infile)
{
/*
    char filetype[25];
	const char *outfilename = "testout.m3d";

	fprintf(stderr, "%s", outfilename);

    int numrows, numslices;
    DbVector3 invec;
    Quat inquat;
    double inphi, inr, inrho;
    int ftype;

    skip_comments(infile);
    fscanf(infile, "%s", filetype);
    if (strcmp(filetype, "QUADFIG") == 0)
        ftype = 0;
    else if (strcmp(filetype, "SLICEFIG") == 0)
        ftype = 1;
    else if (strcmp(filetype, "SLICEFIG2") == 0)
        ftype = 2;
    else {
        fprintf(stderr, "ERROR: input file not of allowable\n");
 	fprintf(stderr, "      but instead is %s\n", filetype);
 	return false;
    }
    // else continue with read

    skip_comments(infile);
    fscanf(infile, "%d %d", &numrows, &numslices);

    DbVector3 offset_fig;
    Quat q_fig;
    double scale_fig;

    // if new Slicefig2, input global spatial state
    if (ftype == 2) {

        skip_comments(infile);
        fscanf(infile, "%lf %lf %lf", &offset_fig.X(), &offset_fig.Y(),
	                              &offset_fig.Z());
        fscanf(infile, "%lf %lf %lf %lf", &q_fig(X), &q_fig(Y), &q_fig(Z),
	                                  &q_fig(W));
        fscanf(infile, "%lf", &scale_fig);

		globalloc.Objpos = offset_fig;
        globalloc.Objorient = q_fig;
        globalloc.Objscale = scale_fig;
    }

    // input figure center, magnification and frame, but ignore for mesh
    skip_comments(infile);
    fscanf(infile, "%lf %lf %lf", &offset_fig.X(), &offset_fig.Y(),
	                          &offset_fig.Z());
    fscanf(infile, "%lf %lf %lf %lf", &q_fig(X), &q_fig(Y), &q_fig(Z),
	                              &q_fig(W));
    fscanf(infile, "%lf", &scale_fig);
	modelcoords.Objpos = offset_fig;
    modelcoords.Objorient = q_fig;
    modelcoords.Objscale = scale_fig;

	// Assign Subdivgrid's fig_grid pointer to an array with the given Diatoms
    fig_grid->dlist = new Diatom[numrows*numslices];

    skip_comments(infile);
    // loop through dlist and enter Diatom values from file
    for (int dindex = 0; dindex < numrows*numslices; dindex++) {
		fscanf(infile, "%lf %lf %lf", &invec.X(), &invec.Y(), &invec.Z());
		fscanf(infile, "%lf %lf %lf %lf",
		       &inquat(X), &inquat(Y), &inquat(Z), &inquat(W));
		fscanf(infile, "%lf %lf %lf", &inphi, &inr, &inrho);

		fig_grid->dlist[dindex].p_val() = invec;
		fig_grid->dlist[dindex].q_val() = inquat;
		fig_grid->dlist[dindex].phi_val() = inphi;
		fig_grid->dlist[dindex].r_val() = inr;
		fig_grid->dlist[dindex].rho_val() = inrho;
		fig_grid->dlist[dindex].eta_val() = DEFAULT_ETA;
		fig_grid->dlist[dindex].clear_vflag();

		// Now, transform from model coords to global
		modify_Diatom(fig_grid->dlist[dindex]);
    }
    fig_grid->rows = numrows;
    fig_grid->cols = numslices;

    for (int row = 0; row < numrows; row++)
        for (int slice = 0; slice < numslices; slice++) {
			int ndex = slice*numrows + row;
			if (row == 0 || row == numrows - 1 || slice == 0 || slice == numslices - 1) {
				if ((row == 0 && slice == 0) || (row == 0 && slice == numslices - 1) 
					|| (row == numrows - 1 && slice == 0)
					|| (row == numrows - 1 && slice == numslices - 1)) 
					fig_grid->dlist[ndex].Diatomtype(CORNER_M);
				else 
					fig_grid->dlist[ndex].Diatomtype(EDGE_M);
            }
			else 
				fig_grid->dlist[ndex].Diatomtype(INTERNAL_M);
		}

    fclose(infile);
*/
    return true;   

}

/********************************************************************************/
/* output_m3d() -- output the Slicefig to .m3d file								*/
/********************************************************************************/
bool Subdivsurf::output_m3d(const char *outfile)
{
/*
	FILE *fp;
	double theta_degrees;
	Quat rotPIover2_x, work_q;

	// create quaternion to rotate frame 90deg between .m3d and .mfig standards
	DbVector3 one(1.0, 0.0, 0.0);
	rotPIover2_x.rot(one, M_PI/2.0);

	if ((fp = fopen(outfile, "r")) != NULL) {
		fprintf(stderr, "ERROR:  file %s exists already...move it!\n", outfile);
		fclose(fp);
		return false;
	}
	if ((fp = fopen(outfile, "w")) == NULL) {
        fprintf(stderr, "ERROR:  couldn't open file %s for write\n", outfile);
		return false;
	}
    else {
		fprintf(fp, "model {\n   figureCount = 1;\n   figure[0] {\n");

		int numrows = fig_grid->rows;
		int numcols = fig_grid->cols;
		fprintf(fp, "      numColumns = %d;\n      numRows = %d;\n",
			    numcols, numrows);
		fprintf(fp, "      type = QuadFigure;\n");

		Diatom workatom;

		for (int row = 0; row < numrows; row++)
			for (int col = 0; col < numcols; col++) {
				int ndex = col*numrows + row;

				workatom = fig_grid->dlist[ndex];
				theta_degrees = workatom.theta_val() * 180.0/M_PI;

				work_q = rotPIover2_x * workatom.q_val();

				fprintf(fp, "      primitive[%d][%d] {\n", row, col);
				fprintf(fp, "         qw = %f;\n", work_q(W));
				fprintf(fp, "         qx = %f;\n", work_q(X));
				fprintf(fp, "         qy = %f;\n", work_q(Y));
				fprintf(fp, "         qz = %f;\n", work_q(Z));
				fprintf(fp, "         r = %f;\n", workatom.r_val());
				fprintf(fp, "         selected = 0;\n");
				fprintf(fp, "         theta = %f;\n", theta_degrees);
				fprintf(fp, "         x = %f;\n", workatom.p_val().X());
				fprintf(fp, "         y = %f;\n", workatom.p_val().Y());
				fprintf(fp, "         z = %f;\n", workatom.p_val().Z());
				fprintf(fp, "      }\n");
			}
		fprintf(fp, "   }\n}\n");
	}

	fprintf(stderr, "file %s successfully written\n", outfile);
	fclose(fp);
*/
	return true;
}




