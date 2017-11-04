#ifndef M3DQUADINTERPOLATER_H
#define M3DQUADINTERPOLATER_H

#include "M3DSpoke.h"
#include "M3DFigure.h"
#include <vnl/vnl_double_1x3.h>

class M3DQuadInterpolater
{
    public:

    M3DQuadInterpolater();
    M3DQuadInterpolater(M3DFigure *figure);
    // Method that returns a single spoke at (u,v)
    //static M3DSpoke* interpolateSpoke(M3DFigure *figure, double u, double v, int side);

    // Methods that return spokes given a context - usually, you'll want to use these
    // With figureId, returns all interpolated spokes for given figureId
    // Adding atomId returns all spokes around specified atom
    // Adding spokeId returns all spokes around specified spoke

    //static vector <M3DSpoke*> getRelevantSpokes(M3DObject* targetObject, int level, int figureId);
    //static vector <M3DSpoke*> getRelevantSpokes(M3DObject* targetObject, int level, int figureId, int atomId);
    //static vector <M3DSpoke*> getRelevantSpokes(M3DObject* targetObject, int level, int figureId, int atomId, int spokeId);

    // Helper functions for hermite interpolation for medial surface

    // Function to draw interpolated boundary surface
    void displayInterpolatedBoundary(M3DFigure *figure, int level);

    // the old function
    M3DSpoke* interpolateQuadSpoke(M3DFigure *figure, double u, double v, int side);


    M3DSpoke interpolateQuadSpoke2(M3DFigure* figure, double u, double v, int side);



    M3DSpoke interpolateSpoke_mean(M3DFigure* figure, double u, double v, int side);
    M3DSpoke interpolateSpoke_hermite(M3DFigure* figure, double u, double v, int side);



    private:

    M3DFigure *_figure;

    double h1(double s) { return 2*(s * s * s) - 3*(s * s) + 1; }
    double h2(double s) { return -2*(s * s * s) + 3*(s * s); }
    double h3(double s) { return (s * s * s) - 2*(s * s) + s; }
    double h4(double s) { return (s * s * s) - (s * s); }

    double h1p(double s) { return 6*(s * s) - 6*(s); }
    double h2p(double s) { return -6*(s * s) + 6*(s); }
    double h3p(double s) { return 3*(s * s) - 4*(s) + 1; }
    double h4p(double s) { return 3*(s * s) - 2*(s); }
    };


#endif // M3DQUADINTERPOLATER_H
