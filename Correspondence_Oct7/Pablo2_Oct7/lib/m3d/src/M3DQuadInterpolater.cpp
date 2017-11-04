#include <math.h>
#include "M3DQuadInterpolater.h"
#include <GL/gl.h>

#include <stdio.h>

M3DQuadInterpolater::M3DQuadInterpolater()
{
}

M3DQuadInterpolater::M3DQuadInterpolater(M3DFigure *figure)
{
    _figure = figure;
}



void M3DQuadInterpolater::displayInterpolatedBoundary(M3DFigure *figure, int level)
{
//    glBegin(GL_QUADS);
//    glNormal3d(0.0,0.0,1.0);
//    glVertex3d(0.0,0.0,0.0);
//    glVertex3d(1.0,0.0,0.0);
//    glVertex3d(0.0,1.0,0.0);
//    glVertex3d(1.0,1.0,0.0);
//    glEnd();
}

/* the old function */
M3DSpoke* M3DQuadInterpolater::interpolateQuadSpoke(M3DFigure *figure, double u, double v, int side) {
    int uBase = ( int ) floor ( u );
    int vBase = ( int ) floor ( v );

    M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );

    if ( uBase == tempFigure->getRowCount() - 1 )
    {
        uBase = uBase - 1;
    }
    if ( vBase == tempFigure->getColumnCount() - 1 )
    {
        vBase = vBase - 1;
    }
    // Get four corner atoms

    M3DPrimitive *atom11 = tempFigure->getPrimitivePtr ( uBase,vBase );
    M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*> ( atom11 );

    M3DPrimitive *atom21 = tempFigure->getPrimitivePtr ( uBase+1,vBase );
    M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*> ( atom21 );

    M3DPrimitive *atom12 = tempFigure->getPrimitivePtr ( uBase,vBase+1 );
    M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*> ( atom12 );

    M3DPrimitive *atom22 = tempFigure->getPrimitivePtr ( uBase+1,vBase+1 );
    M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*> ( atom22 );

    Vector3D x11 = quadAtom11->getX();
    Vector3D x21 = quadAtom21->getX();
    Vector3D x12 = quadAtom12->getX();
    Vector3D x22 = quadAtom22->getX();

    // Corner 11

    Vector3D l1 = x21 - x11;
    Vector3D l2 = x12 - x11;

    double mag1 = l1.norm();
    double mag2 = l2.norm();

    double dot = (l1 * l2) / (mag1 * mag2);
    double t11 = acos(dot);

    // Corner 12

    l1 = x22 - x12;
    l2 = x11 - x12;

    mag1 = l1.norm();
    mag2 = l2.norm();

    dot = (l1 * l2) / (mag1 * mag2);
    double t12 = acos(dot);

    // Corner 21

    l1 = x22 - x21;
    l2 = x11 - x21;

    mag1 = l1.norm();
    mag2 = l2.norm();

    dot = (l1 * l2) / (mag1 * mag2);
    double t21 = acos(dot);

    // Corner 22

    l1 = x21 - x22;
    l2 = x12 - x22;

    mag1 = l1.norm();
    mag2 = l2.norm();

    dot = (l1 * l2) / (mag1 * mag2);
    double t22 = acos(dot);

    M3DSpoke spoke;

    if ((t11 > 2.8) || (t12 > 2.8) || (t21 > 2.8) || (t22 > 2.8))
    {
        //std::cout << "u: " << uBase << ", v: " << vBase << ", using hermite interpolation" << std::endl;
        spoke = interpolateSpoke_hermite( figure, u, v, side );
    }
    else
    {
        //std::cout << "u: " << uBase << ", v: " << vBase << ", using spoke interpolation" << std::endl;
        spoke = interpolateSpoke_mean( figure, u, v, side );
    }

    M3DSpoke* interpolatedSpoke = new M3DSpoke(spoke.getX(), spoke.getU(), spoke.getR());

    return interpolatedSpoke;
}


/* fix the memory leak problem (without use new, and delete all the new [] variables). */
M3DSpoke M3DQuadInterpolater::interpolateQuadSpoke2(M3DFigure *figure, double u, double v, int side)
{
    int uBase = ( int ) floor ( u );
    int vBase = ( int ) floor ( v );

    M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );

    if ( uBase == tempFigure->getRowCount() - 1 )
    {
        uBase = uBase - 1;
    }
    if ( vBase == tempFigure->getColumnCount() - 1 )
    {
        vBase = vBase - 1;
    }
    // Get four corner atoms

    M3DPrimitive *atom11 = tempFigure->getPrimitivePtr ( uBase,vBase );
    M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*> ( atom11 );

    M3DPrimitive *atom21 = tempFigure->getPrimitivePtr ( uBase+1,vBase );
    M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*> ( atom21 );

    M3DPrimitive *atom12 = tempFigure->getPrimitivePtr ( uBase,vBase+1 );
    M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*> ( atom12 );

    M3DPrimitive *atom22 = tempFigure->getPrimitivePtr ( uBase+1,vBase+1 );
    M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*> ( atom22 );

    Vector3D x11 = quadAtom11->getX();
    Vector3D x21 = quadAtom21->getX();
    Vector3D x12 = quadAtom12->getX();
    Vector3D x22 = quadAtom22->getX();

    // Corner 11

    Vector3D l1 = x21 - x11;
    Vector3D l2 = x12 - x11;

    double mag1 = l1.norm();
    double mag2 = l2.norm();

    double dot = (l1 * l2) / (mag1 * mag2);
    double t11 = acos(dot);

    // Corner 12

    l1 = x22 - x12;
    l2 = x11 - x12;

    mag1 = l1.norm();
    mag2 = l2.norm();

    dot = (l1 * l2) / (mag1 * mag2);
    double t12 = acos(dot);

    // Corner 21

    l1 = x22 - x21;
    l2 = x11 - x21;

    mag1 = l1.norm();
    mag2 = l2.norm();

    dot = (l1 * l2) / (mag1 * mag2);
    double t21 = acos(dot);

    // Corner 22

    l1 = x21 - x22;
    l2 = x12 - x22;

    mag1 = l1.norm();
    mag2 = l2.norm();

    dot = (l1 * l2) / (mag1 * mag2);
    double t22 = acos(dot);

    //M3DSpoke* spoke;
    M3DSpoke spoke;

    if ((t11 > 2.8) || (t12 > 2.8) || (t21 > 2.8) || (t22 > 2.8))
    {
        //std::cout << "u: " << uBase << ", v: " << vBase << ", using hermite interpolation" << std::endl;
        spoke = interpolateSpoke_hermite( figure, u, v, side );
    }
    else
    {
        //std::cout << "u: " << uBase << ", v: " << vBase << ", using spoke interpolation" << std::endl;
        spoke = interpolateSpoke_mean( figure, u, v, side );
    }

    return spoke;
}



M3DSpoke M3DQuadInterpolater::interpolateSpoke_mean( M3DFigure *figure, double u, double v, int side)
{
    //std::cout << "interp mean" << std::endl;
    int uBase = ( int ) floor ( u );
    int vBase = ( int ) floor ( v );

    u = u - uBase;
    v = v - vBase;

    M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );

    if ( uBase == tempFigure->getRowCount() - 1 )
    {
        uBase = uBase - 1;
        u=1;
    }
    if ( vBase == tempFigure->getColumnCount() - 1 )
    {
        vBase = vBase - 1;
        v=1;
    }
    // Get four corner atoms

    M3DPrimitive *atom11 = tempFigure->getPrimitivePtr ( uBase,vBase );
    M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*> ( atom11 );

    M3DPrimitive *atom21 = tempFigure->getPrimitivePtr ( uBase+1,vBase );
    M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*> ( atom21 );

    M3DPrimitive *atom12 = tempFigure->getPrimitivePtr ( uBase,vBase+1 );
    M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*> ( atom12 );

    M3DPrimitive *atom22 = tempFigure->getPrimitivePtr ( uBase+1,vBase+1 );
    M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*> ( atom22 );

    Vector3D x11 = quadAtom11->getX();
    Vector3D x21 = quadAtom21->getX();
    Vector3D x12 = quadAtom12->getX();
    Vector3D x22 = quadAtom22->getX();

    //std::cout << "x11: " << x11.getX() << ", " << x11.getY() << ", " << x11.getZ() << std::endl;
    //std::cout << "x21: " << x21.getX() << ", " << x21.getY() << ", " << x21.getZ() << std::endl;
    //std::cout << "x12: " << x12.getX() << ", " << x12.getY() << ", " << x12.getZ() << std::endl;
    //std::cout << "x22: " << x22.getX() << ", " << x22.getY() << ", " << x22.getZ() << std::endl;

    double r11, r12, r21, r22, ru11, ru12, ru21, ru22, rv11, rv12, rv21, rv22;

    Vector3D U0_11, U0_12, U0_21, U0_22, B0_11, B0_12, B0_21, B0_22, S0_11, S0_12,S0_21,S0_22;

    if ( side == 0 )
    {
        U0_11 = quadAtom11->getU0();
        U0_21 = quadAtom21->getU0();
        U0_12 = quadAtom12->getU0();
        U0_22 = quadAtom22->getU0();

        B0_11 = quadAtom11->getX() + quadAtom11->getR0() * quadAtom11->getU0();
        B0_21 = quadAtom21->getX() + quadAtom21->getR0() * quadAtom21->getU0();
        B0_12 = quadAtom12->getX() + quadAtom12->getR0() * quadAtom12->getU0();
        B0_22 = quadAtom22->getX() + quadAtom22->getR0() * quadAtom22->getU0();

        S0_11 = U0_11 * quadAtom11->getR0();
        S0_21 = U0_21 * quadAtom21->getR0();
        S0_12 = U0_12 * quadAtom12->getR0();
        S0_22 = U0_22 * quadAtom22->getR0();

        r11 = quadAtom11->getR0();
        r21 = quadAtom21->getR0();
        r12 = quadAtom12->getR0();
        r22 = quadAtom22->getR0();

//		std::cout << "r11: " << quadAtom11->getR0() << std::endl;
//		std::cout << "r21: " << quadAtom21->getR0() << std::endl;
//		std::cout << "r12: " << quadAtom12->getR0() << std::endl;
//		std::cout << "r22: " << quadAtom22->getR0() << std::endl;
    }
    else
    {
        U0_11 = quadAtom11->getU1();
        U0_21 = quadAtom21->getU1();
        U0_12 = quadAtom12->getU1();
        U0_22 = quadAtom22->getU1();

        B0_11 = quadAtom11->getX() + quadAtom11->getR1() * quadAtom11->getU1();
        B0_21 = quadAtom21->getX() + quadAtom21->getR1() * quadAtom21->getU1();
        B0_12 = quadAtom12->getX() + quadAtom12->getR1() * quadAtom12->getU1();
        B0_22 = quadAtom22->getX() + quadAtom22->getR1() * quadAtom22->getU1();

        S0_11 = U0_11 * quadAtom11->getR1();
        S0_21 = U0_21 * quadAtom21->getR1();
        S0_12 = U0_12 * quadAtom12->getR1();
        S0_22 = U0_22 * quadAtom22->getR1();

        r11 = quadAtom11->getR1();
        r21 = quadAtom21->getR1();
        r12 = quadAtom12->getR1();
        r22 = quadAtom22->getR1();
    }

    //std::cout << "U0_11: " << U0_11.getX() << ", " << U0_11.getY() << ", " << U0_11.getZ() << std::endl;
    //std::cout << "U0_21: " << U0_21.getX() << ", " << U0_21.getY() << ", " << U0_21.getZ() << std::endl;
    //std::cout << "U0_12: " << U0_12.getX() << ", " << U0_12.getY() << ", " << U0_12.getZ() << std::endl;
    //std::cout << "U0_22: " << U0_22.getX() << ", " << U0_22.getY() << ", " << U0_22.getZ() << std::endl;

    Vector3D du11, du21, du12, du22, dv11, dv21, dv12, dv22;

    Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;
    Vector3D dB0du11, dB0dv11, dB0du21, dB0dv21, dB0du12, dB0dv12, dB0du22, dB0dv22;
    Vector3D dS0du11, dS0dv11, dS0du21, dS0dv21, dS0du12, dS0dv12, dS0du22, dS0dv22;

    dB0du11 = B0_21 - B0_11;
    dB0dv11 = B0_12 - B0_11;

    dB0du21 = B0_21 - B0_11;
    dB0dv21 = B0_22 - B0_21;

    dB0du12 = B0_22 - B0_12;
    dB0dv12 = B0_12 - B0_11;

    dB0du22 = B0_22 - B0_12;
    dB0dv22 = B0_22 - B0_21;

    // Get derivatives for each corner, starting at 11, then 21, 12, and 22

    // Corner 11
    // du

    if ( uBase == 0 )
    {
        du11 = x21 - x11;
        ru11 = r21 - r11;
        dU0du11 = U0_21 - U0_11;
        dS0du11 = S0_21 - S0_11;
        dB0du11 = B0_21 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase-1, vBase ) );
        du11 = ( x21 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            ru11 = ( r21 - tempAtom->getR0() ) / 2;
            dU0du11 = ( U0_21 - tempAtom->getU0() ) / 2;
            dS0du11 = ( S0_21 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dB0du11 = ( B0_21 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            ru11 = ( r21 - tempAtom->getR1() ) / 2;
            dU0du11 = ( U0_21 - tempAtom->getU1() ) / 2;
            dS0du11 = ( S0_21 - ( tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
            dB0du11 = ( B0_21 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }

        //dU0du11.print();
    }

    //dv

    if ( vBase == 0 )
    {
        dv11 = x12 - x11;
        rv11 = r12 - r11;
        dU0dv11 = U0_12 - U0_11;
        dS0dv11 = S0_12 - S0_11;
        dB0dv11 = B0_12 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase, vBase-1 ) );
        dv11 = ( x12 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            rv11 = ( r12 - tempAtom->getR0() ) / 2;
            dU0dv11 = ( U0_12 - tempAtom->getU0() ) / 2;
            dS0dv11 = ( U0_12 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dB0dv11 = ( B0_12 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            rv11 = ( r12 - tempAtom->getR1() ) / 2;
            dU0dv11 = ( U0_12 - tempAtom->getU1() ) / 2;
            dS0dv11 = ( U0_12 - ( tempAtom->getR0() * tempAtom->getU1() ) ) / 2;
            dB0dv11 = ( B0_12 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }
    }

    // Corner 21

    // du

    //std::cout << tempFigure->getRowCount() - 1 << std::endl;
    //std::cout << tempFigure->getColumnCount() - 1  << std::endl;
    if ( uBase + 1 == tempFigure->getRowCount() - 1 )
    {
        du21 = x21 - x11;
        ru21 = r21 - r11;
        dS0du21 = S0_21 - S0_11;
        dU0du21 = U0_21 - U0_11;
        dB0du21 = B0_21 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+2, vBase ) );
        du21 = ( tempAtom->getX() - x11 ) / 2;

        if ( side == 0 )
        {
            ru21 = ( tempAtom->getR0() - r11 ) / 2;
            dS0du21 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_11 ) / 2;
            dU0du21 = ( tempAtom->getU0() - U0_11 ) / 2;
            dB0du21 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_11 ) / 2;
        }
        else
        {
            ru21 = ( tempAtom->getR1() - r11 ) / 2;
            dS0du21 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_11 ) / 2;
            dU0du21 = ( tempAtom->getU1() - U0_11 ) / 2;
            dB0du21 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_11 ) / 2;
        }
    }

    //dv

    if ( vBase == 0 )
    {
        dv21 = x22 - x21;
        rv21 = r22 - r21;
        dS0dv21 = S0_22 - S0_21;
        dU0dv21 = U0_22 - U0_21;
        dB0dv21 = B0_22 - B0_21;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+1, vBase-1 ) );
        dv21 = ( x22 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            rv21 = ( r22 - tempAtom->getR0() ) / 2;
            dS0dv21 = ( S0_22 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dU0dv21 = ( U0_22 - tempAtom->getU0() ) / 2;
            dB0dv21 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            rv21 = ( r22 - tempAtom->getR1() ) / 2;
            dS0dv21 = ( S0_22 - ( tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
            dU0dv21 = ( U0_22 - tempAtom->getU1() ) / 2;
            dB0dv21 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }
    }

    // Corner 12

    // du

    if ( uBase == 0 )
    {
        du12 = x22 - x12;
        ru12 = r22 - r12;
        dS0du12 = S0_22 - S0_12;
        dU0du12 = U0_22 - U0_12;
        dB0du12 = B0_22 - B0_12;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase-1,vBase+1 ) );
        du12 = ( x22 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            ru12 = ( r22 - tempAtom->getR0() ) / 2;
            dS0du12 = ( S0_22 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dU0du12 = ( U0_22 - tempAtom->getU0() ) / 2;
            dB0du12 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            ru12 = ( r22 - tempAtom->getR1() ) / 2;
            dS0du12 = ( S0_22 - ( tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
            dU0du12 = ( U0_22 - tempAtom->getU1() ) / 2;
            dB0du12 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }
    }

    if ( vBase + 1 == tempFigure->getColumnCount() - 1 )
    {
        dv12 = x12 - x11;
        rv12 = r12 - r11;
        dS0dv12 = S0_12 - S0_11;
        dU0dv12 = U0_12 - U0_11;
        dB0dv12 = B0_12 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase, vBase+2 ) );
        dv12 = ( tempAtom->getX() - x11 ) / 2;

        if ( side == 0 )
        {
            rv12 = ( tempAtom->getR0() - r11 ) / 2;
            dS0dv12 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_11 ) / 2;
            dU0dv12 = ( tempAtom->getU0() - U0_11 ) / 2;
            dB0dv12 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_11 ) / 2;
        }
        else
        {
            rv12 = ( tempAtom->getR1() - r11 ) / 2;
            dS0dv12 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_11 ) / 2;
            dU0dv12 = ( tempAtom->getU1() - U0_11 ) / 2;
            dB0dv12 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_11 ) / 2;
        }
    }

    // Corner 22

    if ( uBase + 1 == tempFigure->getRowCount() - 1 )
    {
        du22 = x22 - x12;
        ru22 = r22 - r12;
        dS0du22 = S0_22 - S0_12;
        dU0du22 = U0_22 - U0_12;
        dB0du22 = B0_22 - B0_12;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+2,vBase+1 ) );
        du22 = ( tempAtom->getX() - x12 ) / 2;

        if ( side == 0 )
        {
            ru22 = ( tempAtom->getR0() - r12 ) / 2;
            dS0du22 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_12 ) / 2;
            dU0du22 = ( tempAtom->getU0() - U0_12 ) / 2;
            dB0du22 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_12 ) / 2;
        }
        else
        {
            ru22 = ( tempAtom->getR1() - r12 ) / 2;
            dS0du22 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_12 ) / 2;
            dU0du22 = ( tempAtom->getU1() - U0_12 ) / 2;
            dB0du22 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_12 ) / 2;
        }
    }

    if ( vBase + 1 == tempFigure->getColumnCount() - 1 )
    {
        dv22 = x22 - x21;
        rv22 = r22 - r21;
        dS0dv22 = S0_22 - S0_21;
        dU0dv22 = U0_22 - U0_21;
        dB0dv22 = B0_22 - B0_21;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+1,vBase+2 ) );
        dv22 = ( tempAtom->getX() - x21 ) / 2;

        if ( side == 0 )
        {
            rv22 = ( tempAtom->getR0() - r22 ) / 2;
            dS0dv22 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_21 ) / 2;
            dU0dv22 = ( tempAtom->getU0() - U0_21 ) / 2;
            dB0dv22 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_21 ) / 2;
        }
        else
        {
            rv22 = ( tempAtom->getR1() - r22 ) / 2;
            dS0dv22 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_21 ) / 2;
            dU0dv22 = ( tempAtom->getU1() - U0_21 ) / 2;
            dB0dv22 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_21 ) / 2;
        }
    }

    du11 = x21 - x11;
    dv11 = x12 - x11;
    dU0du11 = U0_21 - U0_11;
    dU0dv11 = U0_12 - U0_11;
    dS0du11 = S0_21 - S0_11;
    dS0dv11 = S0_12 - S0_11;
    ru11 = r21 - r11;
    rv11 = r12 - r11;
//
//
    du21 = x21 - x11;
    dv21 = x22 - x21;
    dU0du21 = U0_21 - U0_11;
    dU0dv21 = U0_22 - U0_21;
    dS0du21 = S0_21 - S0_11;
    dS0dv21 = S0_22 - S0_21;
    ru21 = r21 - r11;
    rv21 = r22 - r21;
//
    du12 = x22 - x12;
    dv12 = x12 - x11;
    dU0du12 = U0_22 - U0_12;
    dU0dv12 = U0_12 - U0_11;
    dS0du12 = S0_22 - S0_12;
    dS0dv12 = S0_12 - S0_11;
    ru12 = r22 - r12;
    rv12 = r12 - r11;
//
    du22 = x22 - x12;
    dv22 = x22 - x21;
    dU0du22 = U0_22 - U0_12;
    dU0dv22 = U0_22 - U0_21;
    dS0du22 = S0_22 - S0_12;
    dS0dv22 = S0_22 - S0_21;
    ru22 = r22 - r12;
    rv22 = r22 - r21;



//	Vector3D du11 = (x21 - x01) / 2;
//	Vector3D du21 = (x31 - x11) / 2;
//	Vector3D du12 = (x22 - x02) / 2;
//	Vector3D du22 = (x32 - x12) / 2;
//
//	Vector3D dv11 = (x12 - x10) / 2;
//	Vector3D dv21 = (x22 - x20) / 2;
//	Vector3D dv12 = (x13 - x11) / 2;
//	Vector3D dv22 = (x23 - x21) / 2;

//	std::cout << "dU0du11: " << dU0du11.getX() << ", " << dU0du11.getY() << ", " << dU0du11.getZ() << std::endl;
//	std::cout << "dU0du21: " << dU0du21.getX() << ", " << dU0du21.getY() << ", " << dU0du21.getZ() << std::endl;
//	std::cout << "dU0du12: " << dU0du12.getX() << ", " << dU0du12.getY() << ", " << dU0du12.getZ() << std::endl;
//	std::cout << "dU0du22: " << dU0du22.getX() << ", " << dU0du22.getY() << ", " << dU0du22.getZ() << std::endl;
//
//	std::cout << "dU0dv11: " << dU0dv11.getX() << ", " << dU0dv11.getY() << ", " << dU0dv11.getZ() << std::endl;
//	std::cout << "dU0dv21: " << dU0dv21.getX() << ", " << dU0dv21.getY() << ", " << dU0dv21.getZ() << std::endl;
//	std::cout << "dU0dv12: " << dU0dv12.getX() << ", " << dU0dv12.getY() << ", " << dU0dv12.getZ() << std::endl;
//	std::cout << "dU0dv22: " << dU0dv22.getX() << ", " << dU0dv22.getY() << ", " << dU0dv22.getZ() << std::endl;

    // Get unit normals for the four corners and project derivatives on to the tangent plane
    Vector3D n11 = quadAtom11->getU0() - quadAtom11->getU1();
    n11.normalize();
    Vector3D n21 = quadAtom21->getU0() - quadAtom21->getU1();
    n21.normalize();
    Vector3D n12 = quadAtom12->getU0() - quadAtom12->getU1();
    n12.normalize();
    Vector3D n22 = quadAtom22->getU0() - quadAtom22->getU1();
    n22.normalize();

//	n11 = du11.cross(dv11);
//	n11.normalize();
//	n12 = du12.cross(dv12);
//	n12.normalize();
//	n21 = du21.cross(dv21);
//	n21.normalize();
//	n22 = du22.cross(dv22);
//	n22.normalize();

    Vector3D du11t = du11 - ( du11 * n11 ) * n11;
    Vector3D du21t = du21 - ( du21 * n21 ) * n21;
    Vector3D du12t = du12 - ( du12 * n12 ) * n12;
    Vector3D du22t = du22 - ( du22 * n22 ) * n22;

    Vector3D dv11t = dv11 - ( dv11 * n11 ) * n11;
    Vector3D dv21t = dv21 - ( dv21 * n21 ) * n21;
    Vector3D dv12t = dv12 - ( dv12 * n12 ) * n12;
    Vector3D dv22t = dv22 - ( dv22 * n22 ) * n22;

    // Build matrices for hermite interpolation of medial sheet
    double hx[16] = { x11.getX(), x21.getX(), du11t.getX(), du21t.getX(),
                      x12.getX(), x22.getX(), du12t.getX(), du22t.getX(), dv11t.getX(),
                      dv21t.getX(), 0, 0, dv12t.getX(), dv22t.getX(), 0, 0
                    };

    double hy[16] = { x11.getY(), x21.getY(), du11t.getY(), du21t.getY(),
                      x12.getY(), x22.getY(), du12t.getY(), du22t.getY(), dv11t.getY(),
                      dv21t.getY(), 0, 0, dv12t.getY(), dv22t.getY(), 0, 0
                    };

    double hz[16] = { x11.getZ(), x21.getZ(), du11t.getZ(), du21t.getZ(),
                      x12.getZ(), x22.getZ(),du12t.getZ(), du22t.getZ(), dv11t.getZ(),
                      dv21t.getZ(), 0, 0, dv12t.getZ(), dv22t.getZ(), 0, 0
                    };

    double r_hermite[16] = { r11, r21, ru11, ru21, r12, r22, ru12, ru22, rv11, rv21, 0, 0, rv12, rv22, 0, 0};

    Matrix r_hermite_mat = Matrix ( 4,4,r_hermite, true );

    Matrix hxmat = Matrix ( 4, 4, hx, true );
    Matrix hymat = Matrix ( 4, 4, hy, true );
    Matrix hzmat = Matrix ( 4, 4, hz, true );

    //U0_11.print();

    double hu[4] = { h1 ( u ), h2 ( u ), h3 ( u ), h4 ( u ) };
    double hv[4] = { h1 ( v ), h2 ( v ), h3 ( v ), h4 ( v ) };
    Matrix humat = Matrix ( 1, 4, hu, true );
    Matrix hvmat = Matrix ( 4, 1, hv, true );

    Matrix xn = humat * hxmat * hvmat;
    Matrix yn = humat * hymat * hvmat;
    Matrix zn = humat * hzmat * hvmat;

    double Bx[16] = { B0_11.getX(), B0_21.getX(), dB0du11.getX(), dB0du21.getX(),
                      B0_12.getX(), B0_22.getX(), dB0du12.getX(), dB0du22.getX(), dB0dv11.getX(),
                      dB0dv21.getX(), 0, 0, dB0dv12.getX(), dB0dv22.getX(), 0, 0
                    };

    double By[16] = { B0_11.getY(), B0_21.getY(), dB0du11.getY(), dB0du21.getY(),
                      B0_12.getY(), B0_22.getY(), dB0du12.getY(), dB0du22.getY(), dB0dv11.getY(),
                      dB0dv21.getY(), 0, 0, dB0dv12.getY(), dB0dv22.getY(), 0, 0
                    };

    double Bz[16] = { B0_11.getZ(), B0_21.getZ(), dB0du11.getZ(), dB0du21.getZ(),
                      B0_12.getZ(), B0_22.getZ(), dB0du12.getZ(), dB0du22.getZ(), dB0dv11.getZ(),
                      dB0dv21.getZ(), 0, 0, dB0dv12.getZ(), dB0dv22.getZ(), 0, 0
                    };

    Matrix bxmat = Matrix ( 4, 4, Bx, true );
    Matrix bymat = Matrix ( 4, 4, By, true );
    Matrix bzmat = Matrix ( 4, 4, Bz, true );

    Matrix Bxn = humat * bxmat * hvmat;
    Matrix Byn = humat * bymat * hvmat;
    Matrix Bzn = humat * bzmat * hvmat;


    // Calculate sRad matrices using finite differences
    /*Vector3D dU0du11 = (quadAtom21->getU0() - quadAtom01->getU0()) / 2;
    Vector3D dU0dv11 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

    Vector3D dU0du21 = (quadAtom31->getU0() - quadAtom11->getU0()) / 2;
    Vector3D dU0dv21 = (quadAtom22->getU0() - quadAtom20->getU0()) / 2;

    Vector3D dU0du12 = (quadAtom22->getU0() - quadAtom02->getU0()) / 2;
    Vector3D dU0dv12 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

    Vector3D dU0du22 = (quadAtom32->getU0() - quadAtom12->getU0()) / 2;
    Vector3D dU0dv22 = (quadAtom23->getU0() - quadAtom21->getU0()) / 2;*/

    /*Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;

    if (side == 0)
    {
        dU0du11 = quadAtom21->getU0() - quadAtom11->getU0();
        dU0dv11 = quadAtom12->getU0() - quadAtom11->getU0();

        dU0du21 = quadAtom21->getU0() - quadAtom11->getU0();
        dU0dv21 = quadAtom22->getU0() - quadAtom21->getU0();

        dU0du12 = quadAtom22->getU0() - quadAtom12->getU0();
        dU0dv12 = quadAtom12->getU0() - quadAtom11->getU0();

        dU0du22 = quadAtom22->getU0() - quadAtom12->getU0();
        dU0dv22 = quadAtom22->getU0() - quadAtom21->getU0();
    }
    else
    {
        dU0du11 = quadAtom21->getU1() - quadAtom11->getU1();
        dU0dv11 = quadAtom12->getU1() - quadAtom11->getU1();

        dU0du21 = quadAtom21->getU1() - quadAtom11->getU1();
        dU0dv21 = quadAtom22->getU1() - quadAtom21->getU1();

        dU0du12 = quadAtom22->getU1() - quadAtom12->getU1();
        dU0dv12 = quadAtom12->getU1() - quadAtom11->getU1();

        dU0du22 = quadAtom22->getU1() - quadAtom12->getU1();
        dU0dv22 = quadAtom22->getU1() - quadAtom21->getU1();
    }*/

    // These form the non-orthogonal medial coordinate system for the spoke derivatives

    //Vector3D U0_11, U0_12, U0_21, U0_22;

    //if (side == 0)
    //{
    //	U0_11 = quadAtom11->getU0();
    //	U0_21 = quadAtom21->getU0();
    //	U0_12 = quadAtom12->getU0();
    //	U0_22 = quadAtom22->getU0();
    //}
    //else
    //{
    //	U0_11 = quadAtom11->getU1();
    //	U0_21 = quadAtom21->getU1();
    //	U0_12 = quadAtom12->getU1();
    //	U0_22 = quadAtom22->getU1();
    //}
    Vector3D du11p = n11.cross ( du11 );
    Vector3D du12p = n12.cross ( du12 );
    Vector3D du21p = n21.cross ( du21 );
    Vector3D du22p = n22.cross ( du22 );

    Vector3D du11tp = n11.cross ( du11t );
    Vector3D du12tp = n12.cross ( du12t );
    Vector3D du21tp = n21.cross ( du21t );
    Vector3D du22tp = n22.cross ( du22t );

    Vector3D du11t_norm = du11t / du11t.norm();
    Vector3D du12t_norm = du12t / du12t.norm();
    Vector3D du21t_norm = du21t / du21t.norm();
    Vector3D du22t_norm = du22t / du22t.norm();

    Vector3D du11pt_norm = du11tp / du11tp.norm();
    Vector3D du12pt_norm = du12tp / du12tp.norm();
    Vector3D du21pt_norm = du21tp / du21tp.norm();
    Vector3D du22pt_norm = du22tp / du22tp.norm();

    double basis11[9] = { du11t.getX(), du11t.getY(), du11t.getZ(), dv11t.getX(), dv11t.getY(), dv11t.getZ(),
                          U0_11.getX(), U0_11.getY(), U0_11.getZ()
                        };

    double basis21[9] = { du21t.getX(), du21t.getY(), du21t.getZ(), dv21t.getX(), dv21t.getY(), dv21t.getZ(),
                          U0_21.getX(), U0_21.getY(), U0_21.getZ()
                        };

    double basis12[9] = { du12t.getX(), du12t.getY(), du12t.getZ(), dv12t.getX(), dv12t.getY(), dv12t.getZ(),
                          U0_12.getX(), U0_12.getY(), U0_12.getZ()
                        };

    double basis22[9] = { du22t.getX(), du22t.getY(), du22t.getZ(), dv22t.getX(), dv22t.getY(), dv22t.getZ(),
                          U0_22.getX(), U0_22.getY(), U0_22.getZ()
                        };

    /*double basis11[9] = { du11t_norm.getX(), du11t_norm.getY(), du11t_norm.getZ(), du11pt_norm.getX(), du11pt_norm.getY(), du11pt_norm.getZ(),
        U0_11.getX(), U0_11.getY(), U0_11.getZ() };

    double basis21[9] = { du21t_norm.getX(), du21t_norm.getY(), du21t_norm.getZ(), du21pt_norm.getX(), du21pt_norm.getY(), du21pt_norm.getZ(),
        U0_21.getX(), U0_21.getY(), U0_21.getZ() };

    double basis12[9] = { du12t_norm.getX(), du12t_norm.getY(), du12t_norm.getZ(), du12pt_norm.getX(), du12pt_norm.getY(), du12pt_norm.getZ(),
        U0_12.getX(), U0_12.getY(), U0_12.getZ() };

    double basis22[9] = { du22t_norm.getX(), du22t_norm.getY(), du22t_norm.getZ(), du22pt_norm.getX(), du22pt_norm.getY(), du22pt_norm.getZ(),
        U0_22.getX(), U0_22.getY(), U0_22.getZ() };*/



    Matrix C11 = Matrix ( 3, 3, basis11, true );
    Matrix C21 = Matrix ( 3, 3, basis21, true );
    Matrix C12 = Matrix ( 3, 3, basis12, true );
    Matrix C22 = Matrix ( 3, 3, basis22, true );

    Matrix Ci11, Ci21, Ci12, Ci22;

    C11.inverse ( Ci11 );
    C21.inverse ( Ci21 );
    C12.inverse ( Ci12 );
    C22.inverse ( Ci22 );

    double dU0du11_coeffs[3] = { dU0du11.getX(), dU0du11.getY(), dU0du11.getZ() };
    double dU0dv11_coeffs[3] = { dU0dv11.getX(), dU0dv11.getY(), dU0dv11.getZ() };
    double dU0du21_coeffs[3] = { dU0du21.getX(), dU0du21.getY(), dU0du21.getZ() };
    double dU0dv21_coeffs[3] = { dU0dv21.getX(), dU0dv21.getY(), dU0dv21.getZ() };
    double dU0du12_coeffs[3] = { dU0du12.getX(), dU0du12.getY(), dU0du12.getZ() };
    double dU0dv12_coeffs[3] = { dU0dv12.getX(), dU0dv12.getY(), dU0dv12.getZ() };
    double dU0du22_coeffs[3] = { dU0du22.getX(), dU0du22.getY(), dU0du22.getZ() };
    double dU0dv22_coeffs[3] = { dU0dv22.getX(), dU0dv22.getY(), dU0dv22.getZ() };

    double dS0du11_coeffs[3] = { dS0du11.getX(), dS0du11.getY(), dS0du11.getZ() };
    double dS0dv11_coeffs[3] = { dS0dv11.getX(), dS0dv11.getY(), dS0dv11.getZ() };
    double dS0du21_coeffs[3] = { dS0du21.getX(), dS0du21.getY(), dS0du21.getZ() };
    double dS0dv21_coeffs[3] = { dS0dv21.getX(), dS0dv21.getY(), dS0dv21.getZ() };
    double dS0du12_coeffs[3] = { dS0du12.getX(), dS0du12.getY(), dS0du12.getZ() };
    double dS0dv12_coeffs[3] = { dS0dv12.getX(), dS0dv12.getY(), dS0dv12.getZ() };
    double dS0du22_coeffs[3] = { dS0du22.getX(), dS0du22.getY(), dS0du22.getZ() };
    double dS0dv22_coeffs[3] = { dS0dv22.getX(), dS0dv22.getY(), dS0dv22.getZ() };

    Matrix Au11 = Matrix ( 3,1,dU0du11_coeffs,true );
    Matrix Av11 = Matrix ( 3,1,dU0dv11_coeffs,true );
    Matrix Au21 = Matrix ( 3,1,dU0du21_coeffs,true );
    Matrix Av21 = Matrix ( 3,1,dU0dv21_coeffs,true );
    Matrix Au12 = Matrix ( 3,1,dU0du12_coeffs,true );
    Matrix Av12 = Matrix ( 3,1,dU0dv12_coeffs,true );
    Matrix Au22 = Matrix ( 3,1,dU0du22_coeffs,true );
    Matrix Av22 = Matrix ( 3,1,dU0dv22_coeffs,true );

    Matrix SAu11 = Matrix ( 3,1,dS0du11_coeffs,true );
    Matrix SAv11 = Matrix ( 3,1,dS0dv11_coeffs,true );
    Matrix SAu21 = Matrix ( 3,1,dS0du21_coeffs,true );
    Matrix SAv21 = Matrix ( 3,1,dS0dv21_coeffs,true );
    Matrix SAu12 = Matrix ( 3,1,dS0du12_coeffs,true );
    Matrix SAv12 = Matrix ( 3,1,dS0dv12_coeffs,true );
    Matrix SAu22 = Matrix ( 3,1,dS0du22_coeffs,true );
    Matrix SAv22 = Matrix ( 3,1,dS0dv22_coeffs,true );

    Matrix Bu11 = Ci11 * Au11;
    Matrix Bv11 = Ci11 * Av11;
    Matrix Bu21 = Ci21 * Au21;
    Matrix Bv21 = Ci21 * Av21;
    Matrix Bu12 = Ci12 * Au12;
    Matrix Bv12 = Ci12 * Av12;
    Matrix Bu22 = Ci22 * Au22;
    Matrix Bv22 = Ci22 * Av22;

    Matrix SBu11 = Ci11 * SAu11;
    Matrix SBv11 = Ci11 * SAv11;
    Matrix SBu21 = Ci21 * SAu21;
    Matrix SBv21 = Ci21 * SAv21;
    Matrix SBu12 = Ci12 * SAu12;
    Matrix SBv12 = Ci12 * SAv12;
    Matrix SBu22 = Ci22 * SAu22;
    Matrix SBv22 = Ci22 * SAv22;

    double b11[4] = { -1*Bu11 ( 0,0 ), -1*Bu11 ( 1,0 ), -1*Bv11 ( 0,0 ), -1*Bv11 ( 1,0 ) };
    double b21[4] = { -1*Bu21 ( 0,0 ), -1*Bu21 ( 1,0 ), -1*Bv21 ( 0,0 ), -1*Bv21 ( 1,0 ) };
    double b12[4] = { -1*Bu12 ( 0,0 ), -1*Bu12 ( 1,0 ), -1*Bv12 ( 0,0 ), -1*Bv12 ( 1,0 ) };
    double b22[4] = { -1*Bu22 ( 0,0 ), -1*Bu22 ( 1,0 ), -1*Bv22 ( 0,0 ), -1*Bv22 ( 1,0 ) };

    double Sb11[4] = { -1*SBu11 ( 0,0 ), -1*SBu11 ( 1,0 ), -1*SBv11 ( 0,0 ), -1*SBv11 ( 1,0 ) };
    double Sb21[4] = { -1*SBu21 ( 0,0 ), -1*SBu21 ( 1,0 ), -1*SBv21 ( 0,0 ), -1*SBv21 ( 1,0 ) };
    double Sb12[4] = { -1*SBu12 ( 0,0 ), -1*SBu12 ( 1,0 ), -1*SBv12 ( 0,0 ), -1*SBv12 ( 1,0 ) };
    double Sb22[4] = { -1*SBu22 ( 0,0 ), -1*SBu22 ( 1,0 ), -1*SBv22 ( 0,0 ), -1*SBv22 ( 1,0 ) };

    Matrix Srad11 = Matrix ( 2,2,b11,true );
    Matrix Srad21 = Matrix ( 2,2,b21,true );
    Matrix Srad12 = Matrix ( 2,2,b12,true );
    Matrix Srad22 = Matrix ( 2,2,b22,true );

    Matrix SSrad11 = Matrix ( 2,2,Sb11,true );
    Matrix SSrad21 = Matrix ( 2,2,Sb21,true );
    Matrix SSrad12 = Matrix ( 2,2,Sb12,true );
    Matrix SSrad22 = Matrix ( 2,2,Sb22,true );

//	std::cout << "SSrad: " << std::endl;
//	SSrad12.print();

    Matrix rSrad11 = Srad11 * quadAtom11->getR0();
    Matrix rSrad21 = Srad21 * quadAtom21->getR0();
    Matrix rSrad12 = Srad12 * quadAtom12->getR0();
    Matrix rSrad22 = Srad22 * quadAtom22->getR0();

//	std::cout << "rSrad: " << std::endl;
//	rSrad12.print();

    double Idat[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    Matrix Imat = Matrix ( 3,3,Idat,true );

    //11

    double pdat11[6] = {du11t.getX(), dv11t.getX(), du11t.getY(), dv11t.getY(), du11t.getZ(), dv11t.getZ() };
//	double pdat11[6] = {du11t.getX(), du11tp.getX(), du11t.getY(), du11tp.getY(), du11t.getZ(), du11tp.getZ() };
    Matrix pmat11 = Matrix ( 2,3,pdat11,true );

    double Udat11[3] = {U0_11.getX(), U0_11.getY(), U0_11.getZ() };
    Matrix Umat11 = Matrix ( 1,3,Udat11,true );

    Matrix Q11 = pmat11 * ( ( Umat11.t() * Umat11 ) - Imat );

    Matrix QQ11 = Q11 * Q11.t();
    Matrix QQi11;
    QQ11.inverse ( QQi11 );

    double rdat11[2] = {ru11, rv11};
    Matrix rmat11 = Matrix ( 2,1,rdat11,true );

    Matrix r2mat11 = pmat11 * Umat11.t();

    double dSdu11dat[6] = { dS0du11.getX(), dS0dv11.getX(), dS0du11.getY(), dS0dv11.getY(), dS0du11.getZ(), dS0dv11.getZ() };
    Matrix dSdu11mat = Matrix ( 2,3,dSdu11dat,true );

    double dUdu11dat[6] = { dU0du11.getX(), dU0dv11.getX(), dU0du11.getY(), dU0dv11.getY(), dU0du11.getZ(), dU0dv11.getZ() };
    Matrix dUdu11mat = Matrix ( 2,3,dUdu11dat,true );

    Matrix rS11 = ( dSdu11mat - ( rmat11*Umat11 ) ) * ( Q11.t() *QQi11 );

    Matrix testmult = Umat11 * Umat11.t();
//	testmult.print();

    //21

    double pdat21[6] = {du21t.getX(), dv21t.getX(), du21t.getY(), dv21t.getY(), du21t.getZ(), dv21t.getZ() };
//	double pdat21[6] = {du21t.getX(), du21tp.getX(), du21t.getY(), du21tp.getY(), du21t.getZ(), du21tp.getZ() };
    Matrix pmat21 = Matrix ( 2,3,pdat21,true );

    double Udat21[3] = {U0_21.getX(), U0_21.getY(), U0_21.getZ() };
    Matrix Umat21 = Matrix ( 1,3,Udat21,true );

    Matrix Q21 = pmat21 * ( ( Umat21.t() * Umat21 ) - Imat );

    Matrix QQ21 = Q21 * Q21.t();
    Matrix QQi21;
    QQ21.inverse ( QQi21 );

    double rdat21[2] = {ru21, rv21};
    Matrix rmat21 = Matrix ( 2,1,rdat21,true );

    Matrix r2mat21 = -1 * pmat21 * Umat21.t();

    double dSdu21dat[6] = { dS0du21.getX(), dS0dv21.getX(), dS0du21.getY(), dS0dv21.getY(), dS0du21.getZ(), dS0dv21.getZ() };
    Matrix dSdu21mat = Matrix ( 2,3,dSdu21dat,true );

    Matrix rS21 = ( dSdu21mat - ( rmat21*Umat21 ) ) * ( Q21.t() *QQi21 );

    //12

    double pdat12[6] = {du12t.getX(), dv12t.getX(), du12t.getY(), dv12t.getY(), du12t.getZ(), dv12t.getZ() };
//	double pdat12[6] = {du12t.getX(), du12tp.getX(), du12t.getY(), du12tp.getY(), du12t.getZ(), du12tp.getZ() };
    Matrix pmat12 = Matrix ( 2,3,pdat12,true );
    double Udat12[3] = {U0_12.getX(), U0_12.getY(), U0_12.getZ() };
    Matrix Umat12 = Matrix ( 1,3,Udat12,true );

    Matrix Q12 = pmat12 * ( ( Umat12.t() * Umat12 ) - Imat );

    Matrix QQ12 = Q12 * Q12.t();
    Matrix QQi12;
    QQ12.inverse ( QQi12 );

    double rdat12[2] = {ru12, rv12};
    Matrix rmat12 = Matrix ( 2,1,rdat12,true );

    Matrix r2mat12 = -1 * pmat12 * Umat12.t();

    double dSdu12dat[6] = { dS0du12.getX(), dS0dv12.getX(), dS0du12.getY(), dS0dv12.getY(), dS0du12.getZ(), dS0dv12.getZ() };
    Matrix dSdu12mat = Matrix ( 2,3,dSdu12dat,true );

    Matrix rS12 = ( dSdu12mat - ( rmat12*Umat12 ) ) * ( Q12.t() *QQi12 );

    //22

    double pdat22[6] = {du22t.getX(), dv22t.getX(), du22t.getY(), dv22t.getY(), du22t.getZ(), dv22t.getZ() };
//	double pdat22[6] = {du22t.getX(), du22tp.getX(), du22t.getY(), du22tp.getY(), du22t.getZ(), du22tp.getZ() };
    Matrix pmat22 = Matrix ( 2,3,pdat22,true );

    double Udat22[3] = {U0_22.getX(), U0_22.getY(), U0_22.getZ() };
    Matrix Umat22 = Matrix ( 1,3,Udat22,true );

    Matrix Q22 = pmat22 * ( ( Umat22.t() * Umat22 ) - Imat );

    Matrix QQ22 = Q22 * Q22.t();
    Matrix QQi22;
    QQ22.inverse ( QQi22 );

    double rdat22[2] = {ru22, rv22};
    Matrix rmat22 = Matrix ( 2,1,rdat22,true );

    Matrix r2mat22 = -1 * pmat22 * Umat22.t();

    double dSdu22dat[6] = { dS0du22.getX(), dS0dv22.getX(), dS0du22.getY(), dS0dv22.getY(), dS0du22.getZ(), dS0dv22.getZ() };
    Matrix dSdu22mat = Matrix ( 2,3,dSdu22dat,true );

    Matrix rS22 = ( dSdu22mat - ( rmat22*Umat22 ) ) * ( Q22.t() *QQi22 );

    //rSrad11 = SSrad11;
    //rSrad21 = SSrad21;
    //rSrad12 = SSrad12;
    //rSrad22 = SSrad22;

    //std::cout << "srads" << std::endl;
    //rSrad11.print();
    //rSrad21.print();
    //rSrad12.print();
    //rSrad22.print();

    Matrix rS12_2 = ( dSdu12mat + ( r2mat12*Umat12 ) ) * ( Q12.t() *QQi12 );

//	std::cout << "rSrad from formula, with pu/pv: " << std::endl;
//	rS12_2.t().print();
//
//	std::cout << "rSrad from formula, with ru/rv: " << std::endl;
//	rS12.t().print();

//	rmat11.print();
//	r2mat11.print();

    rSrad11 = rS11.t();
    rSrad12 = rS12.t();
    rSrad21 = rS21.t();
    rSrad22 = rS22.t();


    Matrix rSrad11RT = rSrad11 * rSrad11.t();
    Matrix rSrad21RT = rSrad21 * rSrad21.t();
    Matrix rSrad12RT = rSrad12 * rSrad12.t();
    Matrix rSrad22RT = rSrad22 * rSrad22.t();

    Matrix rSrad11TR = rSrad11.t() * rSrad11;
    Matrix rSrad21TR = rSrad21.t() * rSrad21;
    Matrix rSrad12TR = rSrad12.t() * rSrad12;
    Matrix rSrad22TR = rSrad22.t() * rSrad22;

//	cout << endl << "rSrad: " << endl;
//	cout << rSrad11(0,0) << ", " << rSrad11(0,1) << ", " << rSrad11(1,0) << ", " << rSrad11(1,1) << endl;
//	cout << rSrad21(0,0) << ", " << rSrad21(0,1) << ", " << rSrad21(1,0) << ", " << rSrad21(1,1) << endl;
//	cout << rSrad12(0,0) << ", " << rSrad12(0,1) << ", " << rSrad12(1,0) << ", " << rSrad12(1,1) << endl;
//	cout << rSrad22(0,0) << ", " << rSrad22(0,1) << ", " << rSrad22(1,0) << ", " << rSrad22(1,1) << endl;

    Vector L11, L21, L12, L22;
    Matrix V11, V21, V12, V22;

    Vector L1_11, L2_11, L1_21, L2_21, L1_12, L2_12, L1_22, L2_22;
    Matrix V1_11, V2_11, V1_21, V2_21, V1_12, V2_12, V1_22, V2_22;

    rSrad11.factorEV ( L11, V11, NON_SYM );
    rSrad21.factorEV ( L21, V21, NON_SYM );
    rSrad12.factorEV ( L12, V12, NON_SYM );
    rSrad22.factorEV ( L22, V22, NON_SYM );

    rSrad11RT.factorEV ( L1_11, V1_11, GENERAL );
    rSrad21RT.factorEV ( L1_21, V1_21, GENERAL );
    rSrad12RT.factorEV ( L1_12, V1_12, GENERAL );
    rSrad22RT.factorEV ( L1_22, V1_22, GENERAL );

    rSrad11TR.factorEV ( L2_11, V2_11, GENERAL );
    rSrad21TR.factorEV ( L2_21, V2_21, GENERAL );
    rSrad12TR.factorEV ( L2_12, V2_12, GENERAL );
    rSrad22TR.factorEV ( L2_22, V2_22, GENERAL );

//	std::cout << "eigvals" << std::endl;
//
//	L1_11.print();
//	L2_11.print();
//
//	L1_12.print();
//	L2_12.print();
//
//	L1_21.print();
//	L2_21.print();
//
//	L1_22.print();
//	L2_22.print();

    // Process corners to find out true eigenvalues

    double d1, d2, nd1, nd2, res1, res2, res3, res4, minres;
    Matrix dPP, dPN, dNP, dNN, test1, test2, test3, test4, temp1, temp2, temp3, temp4, EVs_11, EVs_21, EVs_12, EVs_22;

    //Corner 11

    d1 = sqrt ( L2_11 ( 0 ) );
    d2 = sqrt ( L2_11 ( 1 ) );

    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_11 * dPP * V2_11;
    test2 = V1_11 * dPN * V2_11;
    test3 = V1_11 * dNP * V2_11;
    test4 = V1_11 * dNN * V2_11;

    temp1 = test1 - rSrad11;
    temp2 = test2 - rSrad11;
    temp3 = test3 - rSrad11;
    temp4 = test4 - rSrad11;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_11 = dPP;

    if ( res2 < minres )
    {
        minres = res2;
        EVs_11 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_11 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_11 = dNN;
    }

    //std::cout << minres << std::endl;
    //Corner 21

    d1 = sqrt ( L2_21 ( 0 ) );
    d2 = sqrt ( L2_21 ( 1 ) );
    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_21 * dPP * V2_21;
    test2 = V1_21 * dPN * V2_21;
    test3 = V1_21 * dNP * V2_21;
    test4 = V1_21 * dNN * V2_21;

    temp1 = test1 - rSrad21;
    temp2 = test2 - rSrad21;
    temp3 = test3 - rSrad21;
    temp4 = test4 - rSrad21;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_21 = dPP;

    if ( res2 < minres )
    {
        minres = res2;
        EVs_21 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_21 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_21 = dNN;
    }

    //Corner 12

    d1 = sqrt ( L2_12 ( 0 ) );
    d2 = sqrt ( L2_12 ( 1 ) );
    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_12 * dPP * V2_12;
    test2 = V1_12 * dPN * V2_12;
    test3 = V1_12 * dNP * V2_12;
    test4 = V1_12 * dNN * V2_12;

    temp1 = test1 - rSrad12;
    temp2 = test2 - rSrad12;
    temp3 = test3 - rSrad12;
    temp4 = test4 - rSrad12;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_12 = dPP;

    if ( res2 < minres )
    {
        minres = res2;
        EVs_12 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_12 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_12 = dNN;
    }

    //std::cout << minres << std::endl;

    //Corner 22

    d1 = sqrt ( L2_22 ( 0 ) );
    d2 = sqrt ( L2_22 ( 1 ) );
    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_22 * dPP * V2_22;
    test2 = V1_22 * dPN * V2_22;
    test3 = V1_22 * dNP * V2_22;
    test4 = V1_22 * dNN * V2_22;

    temp1 = test1 - rSrad22;
    temp2 = test2 - rSrad22;
    temp3 = test3 - rSrad22;
    temp4 = test4 - rSrad22;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_22 = dPP;



    if ( res2 < minres )
    {
        minres = res2;
        EVs_22 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_22 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_22 = dNN;
    }

//	EVs_11.print();
//	EVs_21.print();
//	EVs_12.print();
//	EVs_22.print();

    //std::cout << minres << std::endl;

//	std::cout << "****BEGIN QUAD****" << std::endl;
//
//	std::cout << "Corner 0,0: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_11.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_11.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_11.print();
//
//	std::cout << "Corner 1,0: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_21.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_21.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_21.print();
//
//	std::cout << "Corner 0,1: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_12.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_12.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_12.print();
//
//	std::cout << "Corner 1,1: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_22.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_22.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_22.print();
//	std::cout << "rSrad: " << std::endl;
//	rSrad22.print();
//
//	std::cout << "**** END QUAD ****" << std::endl;


    double Lam1_11, Lam2_11, Lam1_21, Lam2_21, Lam1_12, Lam2_12, Lam1_22, Lam2_22;
    Vector e1_11, e2_11, e1_21, e2_21, e1_12, e2_12, e1_22, e2_22;

    Lam1_11 = EVs_11 ( 0,0 );
    Lam2_11 = EVs_11 ( 1,1 );
    Lam1_21 = EVs_21 ( 0,0 );
    Lam2_21 = EVs_21 ( 1,1 );
    Lam1_12 = EVs_12 ( 0,0 );
    Lam2_12 = EVs_12 ( 1,1 );
    Lam1_22 = EVs_22 ( 0,0 );
    Lam2_22 = EVs_22 ( 1,1 );

    e1_11 = V1_11;
    e2_11 = V2_11;
    e1_21 = V1_21;
    e2_21 = V2_21;
    e1_12 = V1_12;
    e2_12 = V2_12;
    e1_22 = V1_22;
    e2_22 = V2_22;

    // Corner 1
//	if (L11(0) > L11(1))
//	{
//		Lam1_11 = L11(0);
//		Lam2_11 = L11(1);
//
//		e1_11 = V11.getColumn(0);
//		e2_11 = V11.getColumn(1);
//	}
//	else
//	{
//		Lam1_11 = L11(1);
//		Lam2_11 = L11(0);
//
//		e1_11 = V11.getColumn(1);
//		e2_11 = V11.getColumn(0);
//	}
//
//	// Corner 2
//	if (L21(0) > L21(1))
//	{
//		Lam1_21 = L21(0);
//		Lam2_21 = L21(1);
//
//		e1_21 = V21.getColumn(0);
//		e2_21 = V21.getColumn(1);
//	}
//	else
//	{
//		Lam1_21 = L21(1);
//		Lam2_21 = L21(0);
//
//		e1_21 = V21.getColumn(1);
//		e2_21 = V21.getColumn(0);
//	}
//
//	//Corner 3
//	if (L12(0) > L12(1))
//	{
//		Lam1_12 = L12(0);
//		Lam2_12 = L12(1);
//
//		e1_12 = V12.getColumn(0);
//		e2_12 = V12.getColumn(1);
//	}
//	else
//	{
//		Lam1_12 = L12(1);
//		Lam2_12 = L12(0);
//
//		e1_12 = V12.getColumn(1);
//		e2_12 = V12.getColumn(0);
//	}
//
//	//Corner 4
//	if (L22(0) > L22(1))
//	{
//		Lam1_22 = L22(0);
//		Lam2_22 = L22(1);
//
//		e1_22 = V22.getColumn(0);
//		e2_22 = V22.getColumn(1);
//	}
//	else
//	{
//		Lam1_22 = L22(1);
//		Lam2_22 = L22(0);
//
//		e1_22 = V22.getColumn(1);
//		e2_22 = V22.getColumn(0);
//	}

//	cout << endl << "Eigenvectors: " << endl;
//	cout << "(" << e1_11(0) << ", " << e1_11(1) << "), (" << e2_11(0) << ", " << e2_11(1) << ")" << endl;
//	cout << "(" << e1_21(0) << ", " << e1_21(1) << "), (" << e2_21(0) << ", " << e2_21(1) << ")" << endl;
//	cout << "(" << e1_12(0) << ", " << e1_12(1) << "), (" << e2_12(0) << ", " << e2_12(1) << ")" << endl;
//	cout << "(" << e1_22(0) << ", " << e1_22(1) << "), (" << e2_22(0) << ", " << e2_22(1) << ")" << endl;

    double testlambda[4] = { Lam1_11, 0, 0, Lam2_11 };

    Matrix testL = Matrix ( 2,2,testlambda,true );
    Matrix testinv;

    Matrix testsrad = V11 * testL * V11.inverse ( testinv );

//	cout << endl << endl;
//	cout << testsrad(0,0) << ", " << testsrad(0,1) << ", " << testsrad(1,0) << ", " << testsrad(1,1) << endl;
//	cout << endl << endl;

//	cout << endl << "Eigenvalues: " << endl;
//	cout << Lam1_11 << ", " << Lam2_11 << endl;
//	cout << Lam1_21 << ", " << Lam2_21 << endl;
//	cout << Lam1_12 << ", " << Lam2_12 << endl;
//	cout << Lam1_22 << ", " << Lam2_22 << endl;

    int num_iter = 10;
    double curru = 0;
    double currv = 0;
    double du = u / num_iter;
    double dv = v / num_iter;

    double logLam1_11 = log ( 1 - Lam1_11 );
    double logLam2_11 = log ( 1 - Lam2_11 );
    double logLam1_21 = log ( 1 - Lam1_21 );
    double logLam2_21 = log ( 1 - Lam2_21 );
    double logLam1_12 = log ( 1 - Lam1_12 );
    double logLam2_12 = log ( 1 - Lam2_12 );
    double logLam1_22 = log ( 1 - Lam1_22 );
    double logLam2_22 = log ( 1 - Lam2_22 );

    double logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    double logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    double avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    double avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    double Lam1 = 1 - exp ( logAvg1 );
    double Lam2 = 1 - exp ( logAvg2 );

    //Lam1 = logAvg1;
    //Lam2 = logAvg2;

    // Uses a form of bilinear interpolation to get the eigenvectors: This should probably be a Frechet mean of the thetas
    // on the unit circle in the future

    /*double theta1_11 = atan2(e1_11(1), e1_11(0));
    double theta2_11 = atan2(e2_11(1), e2_11(0));
    double theta1_21 = atan2(e1_21(1), e1_21(0));
    double theta2_21 = atan2(e2_21(1), e2_21(0));
    double theta1_12 = atan2(e1_12(1), e1_12(0));
    double theta2_12 = atan2(e2_12(1), e2_12(0));
    double theta1_22 = atan2(e1_22(1), e1_22(0));
    double theta2_22 = atan2(e2_22(1), e2_22(0));*/

    double PI = 3.14159265;

    double avgx1, avgy1, avgx2, avgy2, avgtheta1, avgtheta2;//, at1_11_21, at2_11_21, at1_12_22, at2_12_22;

//	std::cout << "eigs" << std::endl;
//	e1_11.print();
//	e1_21.print();
//	e1_12.print();
//	e1_22.print();

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    // Go from 11 to 21 (u direction)
    //double thetadist1_11_21 = theta1_21 - theta1_11;
    //if (thetadist1_11_21 > PI)
    //{
    //	theta1_21 -= 2*PI;
    //	thetadist1_11_21 = theta1_21 - theta1_11;
    //}
    //else if (thetadist1_11_21 < -PI)
    //{
    //	theta1_21 += 2*PI;
    //	thetadist1_11_21 = theta1_21 - theta1_11;
    //}
    //at1_11_21 = theta1_11 + u*thetadist1_11_21;

    //double thetadist2_11_21 = theta2_21 - theta2_11;
    //if (thetadist2_11_21 > PI)
    //{
    //	theta2_21 -= 2*PI;
    //	thetadist2_11_21 = theta2_21 - theta2_11;
    //}
    //else if (thetadist2_11_21 < -PI)
    //{
    //	theta2_21 += 2*PI;
    //	thetadist2_11_21 = theta2_21 - theta2_11;
    //}
    //at2_11_21 = theta2_11 + u*thetadist2_11_21;

    //// Go from 12 to 22 (u direction)
    //double thetadist1_12_22 = theta1_22 - theta1_12;
    //if (thetadist1_12_22 > PI)
    //{
    //	theta1_22 -= 2*PI;
    //	thetadist1_12_22 = theta1_22 - theta1_12;
    //}
    //else if (thetadist1_12_22 < -PI)
    //{
    //	theta1_22 += 2*PI;
    //	thetadist1_12_22 = theta1_22 - theta1_12;
    //}
    //at1_12_22 = theta1_12 + u*thetadist1_12_22;

    //double thetadist2_12_22 = theta2_22 - theta2_12;
    //if (thetadist2_12_22 > PI)
    //{
    //	theta2_22 -= 2*PI;
    //	thetadist2_12_22 = theta2_22 - theta2_12;
    //}
    //else if (thetadist2_12_22 < -PI)
    //{
    //	theta2_22 += 2*PI;
    //	thetadist2_12_22 = theta2_22 - theta2_12;
    //}
    //at2_12_22 = theta2_12 + u*thetadist2_12_22;

    //// Now go v direction
    //double thetadist1_v = at1_12_22 - at1_11_21;
    //if (thetadist1_v > PI)
    //{
    //	at1_12_22 -= 2*PI;
    //	thetadist1_v = at1_12_22 - at1_11_21;
    //}
    //else if (thetadist1_v < -PI)
    //{
    //	at1_12_22 += 2*PI;
    //	thetadist1_v = at1_12_22 - at1_11_21;
    //}
    //avgtheta1 = at1_11_21 + v*thetadist1_v;

    //double thetadist2_v = at2_12_22 - at2_11_21;
    //if (thetadist2_v > PI)
    //{
    //	at2_12_22 -= 2*PI;
    //	thetadist2_v = at2_12_22 - at2_11_21;
    //}
    //else if (thetadist2_v < -PI)
    //{
    //	at2_12_22 += 2*PI;
    //	thetadist2_v = at2_12_22 - at2_11_21;
    //}
    //avgtheta2 = at2_11_21 + v*thetadist2_v;

    double neweigenv[4] = { cos ( avgtheta1 ), sin ( avgtheta1 ), cos ( avgtheta2 ), sin ( avgtheta2 ) };
    Matrix newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    Matrix newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    Matrix NewV = Matrix ( 2,2,neweigenv,true );

    double newlambda[4] = { Lam1, 0, 0, Lam2 };

    Matrix NewL = Matrix ( 2,2,newlambda,true );

//	NewL.print();
    /*cout << V11(0,0) << ", " << V11(1,0)  << ", " << V11(0,1) << ", " << V11(1,1) << endl;
    cout << V21(0,0) << ", " << V21(1,0)  << ", " << V21(0,1) << ", " << V21(1,1) << endl;
    cout << V12(0,0) << ", " << V12(1,0)  << ", " << V12(0,1) << ", " << V12(1,1) << endl;
    cout << V22(0,0) << ", " << V22(1,0)  << ", " << V22(0,1) << ", " << V22(1,1) << endl;
    cout << NewV(0,0) << ", " << NewV(1,0) << ", " << NewV(0,1) << ", " << NewV(1,1) << endl;*/

    Matrix NewVi;
    NewV.inverse ( NewVi );

    //Matrix NewrSrad = NewV * NewL * NewVi;

    Matrix NewrSrad = newLeft * NewL * newRight;
//	NewrSrad.print();

    if ( side == 0 )
    {
        U0_11 = quadAtom11->getU0();
    }
    else
    {
        U0_11 = quadAtom11->getU1();
    }

    //Vector3D du11t_norm = du11t / du11t.norm();
    //Vector3D dv11t_norm = dv11t / dv11t.norm();



    double Pdata[6] = { du11t.getX(), dv11t.getX(), du11t.getY(), dv11t.getY(), du11t.getZ(), dv11t.getZ() };
    Matrix P = Matrix ( 2,3,Pdata,true );

//	du11.print();
//	dv11.print();
//	P.print();

    double Udata[3] = { U0_11.getX(), U0_11.getY(), U0_11.getZ() };
    Matrix U = Matrix ( 1,3,Udata,true );

    //double ru11 = du11 * U0_11;
    //double rv11 = dv11 * U0_11;

    //double ru21 = du21 * U0_21;
    //double rv21 = dv21 * U0_21;

    //double ru12 = du12 * U0_12;
    //double rv12 = dv12 * U0_12;

    //double ru22 = du22 * U0_22;
    //double rv22 = dv22 * U0_22;

//	cout << endl << "Ru & Rv: " << endl;
//	cout << ru11 << ", " << rv11 << endl;
//	cout << ru21 << ", " << rv21 << endl;
//	cout << ru12 << ", " << rv12 << endl;
//	cout << ru22 << ", " << rv22 << endl;

    double dSdu11data[6] = { dS0du11.getX(), dS0dv11.getX(), dS0du11.getY(), dS0dv11.getY(), dS0du11.getZ(), dS0dv11.getZ() };
    Matrix dSdu11 = Matrix ( 2,3,dSdu11data,true );

    double rdata[2] = { ru11, rv11 };
    Matrix rmat = Matrix ( 2,1,rdata,true );

    double Idata[9] = {1,0,0,0,1,0,0,0,1};
    Matrix I = Matrix ( 3,3,Idata,true );

    Matrix Q = P * ( ( U.t() * U ) - I );

    Matrix QQi;

    Matrix QQ = Q * Q.t();
    QQ.inverse ( QQi );

    Matrix UtU, UtUp, UtUi;

    UtU = U.t() * U;
    Q = P * ( UtU - I );
    UtUp = I + UtU;
    UtUp.inverse(UtUi);

    Matrix UtUpu = I + du*UtU;
    Matrix UtUpv = I + dv*UtU;

    Matrix UtUiu, UtUiv;
    UtUpu.inverse(UtUiu);
    UtUpv.inverse(UtUiv);


    Matrix R = -1 * P * U.t();

    Matrix dS = ( rSrad11.t() * Q ) + ( rmat * U );

    Matrix testrSrad11 = ( dSdu11 - ( rmat*U ) ) * ( Q.t() *QQi );
//	std::cout << "rs11: " << std::endl;
//	rSrad11.print();
//	std::cout << "trs11: " << std::endl;
//	testrSrad11.t().print();

    //Q.print();
    //Q11.print();

    Vector dSdu = dS.getRow ( 0 );
    Vector dSdv = dS.getRow ( 1 );

    //Vector3D dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
    //Vector3D dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );

    Vector Suu = ((rSrad11.t().getColumn(0) * Q.getRow(0)) -
                                    2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
    Vector Suv = ((rSrad11.t().getColumn(1) * Q.getRow(1)) -
                            2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

//		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
//						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
//		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
//						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

    Vector3D dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
    Vector3D dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

    Vector3D S0;

    if ( side == 0 )
    {
        S0 = U0_11 * quadAtom11->getR0();
    }
    else
    {
        S0 = U0_11 * quadAtom11->getR1();
    }


    Vector3D NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru = curru + du;
    currv = currv + dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    Vector3D newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );


    Matrix temprSrad = 0.25* ( rSrad11 + rSrad12 + rSrad21 + rSrad22 );

    //U0_11.print();
    for ( int i = 1; i < num_iter; i++ )
    {

//		std::cout << i << std::endl;
        // To find rSrad, need to calculate lambdas and eigenvectors
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        Q = P * ( ( U.t() * U ) - I );

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();

        UtU = U.t() * U;
        Q = P * ( UtU - I );
        UtUp = I + UtU;
        UtUp.inverse(UtUi);

        Matrix UtUpu = I + UtU;
        Matrix UtUpv = I + UtU;

        Matrix UtUiu, UtUiv;
        UtUpu.inverse(UtUiu);
        UtUpv.inverse(UtUiv);

        //Matrix Q2 = NewrSrad.t() * P



//        dS = ( NewrSrad.t() * Q ) + ( R * U );
//
//        dSdu = dS.getRow ( 0 );
//        dSdv = dS.getRow ( 1 );
//
//        dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
//        dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );
//
//        Vector3D TempSpoke = r * NewSpoke;
//
//        NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );
        //NewSpoke.print();

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
        Matrix S = Matrix(1,3,Spokedata,true);


        Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
                                2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
        Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
                                2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

//		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
//						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
//		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
//						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

        dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
        dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

        Vector3D TempSpoke = r*NewSpoke;
        NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);

        curru = curru + du;
        currv = currv + dv;

    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );

    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );
    Vector3D newbound = Vector3D ( Bxn ( 0,0 ), Byn ( 0,0 ), Bzn ( 0,0 ) );

    //cout << "{" << quadAtom11->getX().getX() + (quadAtom11->getR0() * quadAtom11->getU0().getX()) << ", " << quadAtom11->getX().getY() + (quadAtom11->getR0() * quadAtom11->getU0().getY()) << ", " << quadAtom11->getX().getZ() + (quadAtom11->getR0() * quadAtom11->getU0().getZ()) << "}" << endl;
    //cout << "{" << quadAtom21->getX().getX() + (quadAtom21->getR0() * quadAtom21->getU0().getX()) << ", " << quadAtom21->getX().getY() + (quadAtom21->getR0() * quadAtom21->getU0().getY()) << ", " << quadAtom21->getX().getZ() + (quadAtom21->getR0() * quadAtom21->getU0().getZ()) << "}" << endl;
    //cout << "{" << quadAtom12->getX().getX() + (quadAtom12->getR0() * quadAtom12->getU0().getX()) << ", " << quadAtom12->getX().getY() + (quadAtom12->getR0() * quadAtom12->getU0().getY()) << ", " << quadAtom12->getX().getZ() + (quadAtom12->getR0() * quadAtom12->getU0().getZ()) << "}" << endl;
    //cout << "{" << quadAtom22->getX().getX() + (quadAtom22->getR0() * quadAtom22->getU0().getX()) << ", " << quadAtom22->getX().getY() + (quadAtom22->getR0() * quadAtom22->getU0().getY()) << ", " << quadAtom22->getX().getZ() + (quadAtom22->getR0() * quadAtom22->getU0().getZ()) << "}" << endl;
    //cout << "{" << newpos.getX() + NewSpoke.getX() << ", " << newpos.getY() + NewSpoke.getY() << ", " << newpos.getZ() + NewSpoke.getZ() << endl;

    //cout << quadAtom11->getR0() << endl;
    //cout << quadAtom21->getR0() << endl;
    //cout << quadAtom12->getR0() << endl;
    //cout << quadAtom22->getR0() << endl;
    //cout << NewSpoke.norm() << endl;

    //cout << quadAtom11->getX().getX() << ", " << quadAtom11->getX().getY() << ", " << quadAtom11->getX().getZ() << endl;
    //cout << quadAtom21->getX().getX() << ", " << quadAtom21->getX().getY() << ", " << quadAtom21->getX().getZ() << endl;
    //cout << quadAtom12->getX().getX() << ", " << quadAtom12->getX().getY() << ", " << quadAtom12->getX().getZ() << endl;
    //cout << quadAtom22->getX().getX() << ", " << quadAtom22->getX().getY() << ", " << quadAtom22->getX().getZ() << endl;

    double NewR = NewSpoke.normalize();
    Vector3D htestD = newbound - newpos;
    double htestR = htestD.normalize();


    //cout << newpos.getX() << ", " << newpos.getY() << ", " << newpos.getZ() << endl;

    //cout << newpos.getX() << ", " << NewR << ", " << NewSpoke.getX() << endl;
    M3DSpoke* hspoke = new M3DSpoke ( newpos, htestD, htestR ); // This is for testing purposes
    //M3DSpoke* spoke11 = new M3DSpoke(newpos, htestD, htestR); // This is for testing purposes
    M3DSpoke* spoke11 = new M3DSpoke ( newpos, NewSpoke, NewR );

    // from corner 21

    curru = 1;
    currv = 0;
    du = ( u-1 ) /num_iter;
    dv = v / num_iter;

    //std::cout << "du: " << du << ", dv: " << dv << std::endl;

    logLam1_11 = log ( 1 - Lam1_11 );
    logLam2_11 = log ( 1 - Lam2_11 );
    logLam1_21 = log ( 1 - Lam1_21 );
    logLam2_21 = log ( 1 - Lam2_21 );
    logLam1_12 = log ( 1 - Lam1_12 );
    logLam2_12 = log ( 1 - Lam2_12 );
    logLam1_22 = log ( 1 - Lam1_22 );
    logLam2_22 = log ( 1 - Lam2_22 );

    logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    Lam1 = 1 - exp ( logAvg1 );
    Lam2 = 1 - exp ( logAvg2 );

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    neweigenv[0] = cos ( avgtheta1 );
    neweigenv[1] = sin ( avgtheta1 );
    neweigenv[2] = cos ( avgtheta2 );
    neweigenv[3] = sin ( avgtheta2 );

    newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    NewV = Matrix ( 2,2,neweigenv,true );

    newlambda[0] = Lam1;
    newlambda[1] = 0;
    newlambda[2] = 0;
    newlambda[3] = Lam2;

    NewL = Matrix ( 2,2,newlambda,true );

    NewV.inverse ( NewVi );

    NewrSrad = newLeft * NewL * newRight;


    Pdata[0] = du21t.getX();
    Pdata[1] = dv21t.getX();
    Pdata[2] = du21t.getY();
    Pdata[3] = dv21t.getY();
    Pdata[4] = du21t.getZ();
    Pdata[5] = dv21t.getZ();
    P = Matrix ( 2,3,Pdata,true );

    Udata[0] = U0_21.getX();
    Udata[1] = U0_21.getY();
    Udata[2] = U0_21.getZ();
    U = Matrix ( 1,3,Udata,true );

    double dSdu21data[6] = { dS0du21.getX(), dS0dv21.getX(), dS0du21.getY(), dS0dv21.getY(), dS0du21.getZ(), dS0dv21.getZ() };
    Matrix dSdu21 = Matrix ( 2,3,dSdu21data,true );

    rdata[0] = ru21;
    rdata[1] = rv21;
    rmat = Matrix ( 2,1,rdata,true );

    Q = P * ( ( U.t() * U ) - I );

    QQi;

    QQ = Q * Q.t();
    QQ.inverse ( QQi );

    R = -1 * P * U.t();

    dS = ( rSrad21.t() * Q ) + ( rmat * U );

    dSdu = dS.getRow ( 0 );
    dSdv = dS.getRow ( 1 );

    dSu = Vector3D ( dSdu ( 0 ),dSdu ( 1 ),dSdu ( 2 ) );
    dSv = Vector3D ( dSdv ( 0 ),dSdv ( 1 ),dSdv ( 2 ) );

    S0 = U0_21 * r21;

    NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru += du;
    currv += dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    for ( int i = 1; i < num_iter; i++ )
    {
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        double pdatam[6] = { npu.getX() - r*NewSpoke.getX(), npv.getX() - r*NewSpoke.getX(),
                        npu.getY() - r*NewSpoke.getY(), npv.getY() - r*NewSpoke.getY(),
                        npu.getZ() - r*NewSpoke.getZ(), npv.getZ() - r*NewSpoke.getZ()
                };
        Matrix Pm = Matrix(2,3,pdatam,true);
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        UtU = U.t() * U;
        Q = P * ( UtU - I );
        UtUp = I + UtU;
        UtUp.inverse(UtUi);

        Matrix UtUpu = I + UtU;
        Matrix UtUpv = I + UtU;

        Matrix UtUiu, UtUiv;
        UtUpu.inverse(UtUiu);
        UtUpv.inverse(UtUiv);

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();



        dS = ( NewrSrad.t() * Q ) + ( R * U );

        dSdu = dS.getRow ( 0 );
        dSdv = dS.getRow ( 1 );



        //Vector3D TempSpoke = r * NewSpoke;

        //NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
        Matrix S = Matrix(1,3,Spokedata,true);


        Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
                                2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
        Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
                                2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

//		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
//						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
//		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
//						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

        dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
        dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

        Vector3D TempSpoke = r*NewSpoke;
        NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);


        //NewSpoke.print();

        curru = curru + du;
        currv = currv + dv;
    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    NewR = NewSpoke.normalize();
    M3DSpoke* spoke21 = new M3DSpoke ( newpos, NewSpoke, NewR );

    // from corner 12

    curru = 0;
    currv = 1;
    du = u /num_iter;
    dv = ( v-1 ) / num_iter;

    //std::cout << "du: " << du << ", dv: " << dv << std::endl;

    logLam1_11 = log ( 1 - Lam1_11 );
    logLam2_11 = log ( 1 - Lam2_11 );
    logLam1_21 = log ( 1 - Lam1_21 );
    logLam2_21 = log ( 1 - Lam2_21 );
    logLam1_12 = log ( 1 - Lam1_12 );
    logLam2_12 = log ( 1 - Lam2_12 );
    logLam1_22 = log ( 1 - Lam1_22 );
    logLam2_22 = log ( 1 - Lam2_22 );

    logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    Lam1 = 1 - exp ( logAvg1 );
    Lam2 = 1 - exp ( logAvg2 );

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    neweigenv[0] = cos ( avgtheta1 );
    neweigenv[1] = sin ( avgtheta1 );
    neweigenv[2] = cos ( avgtheta2 );
    neweigenv[3] = sin ( avgtheta2 );

    newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    NewV = Matrix ( 2,2,neweigenv,true );

    newlambda[0] = Lam1;
    newlambda[3] = Lam2;

    NewL = Matrix ( 2,2,newlambda,true );

    NewV.inverse ( NewVi );

    NewrSrad = newLeft * NewL * newRight;

    Pdata[0] = du12t.getX();
    Pdata[1] = dv12t.getX();
    Pdata[2] = du12t.getY();
    Pdata[3] = dv12t.getY();
    Pdata[4] = du12t.getZ();
    Pdata[5] = dv12t.getZ();
    P = Matrix ( 2,3,Pdata,true );

    Udata[0] = U0_12.getX();
    Udata[1] = U0_12.getY();
    Udata[2] = U0_12.getZ();
    U = Matrix ( 1,3,Udata,true );

    double dSdu12data[6] = { dS0du12.getX(), dS0dv12.getX(), dS0du12.getY(), dS0dv12.getY(), dS0du12.getZ(), dS0dv12.getZ() };
    Matrix dSdu12 = Matrix ( 2,3,dSdu12data,true );

    rdata[0] = ru12;
    rdata[1] = rv12;
    rmat = Matrix ( 2,1,rdata,true );

    Q = P * ( ( U.t() * U ) - I );

    QQi;

    QQ = Q * Q.t();
    QQ.inverse ( QQi );

    R = -1 * P * U.t();

    dS = ( rSrad12.t() * Q ) + ( rmat * U );

    dSdu = dS.getRow ( 0 );
    dSdv = dS.getRow ( 1 );

    dSu = Vector3D ( dSdu ( 0 ),dSdu ( 1 ),dSdu ( 2 ) );
    dSv = Vector3D ( dSdv ( 0 ),dSdv ( 1 ),dSdv ( 2 ) );

    S0 = U0_12 * r12;

    NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru += du;
    currv += dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    for ( int i = 1; i < num_iter; i++ )
    {
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        Q = P * ( ( U.t() * U ) - I );

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();

        UtU = U.t() * U;
                Q = P * ( UtU - I );
                UtUp = I + UtU;
                UtUp.inverse(UtUi);

                Matrix UtUpu = I + UtU;
                Matrix UtUpv = I + UtU;

                Matrix UtUiu, UtUiv;
                UtUpu.inverse(UtUiu);
                UtUpv.inverse(UtUiv);

        dS = ( NewrSrad.t() * Q ) + ( R * U );

        dSdu = dS.getRow ( 0 );
        dSdv = dS.getRow ( 1 );

        dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
        dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );

        //Vector3D TempSpoke = r * NewSpoke;

        //NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );
        //NewSpoke.print();

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
                Matrix S = Matrix(1,3,Spokedata,true);


                Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
                                        2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
                Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
                                        2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

        //		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
        //						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
        //		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
        //						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

                dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
                dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

                Vector3D TempSpoke = r*NewSpoke;
                NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);

        curru = curru + du;
        currv = currv + dv;
    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    NewR = NewSpoke.normalize();
    M3DSpoke* spoke12 = new M3DSpoke ( newpos, NewSpoke, NewR );

    // from corner 22

    curru = 1;
    currv = 1;
    du = ( u-1 ) /num_iter;
    dv = ( v-1 ) / num_iter;

    //std::cout << "du: " << du << ", dv: " << dv << std::endl;

    logLam1_11 = log ( 1 - Lam1_11 );
    logLam2_11 = log ( 1 - Lam2_11 );
    logLam1_21 = log ( 1 - Lam1_21 );
    logLam2_21 = log ( 1 - Lam2_21 );
    logLam1_12 = log ( 1 - Lam1_12 );
    logLam2_12 = log ( 1 - Lam2_12 );
    logLam1_22 = log ( 1 - Lam1_22 );
    logLam2_22 = log ( 1 - Lam2_22 );

    logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    Lam1 = 1 - exp ( logAvg1 );
    Lam2 = 1 - exp ( logAvg2 );

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    neweigenv[0] = cos ( avgtheta1 );
    neweigenv[1] = sin ( avgtheta1 );
    neweigenv[2] = cos ( avgtheta2 );
    neweigenv[3] = sin ( avgtheta2 );

    newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    NewV = Matrix ( 2,2,neweigenv,true );

    newlambda[0] = Lam1;
    newlambda[3] = Lam2;

    NewL = Matrix ( 2,2,newlambda,true );

    NewV.inverse ( NewVi );

    NewrSrad = newLeft * NewL * newRight;


    Pdata[0] = du22t.getX();
    Pdata[1] = dv22t.getX();
    Pdata[2] = du22t.getY();
    Pdata[3] = dv22t.getY();
    Pdata[4] = du22t.getZ();
    Pdata[5] = dv22t.getZ();
    P = Matrix ( 2,3,Pdata,true );

    Udata[0] = U0_22.getX();
    Udata[1] = U0_22.getY();
    Udata[2] = U0_22.getZ();
    U = Matrix ( 1,3,Udata,true );

    double dSdu22data[6] = { dS0du22.getX(), dS0dv22.getX(), dS0du22.getY(), dS0dv22.getY(), dS0du22.getZ(), dS0dv22.getZ() };
    Matrix dSdu22 = Matrix ( 2,3,dSdu22data,true );

    rdata[0] = ru22;
    rdata[1] = rv22;
    rmat = Matrix ( 2,1,rdata,true );

    Q = P * ( ( U.t() * U ) - I );

    QQi;

    QQ = Q * Q.t();
    QQ.inverse ( QQi );

    R = -1 * P * U.t();

    dS = ( rSrad22.t() * Q ) + ( rmat * U );

    dSdu = dS.getRow ( 0 );
    dSdv = dS.getRow ( 1 );

    dSu = Vector3D ( dSdu ( 0 ),dSdu ( 1 ),dSdu ( 2 ) );
    dSv = Vector3D ( dSdv ( 0 ),dSdv ( 1 ),dSdv ( 2 ) );

    S0 = U0_22 * r22;

    NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru += du;
    currv += dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    for ( int i = 1; i < num_iter; i++ )
    {
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        Q = P * ( ( U.t() * U ) - I );

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();

        UtU = U.t() * U;
                Q = P * ( UtU - I );
                UtUp = I + UtU;
                UtUp.inverse(UtUi);

                Matrix UtUpu = I + UtU;
                Matrix UtUpv = I + UtU;

                Matrix UtUiu, UtUiv;
                UtUpu.inverse(UtUiu);
                UtUpv.inverse(UtUiv);

        dS = ( NewrSrad.t() * Q ) + ( R * U );

        dSdu = dS.getRow ( 0 );
        dSdv = dS.getRow ( 1 );

        dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
        dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );

        //Vector3D TempSpoke = r * NewSpoke;

        //NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );
        //NewSpoke.print();

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
                Matrix S = Matrix(1,3,Spokedata,true);


                Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
                                        2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
                Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
                                        2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

        //		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
        //						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
        //		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
        //						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

                dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
                dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

                Vector3D TempSpoke = r*NewSpoke;
                NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);

        curru = curru + du;
        currv = currv + dv;
    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    NewR = NewSpoke.normalize();
    M3DSpoke* spoke22 = new M3DSpoke ( newpos, NewSpoke, NewR );

    // Mean of spokes:

    Vector3D s11 = spoke11->getS();
    Vector3D s21 = spoke21->getS();
    Vector3D s12 = spoke12->getS();
    Vector3D s22 = spoke22->getS();

    //spoke11->draw();
    //spoke12->draw();
    //spoke21->draw();
    //spoke22->draw();

    Vector3D meanspoke = ( ( 1-u ) * ( 1-v ) *s11 + ( u ) * ( 1-v ) *s21 + ( 1-u ) * ( v ) *s12 + ( u ) * ( v ) *s22 );
    NewR = meanspoke.normalize();

//    M3DSpoke* mean = new M3DSpoke ( newpos, meanspoke, NewR );
    M3DSpoke mean( newpos, meanspoke, NewR );
    //mean->draw();

    //std::cout << "*****ENDING   INTERPOLATION*****" << std::endl;
    //return hspoke;
//         std::cout << "one" << std::endl;

    delete spoke22;
    delete spoke12;
    delete spoke21;
    delete spoke11;
    delete hspoke;

    return mean;
    //return spoke11;

}



M3DSpoke M3DQuadInterpolater::interpolateSpoke_hermite ( M3DFigure *figure, double u, double v, int side)
{
    int ubase = ( int ) floor ( u );
    int vbase = ( int ) floor ( v );

    u = u - ubase;
    v = v - vbase;

    M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );

    if ( ubase == tempFigure->getRowCount() - 1 )
    {
        ubase = ubase - 1;
        u=1;
    }
    if ( vbase == tempFigure->getColumnCount() - 1 )
    {
        vbase = vbase - 1;
        v=1;
    }
    // Get four corner atoms

    M3DPrimitive *atom11 = tempFigure->getPrimitivePtr ( ubase,vbase );
    M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*> ( atom11 );

    M3DPrimitive *atom21 = tempFigure->getPrimitivePtr ( ubase+1,vbase );
    M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*> ( atom21 );

    M3DPrimitive *atom12 = tempFigure->getPrimitivePtr ( ubase,vbase+1 );
    M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*> ( atom12 );

    M3DPrimitive *atom22 = tempFigure->getPrimitivePtr ( ubase+1,vbase+1 );
    M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*> ( atom22 );

    Vector3D x11 = quadAtom11->getX();
    Vector3D x21 = quadAtom21->getX();
    Vector3D x12 = quadAtom12->getX();
    Vector3D x22 = quadAtom22->getX();


    double r11, r12, r21, r22, ru11, ru12, ru21, ru22, rv11, rv12, rv21, rv22;

    Vector3D U0_11, U0_12, U0_21, U0_22, B0_11, B0_12, B0_21, B0_22, S0_11, S0_12,S0_21,S0_22;

    if ( side == 0 )
    {
        U0_11 = quadAtom11->getU0();
        U0_21 = quadAtom21->getU0();
        U0_12 = quadAtom12->getU0();
        U0_22 = quadAtom22->getU0();

        B0_11 = quadAtom11->getX() + quadAtom11->getR0() * quadAtom11->getU0();
        B0_21 = quadAtom21->getX() + quadAtom21->getR0() * quadAtom21->getU0();
        B0_12 = quadAtom12->getX() + quadAtom12->getR0() * quadAtom12->getU0();
        B0_22 = quadAtom22->getX() + quadAtom22->getR0() * quadAtom22->getU0();

        S0_11 = U0_11 * quadAtom11->getR0();
        S0_21 = U0_21 * quadAtom21->getR0();
        S0_12 = U0_12 * quadAtom12->getR0();
        S0_22 = U0_22 * quadAtom22->getR0();

        r11 = quadAtom11->getR0();
        r21 = quadAtom21->getR0();
        r12 = quadAtom12->getR0();
        r22 = quadAtom22->getR0();

    }
    else
    {
        U0_11 = quadAtom11->getU1();
        U0_21 = quadAtom21->getU1();
        U0_12 = quadAtom12->getU1();
        U0_22 = quadAtom22->getU1();

        B0_11 = quadAtom11->getX() + quadAtom11->getR1() * quadAtom11->getU1();
        B0_21 = quadAtom21->getX() + quadAtom21->getR1() * quadAtom21->getU1();
        B0_12 = quadAtom12->getX() + quadAtom12->getR1() * quadAtom12->getU1();
        B0_22 = quadAtom22->getX() + quadAtom22->getR1() * quadAtom22->getU1();

        S0_11 = U0_11 * quadAtom11->getR1();
        S0_21 = U0_21 * quadAtom21->getR1();
        S0_12 = U0_12 * quadAtom12->getR1();
        S0_22 = U0_22 * quadAtom22->getR1();

        r11 = quadAtom11->getR1();
        r21 = quadAtom21->getR1();
        r12 = quadAtom12->getR1();
        r22 = quadAtom22->getR1();
    }

    Vector3D du11, du21, du12, du22, dv11, dv21, dv12, dv22;

    Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;
    Vector3D dB0du11, dB0dv11, dB0du21, dB0dv21, dB0du12, dB0dv12, dB0du22, dB0dv22;
    Vector3D dS0du11, dS0dv11, dS0du21, dS0dv21, dS0du12, dS0dv12, dS0du22, dS0dv22;

    dB0du11 = B0_21 - B0_11;
    dB0dv11 = B0_12 - B0_11;

    dB0du21 = B0_21 - B0_11;
    dB0dv21 = B0_22 - B0_21;

    dB0du12 = B0_22 - B0_12;
    dB0dv12 = B0_12 - B0_11;

    dB0du22 = B0_22 - B0_12;
    dB0dv22 = B0_22 - B0_21;

    du11 = x21 - x11;
    dv11 = x12 - x11;
    dU0du11 = U0_21 - U0_11;
    dU0dv11 = U0_12 - U0_11;
    dS0du11 = S0_21 - S0_11;
    dS0dv11 = S0_12 - S0_11;
    ru11 = r21 - r11;
    rv11 = r12 - r11;
//
//
    du21 = x21 - x11;
    dv21 = x22 - x21;
    dU0du21 = U0_21 - U0_11;
    dU0dv21 = U0_22 - U0_21;
    dS0du21 = S0_21 - S0_11;
    dS0dv21 = S0_22 - S0_21;
    ru21 = r21 - r11;
    rv21 = r22 - r21;
//
    du12 = x22 - x12;
    dv12 = x12 - x11;
    dU0du12 = U0_22 - U0_12;
    dU0dv12 = U0_12 - U0_11;
    dS0du12 = S0_22 - S0_12;
    dS0dv12 = S0_12 - S0_11;
    ru12 = r22 - r12;
    rv12 = r12 - r11;
//
    du22 = x22 - x12;
    dv22 = x22 - x21;
    dU0du22 = U0_22 - U0_12;
    dU0dv22 = U0_22 - U0_21;
    dS0du22 = S0_22 - S0_12;
    dS0dv22 = S0_22 - S0_21;
    ru22 = r22 - r12;
    rv22 = r22 - r21;

    if ( ubase != 0 )
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( ubase-1, vbase ) );
    }


    // Get unit normals for the four corners and project derivatives on to the tangent plane
    Vector3D n11 = quadAtom11->getU0() - quadAtom11->getU1();
    n11.normalize();
    Vector3D n21 = quadAtom21->getU0() - quadAtom21->getU1();
    n21.normalize();
    Vector3D n12 = quadAtom12->getU0() - quadAtom12->getU1();
    n12.normalize();
    Vector3D n22 = quadAtom22->getU0() - quadAtom22->getU1();
    n22.normalize();

    Vector3D du11t = du11 - ( du11 * n11 ) * n11;
    Vector3D du21t = du21 - ( du21 * n21 ) * n21;
    Vector3D du12t = du12 - ( du12 * n12 ) * n12;
    Vector3D du22t = du22 - ( du22 * n22 ) * n22;

    Vector3D dv11t = dv11 - ( dv11 * n11 ) * n11;
    Vector3D dv21t = dv21 - ( dv21 * n21 ) * n21;
    Vector3D dv12t = dv12 - ( dv12 * n12 ) * n12;
    Vector3D dv22t = dv22 - ( dv22 * n22 ) * n22;

    // Build matrices for hermite interpolation of medial sheet
    double hx[16] = { x11.getX(), x21.getX(), du11t.getX(), du21t.getX(),
                      x12.getX(), x22.getX(), du12t.getX(), du22t.getX(), dv11t.getX(),
                      dv21t.getX(), 0, 0, dv12t.getX(), dv22t.getX(), 0, 0
                    };

    double hy[16] = { x11.getY(), x21.getY(), du11t.getY(), du21t.getY(),
                      x12.getY(), x22.getY(), du12t.getY(), du22t.getY(), dv11t.getY(),
                      dv21t.getY(), 0, 0, dv12t.getY(), dv22t.getY(), 0, 0
                    };

    double hz[16] = { x11.getZ(), x21.getZ(), du11t.getZ(), du21t.getZ(),
                      x12.getZ(), x22.getZ(),du12t.getZ(), du22t.getZ(), dv11t.getZ(),
                      dv21t.getZ(), 0, 0, dv12t.getZ(), dv22t.getZ(), 0, 0
                    };

    double r_hermite[16] = { r11, r21, ru11, ru21, r12, r22, ru12, ru22, rv11, rv21, 0, 0, rv12, rv22, 0, 0};

    Matrix r_hermite_mat = Matrix ( 4,4,r_hermite, true );

    Matrix hxmat = Matrix ( 4, 4, hx, true );
    Matrix hymat = Matrix ( 4, 4, hy, true );
    Matrix hzmat = Matrix ( 4, 4, hz, true );

    //U0_11.print();

    double hu[4] = { h1 ( u ), h2 ( u ), h3 ( u ), h4 ( u ) };
    double hv[4] = { h1 ( v ), h2 ( v ), h3 ( v ), h4 ( v ) };
    Matrix humat = Matrix ( 1, 4, hu, true );
    Matrix hvmat = Matrix ( 4, 1, hv, true );

    Matrix xn = humat * hxmat * hvmat;
    Matrix yn = humat * hymat * hvmat;
    Matrix zn = humat * hzmat * hvmat;

    double Bx[16] = { B0_11.getX(), B0_21.getX(), dB0du11.getX(), dB0du21.getX(),
                      B0_12.getX(), B0_22.getX(), dB0du12.getX(), dB0du22.getX(), dB0dv11.getX(),
                      dB0dv21.getX(), 0, 0, dB0dv12.getX(), dB0dv22.getX(), 0, 0
                    };

    double By[16] = { B0_11.getY(), B0_21.getY(), dB0du11.getY(), dB0du21.getY(),
                      B0_12.getY(), B0_22.getY(), dB0du12.getY(), dB0du22.getY(), dB0dv11.getY(),
                      dB0dv21.getY(), 0, 0, dB0dv12.getY(), dB0dv22.getY(), 0, 0
                    };

    double Bz[16] = { B0_11.getZ(), B0_21.getZ(), dB0du11.getZ(), dB0du21.getZ(),
                      B0_12.getZ(), B0_22.getZ(), dB0du12.getZ(), dB0du22.getZ(), dB0dv11.getZ(),
                      dB0dv21.getZ(), 0, 0, dB0dv12.getZ(), dB0dv22.getZ(), 0, 0
                    };

    Matrix bxmat = Matrix ( 4, 4, Bx, true );
    Matrix bymat = Matrix ( 4, 4, By, true );
    Matrix bzmat = Matrix ( 4, 4, Bz, true );

    Matrix Bxn = humat * bxmat * hvmat;
    Matrix Byn = humat * bymat * hvmat;
    Matrix Bzn = humat * bzmat * hvmat;

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    Vector3D newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );
    Vector3D newbound = Vector3D ( Bxn ( 0,0 ), Byn ( 0,0 ), Bzn ( 0,0 ) );


    Vector3D htestD = newbound - newpos;
    double htestR = htestD.normalize();

    //M3DSpoke* hspoke = new M3DSpoke ( newpos, htestD, htestR ); // This is for testing purposes

    M3DSpoke hSpoke(newpos, htestD, htestR );

    return hSpoke;
}
