/**************************************************************************************************************

	int tri_tri_intersect_with_isectline_w_UVTRad(...)
	if the intersection exists, which is composed of 2 vertices, this routine 
	should return the UVT of the vertices in BOTH coordinate systems of the sub-figure and main-figure

***************************************************************************************************************/

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 * updated: 2001-06-20 (added line of intersection)
 *
 * int tri_tri_intersect(MyFloat V0[3],MyFloat V1[3],MyFloat V2[3],
 *                       MyFloat U0[3],MyFloat U1[3],MyFloat U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 * Here is a version withouts divisions (a little faster)
 * int NoDivTriTriIsect(MyFloat V0[3],MyFloat V1[3],MyFloat V2[3],
 *                      MyFloat U0[3],MyFloat U1[3],MyFloat U2[3]);
 * 
 * This version computes the line of intersection as well (if they are not coplanar):
 * int tri_tri_intersect_with_isectline(MyFloat V0[3],MyFloat V1[3],MyFloat V2[3], 
 *				        MyFloat U0[3],MyFloat U1[3],MyFloat U2[3],int *coplanar,
 *				        MyFloat isectpt1[3],MyFloat isectpt2[3]);
 * coplanar returns whether the tris are coplanar
 * isectpt1, isectpt2 are the endpoints of the line of intersection
 */

#include <math.h>

#ifndef MyFloat
#define MyFloat double
#endif

#define FABS(x) ((MyFloat)fabs(x))        /* implement as is fastest on your machine */

/* if USE_EPSILON_TEST is true then we do a check: 
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
#define USE_EPSILON_TEST TRUE  

// BUG FIXED
// in fact i really hesitated to call this a bug because it's more or less a numerical
// problem, so far it's temporarily 'solved'
//define EPSILON 0.000001
#define EPSILON 0.000000001

/* some macros */
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; dest[1]=v1[1]-v2[1]; dest[2]=v1[2]-v2[2]; 

#define ADD(dest,v1,v2) dest[0]=v1[0]+v2[0]; dest[1]=v1[1]+v2[1]; dest[2]=v1[2]+v2[2]; 

#define MULT(dest,v,factor) dest[0]=factor*v[0]; dest[1]=factor*v[1]; dest[2]=factor*v[2];

#define SET(dest,src) dest[0]=src[0]; dest[1]=src[1]; dest[2]=src[2]; 

/* sort so that a<=b */
#define SORT(a,b)       \
             if(a>b)    \
             {          \
               MyFloat c; \
               c=a;     \
               a=b;     \
               b=c;     \
             }

#define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
              isect1=VV0+(VV2-VV0)*D0/(D0-D2);


#define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
  }



/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */ 
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }                                

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  MyFloat Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  MyFloat a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}

int coplanar_tri_tri(MyFloat N[3],MyFloat V0[3],MyFloat V1[3],MyFloat V2[3],
                     MyFloat U0[3],MyFloat U1[3],MyFloat U2[3])
{
   MyFloat A[3];
   short i0,i1;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
   A[0]=fabs(N[0]);
   A[1]=fabs(N[1]);
   A[2]=fabs(N[2]);
   if(A[0]>A[1])
   {
      if(A[0]>A[2])  
      {
          i0=1;      /* A[0] is greatest */
          i1=2;
      }
      else
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
   }
   else   /* A[0]<=A[1] */
   {
      if(A[2]>A[1])
      {
          i0=0;      /* A[2] is greatest */
          i1=1;                                           
      }
      else
      {
          i0=0;      /* A[1] is greatest */
          i1=2;
      }
    }               

    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);

    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    POINT_IN_TRI(V0,U0,U1,U2);
    POINT_IN_TRI(U0,V0,V1,V2);

    return 0;
}


int tri_tri_intersect(MyFloat V0[3],MyFloat V1[3],MyFloat V2[3],
                      MyFloat U0[3],MyFloat U1[3],MyFloat U2[3])
{
  MyFloat E1[3],E2[3];
  MyFloat N1[3],N2[3],d1,d2;
  MyFloat du0,du1,du2,dv0,dv1,dv2;
  MyFloat D[3];
  MyFloat isect1[2], isect2[2];
  MyFloat du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  MyFloat vp0,vp1,vp2;
  MyFloat up0,up1,up2;
  MyFloat b,c,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(fabs(du0)<EPSILON) du0=0.0;
  if(fabs(du1)<EPSILON) du1=0.0;
  if(fabs(du2)<EPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
  if(fabs(dv0)<EPSILON) dv0=0.0;
  if(fabs(dv1)<EPSILON) dv1=0.0;
  if(fabs(dv2)<EPSILON) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
  max=fabs(D[0]);
  index=0;
  b=fabs(D[1]);
  c=fabs(D[2]);
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0[index];
  vp1=V1[index];
  vp2=V2[index];

  up0=U0[index];
  up1=U1[index];
  up2=U2[index];

  /* compute interval for triangle 1 */
  COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);

  /* compute interval for triangle 2 */
  COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);

  SORT(isect1[0],isect1[1]);
  SORT(isect2[0],isect2[1]);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}


#define NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1) \
{ \
        if(D0D1>0.0f) \
        { \
                /* here we know that D0D2<=0.0 */ \
            /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
        } \
        else if(D0D2>0.0f)\
        { \
                /* here we know that d0d1<=0.0 */ \
            A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
        } \
        else if(D1*D2>0.0f || D0!=0.0f) \
        { \
                /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
                A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
        } \
        else if(D1!=0.0f) \
        { \
                A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
        } \
        else if(D2!=0.0f) \
        { \
                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
        } \
        else \
        { \
                /* triangles are coplanar */ \
                return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
        } \
}



int NoDivTriTriIsect(MyFloat V0[3],MyFloat V1[3],MyFloat V2[3],
                     MyFloat U0[3],MyFloat U1[3],MyFloat U2[3])
{
  MyFloat E1[3],E2[3];
  MyFloat N1[3],N2[3],d1,d2;
  MyFloat du0,du1,du2,dv0,dv1,dv2;
  MyFloat D[3];
  MyFloat isect1[2], isect2[2];
  MyFloat du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  MyFloat vp0,vp1,vp2;
  MyFloat up0,up1,up2;
  MyFloat bb,cc,max;
  MyFloat a,b,c,x0,x1;
  MyFloat d,e,f,y0,y1;
  MyFloat xx,yy,xxyy,tmp;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(FABS(du0)<EPSILON) du0=0.0;
  if(FABS(du1)<EPSILON) du1=0.0;
  if(FABS(du2)<EPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
  if(FABS(dv0)<EPSILON) dv0=0.0;
  if(FABS(dv1)<EPSILON) dv1=0.0;
  if(FABS(dv2)<EPSILON) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
  max=(MyFloat)FABS(D[0]);
  index=0;
  bb=(MyFloat)FABS(D[1]);
  cc=(MyFloat)FABS(D[2]);
  if(bb>max) max=bb,index=1;
  if(cc>max) max=cc,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0[index];
  vp1=V1[index];
  vp2=V2[index];

  up0=U0[index];
  up1=U1[index];
  up2=U2[index];

  /* compute interval for triangle 1 */
  NEWCOMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,a,b,c,x0,x1);

  /* compute interval for triangle 2 */
  NEWCOMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,d,e,f,y0,y1);

  xx=x0*x1;
  yy=y0*y1;
  xxyy=xx*yy;

  tmp=a*xxyy;
  isect1[0]=tmp+b*x1*yy;
  isect1[1]=tmp+c*x0*yy;

  tmp=d*xxyy;
  isect2[0]=tmp+e*xx*y1;
  isect2[1]=tmp+f*xx*y0;

  SORT(isect1[0],isect1[1]);
  SORT(isect2[0],isect2[1]);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}

/* sort so that a<=b */
#define SORT2(a,b,smallest)       \
             if(a>b)       \
             {             \
               MyFloat c;    \
               c=a;        \
               a=b;        \
               b=c;        \
               smallest=1; \
             }             \
             else smallest=0;


inline void isect2(MyFloat VTX0[3],MyFloat VTX1[3],MyFloat VTX2[3],MyFloat VV0,MyFloat VV1,MyFloat VV2,
	    MyFloat D0,MyFloat D1,MyFloat D2,MyFloat *isect0,MyFloat *isect1,MyFloat isectpoint0[3],MyFloat isectpoint1[3]) 
{
  MyFloat tmp=D0/(D0-D1);          
  MyFloat diff[3];
  *isect0=VV0+(VV1-VV0)*tmp;         
  SUB(diff,VTX1,VTX0);              
  MULT(diff,diff,tmp);               
  ADD(isectpoint0,diff,VTX0);        
  tmp=D0/(D0-D2);                    
  *isect1=VV0+(VV2-VV0)*tmp;          
  SUB(diff,VTX2,VTX0);                   
  MULT(diff,diff,tmp);                 
  ADD(isectpoint1,VTX0,diff);          
}


#if 0
#define ISECT2(VTX0,VTX1,VTX2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1) \
              tmp=D0/(D0-D1);                    \
              isect0=VV0+(VV1-VV0)*tmp;          \
	      SUB(diff,VTX1,VTX0);               \
	      MULT(diff,diff,tmp);               \
              ADD(isectpoint0,diff,VTX0);        \
              tmp=D0/(D0-D2);                    
/*              isect1=VV0+(VV2-VV0)*tmp;          \ */
/*              SUB(diff,VTX2,VTX0);               \     */
/*              MULT(diff,diff,tmp);               \   */
/*              ADD(isectpoint1,VTX0,diff);           */
#endif

inline int compute_intervals_isectline(MyFloat VERT0[3],MyFloat VERT1[3],MyFloat VERT2[3],
				       MyFloat VV0,MyFloat VV1,MyFloat VV2,MyFloat D0,MyFloat D1,MyFloat D2,
				       MyFloat D0D1,MyFloat D0D2,MyFloat *isect0,MyFloat *isect1,
				       MyFloat isectpoint0[3],MyFloat isectpoint1[3])
{
  if(D0D1>0.0f)                                        
  {                                                    
    /* here we know that D0D2<=0.0 */                  
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1);
  } 
  else if(D0D2>0.0f)                                   
    {                                                   
    /* here we know that d0d1<=0.0 */             
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1);
  }                                                  
  else if(D1*D2>0.0f || D0!=0.0f)   
  {                                   
    /* here we know that d0d1<=0.0 or that D0!=0.0 */
    isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1);   
  }                                                  
  else if(D1!=0.0f)                                  
  {                                               
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1); 
  }                                         
  else if(D2!=0.0f)                                  
  {                                                   
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1);     
  }                                                 
  else                                               
  {                                                   
    /* triangles are coplanar */    
    return 1;
  }
  return 0;
}

#define COMPUTE_INTERVALS_ISECTLINE(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1,isectpoint0,isectpoint1) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     
#if 0
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
    coplanar=1;                                         \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
  }
#endif

int tri_tri_intersect_with_isectline(MyFloat V0[3],MyFloat V1[3],MyFloat V2[3],
				     MyFloat U0[3],MyFloat U1[3],MyFloat U2[3],int *coplanar,
				     MyFloat isectpt1[3],MyFloat isectpt2[3])
{
  MyFloat E1[3],E2[3];
  MyFloat N1[3],N2[3],d1,d2;
  MyFloat du0,du1,du2,dv0,dv1,dv2;
  MyFloat D[3];
  MyFloat isect1[2], isect2[2];
  MyFloat isectpointA1[3],isectpointA2[3];
  MyFloat isectpointB1[3],isectpointB2[3];
  MyFloat du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  MyFloat vp0,vp1,vp2;
  MyFloat up0,up1,up2;
  MyFloat b,c,max;
  // MyFloat tmp,diff[3];
  int smallest1,smallest2;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(fabs(du0)<EPSILON) du0=0.0;
  if(fabs(du1)<EPSILON) du1=0.0;
  if(fabs(du2)<EPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
  if(fabs(dv0)<EPSILON) dv0=0.0;
  if(fabs(dv1)<EPSILON) dv1=0.0;
  if(fabs(dv2)<EPSILON) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
  max=fabs(D[0]);
  index=0;
  b=fabs(D[1]);
  c=fabs(D[2]);
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0[index];
  vp1=V1[index];
  vp2=V2[index];

  up0=U0[index];
  up1=U1[index];
  up2=U2[index];

  /* compute interval for triangle 1 */
  *coplanar=compute_intervals_isectline(V0,V1,V2,vp0,vp1,vp2,dv0,dv1,dv2,
				       dv0dv1,dv0dv2,&isect1[0],&isect1[1],isectpointA1,isectpointA2);
  if(*coplanar) return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);     


  /* compute interval for triangle 2 */
  compute_intervals_isectline(U0,U1,U2,up0,up1,up2,du0,du1,du2,
			      du0du1,du0du2,&isect2[0],&isect2[1],isectpointB1,isectpointB2);

  SORT2(isect1[0],isect1[1],smallest1);
  SORT2(isect2[0],isect2[1],smallest2);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;

  /* at this point, we know that the triangles intersect */

  if(isect2[0]<isect1[0])
  {
    if(smallest1==0) { SET(isectpt1,isectpointA1); }
    else { SET(isectpt1,isectpointA2); }

    if(isect2[1]<isect1[1])
    {
      if(smallest2==0) { SET(isectpt2,isectpointB2); }
      else { SET(isectpt2,isectpointB1); }
    }
    else
    {
      if(smallest1==0) { SET(isectpt2,isectpointA2); }
      else { SET(isectpt2,isectpointA1); }
    }
  }
  else
  {
    if(smallest2==0) { SET(isectpt1,isectpointB1); }
    else { SET(isectpt1,isectpointB2); }

    if(isect2[1]>isect1[1])
    {
      if(smallest1==0) { SET(isectpt2,isectpointA2); }
      else { SET(isectpt2,isectpointA1); }      
    }
    else
    {
      if(smallest2==0) { SET(isectpt2,isectpointB2); }
      else { SET(isectpt2,isectpointB1); } 
    }
  }

  return 1;
}

MyFloat Distance(Bpoint v0, Bpoint v1)
{
	MyFloat t[3];
	for(int i=0; i<3; i++)
		t[i]=v0.pnt[i]-v1.pnt[i];
	return sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
}

MyFloat Distance(Bpoint2 *v0, Bpoint v1)
{
	MyFloat t[3];
	for(int i=0; i<3; i++)
		t[i]=v0->pnt[i]-v1.pnt[i];
	return sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
}

// if the intersection exists, which is composed of 2 vertices
// this routine should return the UVT of the vertices in BOTH of the figures' coordinates
int tri_tri_intersect_with_isectline_w_UVTRad(Bpoint v0, Bpoint v1, Bpoint v2, Bpoint u0, Bpoint u1, Bpoint u2, \
											  int *coplanar, Bpoint2 *pt1, Bpoint2 *pt2)//MyFloat isectpt1[3], MyFloat isectpt2[3])
{
  MyFloat V0[3], V1[3], V2[3];
  MyFloat U0[3], U1[3], U2[3];
  MyFloat isectpt1[3], isectpt2[3];

  MyFloat E1[3],E2[3];
  MyFloat N1[3],N2[3],d1,d2;
  MyFloat du0,du1,du2,dv0,dv1,dv2;
  MyFloat D[3];
  MyFloat isect1[2], isect2[2];
  MyFloat isectpointA1[3],isectpointA2[3];
  MyFloat isectpointB1[3],isectpointB2[3];
  MyFloat du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  MyFloat vp0,vp1,vp2;
  MyFloat up0,up1,up2;
  MyFloat b,c,max;
  // MyFloat tmp,diff[3];
  int smallest1,smallest2;


  int i;
  for(i=0; i<3; i++)
  {
	  V0[i]=v0.pnt[i];
	  V1[i]=v1.pnt[i];
	  V2[i]=v2.pnt[i];
	  U0[i]=u0.pnt[i];
	  U1[i]=u1.pnt[i];
	  U2[i]=u2.pnt[i];
  }

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(fabs(du0)<EPSILON) du0=0.0;
  if(fabs(du1)<EPSILON) du1=0.0;
  if(fabs(du2)<EPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
  if(fabs(dv0)<EPSILON) dv0=0.0;
  if(fabs(dv1)<EPSILON) dv1=0.0;
  if(fabs(dv2)<EPSILON) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
  max=fabs(D[0]);
  index=0;
  b=fabs(D[1]);
  c=fabs(D[2]);
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0[index];
  vp1=V1[index];
  vp2=V2[index];

  up0=U0[index];
  up1=U1[index];
  up2=U2[index];

  /* compute interval for triangle 1 */
  *coplanar=compute_intervals_isectline(V0,V1,V2,vp0,vp1,vp2,dv0,dv1,dv2,
				       dv0dv1,dv0dv2,&isect1[0],&isect1[1],isectpointA1,isectpointA2);
  if(*coplanar) return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);     


  /* compute interval for triangle 2 */
  compute_intervals_isectline(U0,U1,U2,up0,up1,up2,du0,du1,du2,
			      du0du1,du0du2,&isect2[0],&isect2[1],isectpointB1,isectpointB2);

  SORT2(isect1[0],isect1[1],smallest1);
  SORT2(isect2[0],isect2[1],smallest2);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;

  /* at this point, we know that the triangles intersect */

  if(isect2[0]<isect1[0])
  {
    if(smallest1==0) { SET(isectpt1,isectpointA1); }
    else { SET(isectpt1,isectpointA2); }

    if(isect2[1]<isect1[1])
    {
      if(smallest2==0) { SET(isectpt2,isectpointB2); }
      else { SET(isectpt2,isectpointB1); }
    }
    else
    {
      if(smallest1==0) { SET(isectpt2,isectpointA2); }
      else { SET(isectpt2,isectpointA1); }
    }
  }
  else
  {
    if(smallest2==0) { SET(isectpt1,isectpointB1); }
    else { SET(isectpt1,isectpointB2); }

    if(isect2[1]>isect1[1])
    {
      if(smallest1==0) { SET(isectpt2,isectpointA2); }
      else { SET(isectpt2,isectpointA1); }      
    }
    else
    {
      if(smallest2==0) { SET(isectpt2,isectpointB2); }
      else { SET(isectpt2,isectpointB1); } 
    }
  }

  for(i=0; i<3; i++)
  {
	  pt1->pnt[i]=isectpt1[i];
	  pt2->pnt[i]=isectpt2[i];
  }


  MyFloat ed0, ed1, ed2, eds;
  MyFloat area0, area;
  MyFloat U, V, W, U_2, V_2, W_2;

  // UVT for first intersection point in the triangle of main-figure
  ed0=Distance(v0, v1);
  ed1=Distance(v1, v2);
  ed2=Distance(v2, v0);
  eds=(ed0+ed1+ed2)/2;
  //area0=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area0=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area0>=0)
	  area0=sqrt(area0);
  else
  {
	  area0=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area0 is INVALID!\n");
  }
  ed0=Distance(pt1, v1);
  ed2=Distance(pt1, v2);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }
  U=area/area0;
  if(U>1) U=1;
  ed0=Distance(pt2, v1);
  ed2=Distance(pt2, v2);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }
  U_2=area/area0;
  if(U_2>1) U_2=1;

  ed0=Distance(pt1, v2);
  ed1=Distance(v2, v0);
  ed2=Distance(pt1, v0);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }
  V=area/area0;
  if(V>1) V=1;
  ed0=Distance(pt2, v2);
  ed1=Distance(v2, v0);
  ed2=Distance(pt2, v0);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }

  V_2=area/area0;
  if(V_2>1) V_2=1;

  W=1-U-V;
  if(W<0) W=0;
  W_2=1-U_2-V_2;
  if(W_2<0) W_2=0;

  pt1->u=U*v0.u+V*v1.u+W*v2.u;
  pt1->v=U*v0.v+V*v1.v+W*v2.v;
  pt1->t=U*v0.t+V*v1.t+W*v2.t;
  pt1->rad=U*v0.rad+V*v1.rad+W*v2.rad;
  pt2->u=U_2*v0.u+V_2*v1.u+W_2*v2.u;
  pt2->v=U_2*v0.v+V_2*v1.v+W_2*v2.v;
  pt2->t=U_2*v0.t+V_2*v1.t+W_2*v2.t;
  pt2->rad=U_2*v0.rad+V_2*v1.rad+W_2*v2.rad;

  // UVT for first intersection point in the triangle of sub-figure
  ed0=Distance(u0, u1);
  ed1=Distance(u1, u2);
  ed2=Distance(u2, u0);
  eds=(ed0+ed1+ed2)/2;
  //area0=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area0=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area0>=0)
	  area0=sqrt(area0);
  else
  {
	  area0=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area0 is INVALID!\n");
  }
  ed0=Distance(pt1, u1);
  ed2=Distance(pt1, u2);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }
  U=area/area0;
  if(U>1) U=1;
  ed0=Distance(pt2, u1);
  ed2=Distance(pt2, u2);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }
  U_2=area/area0;
  if(U_2>1) U_2=1;

  ed0=Distance(pt1, u2);
  ed1=Distance(u2, u0);
  ed2=Distance(pt1, u0);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }
  V=area/area0;
  if(V>1) V=1;
  ed0=Distance(pt2, u2);
  ed1=Distance(u2, u0);
  ed2=Distance(pt2, u0);
  eds=(ed0+ed1+ed2)/2;
  //area=sqrt(eds*(eds-ed0)*(eds-ed1)*(eds-ed2));
  area=eds*(eds-ed0)*(eds-ed1)*(eds-ed2);
  if(area>=0)
	  area=sqrt(area);
  else
  {
	  area=0;
	  //printf("In tri_tri_intersect_with_isectline_w_UVTRad(), the area is INVALID!\n");
  }
  V_2=area/area0;
  if(V_2>1) V_2=1;

  W=1-U-V;
  if(W<0) W=0;
  W_2=1-U_2-V_2;
  if(W_2<0) W_2=0;

  pt1->u_2=U*u0.u+V*u1.u+W*u2.u;
  pt1->v_2=U*u0.v+V*u1.v+W*u2.v;
  pt1->t_2=U*u0.t+V*u1.t+W*u2.t;
  pt1->rad_2=U*u0.rad+V*u1.rad+W*u2.rad;
  pt2->u_2=U_2*u0.u+V_2*u1.u+W_2*u2.u;
  pt2->v_2=U_2*u0.v+V_2*u1.v+W_2*u2.v;
  pt2->t_2=U_2*u0.t+V_2*u1.t+W_2*u2.t;
  pt2->rad_2=U_2*u0.rad+V_2*u1.rad+W_2*u2.rad;

  /*
  pt1->pnt[0]=isectpt1[0];
  pt1->pnt[1]=isectpt1[1];
  pt1->pnt[2]=isectpt1[2];
  pt2->pnt[0]=isectpt2[0];
  pt2->pnt[1]=isectpt2[1];
  pt2->pnt[2]=isectpt2[2];
  */

  return 1;
}

