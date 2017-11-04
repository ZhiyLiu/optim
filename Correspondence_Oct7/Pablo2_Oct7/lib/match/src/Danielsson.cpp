#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif
#include <iostream>

// BINARY corrects for non-cubic voxels but grayscale does not,
// apparently because distance maps are for optimization in BINARY
// and for multi-figure blending in non-BINARY (always cubic).
// I'd prefer a switch that tells if non-cubic corrections should be
// applied or not, rather than 2 APIs.
// -GST 20080111

#ifdef BINARY
// Req'd for anisotropic voxels
#define DIS(x,y,z) ((int) (sqtab_scaledx[x] + sqtab_scaledy[y] + sqtab_scaledz[z]))
#else
#define SQ(x) sqtab[x]
#define DIS(x,y,z) (SQ(x)+SQ(y)+SQ(z))
#endif

// finds the maximum of 3 numbers
#define MAX(x,y) (((x) > (y)) ? (x) : (y)), (z)
#define MAX3(x,y,z) MAX((((x) > (y)) ? (x) : (y)), (z))


/* Extension of  Danielsson algorithm  to 3D */
/* CGIP 14, 227-248			     */
/* From Martin Styner			     */
/* Modified by Sarang Joshi to handle char input */


#ifdef BINARY
int edt3ddan(char * bdryIm,		// 1 at boundary of object; 0 elsewhere
			 int xdim, int ydim, int zdim,	// for bdryIm
			 double spacingx, double spacingy, double spacingz,	// for bdryIm
			 unsigned char metric,	// 6 or 26: # of neighbor distances to check
			 short ** odx,	// Output Distance map images, per direction
			 short ** ody,
			 short ** odz)
#else
int edt3ddan(char * bdryIm,
			 int xdim, int ydim, int zdim,
			 unsigned char metric,
			 short ** odx, short ** ody, short ** odz)
#endif
{
  int x,y,z,index,indexc, xydim;
  int *sqtab;	// lookup table of squared signed integers: sqtab[i] = i*i
  int *sq;
  int maxdim;	// max possible voxel distance in any direction
  short *dx, *dy, *dz;
  int stop, upy;
  //int part;	// DEBUG
  int pass = 1;		// print status
  int dist, testDist;
#ifdef BINARY
  double spacingx2 = spacingx * spacingx;
  double spacingy2 = spacingy * spacingy;
  double spacingz2 = spacingz * spacingz;
#endif

  maxdim = 4 * MAX3(xdim,ydim,zdim);

#ifdef BINARY
  // No longer req'd, but I'm not sure what free sq will do if sq never got malloc'd
#endif
  if(!(sq=sqtab=(int *)malloc((4*maxdim+1)*sizeof(int)))) {
    return(-1);
  }
  sqtab = &sqtab[2*maxdim];
  for(index=(-2*maxdim);index<=(2*maxdim);index++) {
     sqtab[index] = index * index;
  }

#ifdef BINARY
  // For anisotropic distances, use 3 square tables, each one scaled so that
  // look ups in x, y, z returns the correct square distance.
  double * sqtab_scaledx = new double[4*maxdim+1];
  double * sqtab_scaledy = new double[4*maxdim+1];
  double * sqtab_scaledz = new double[4*maxdim+1];
  sqtab_scaledx = &sqtab_scaledx[2*maxdim];		// Shift the arrays forward by half their length
  sqtab_scaledy = &sqtab_scaledy[2*maxdim];		// so that we can index negative numbers
  sqtab_scaledz = &sqtab_scaledz[2*maxdim];

  for ((index=-2*maxdim); index<=(2*maxdim); index++) {
		double index_sq = index*index;
		sqtab_scaledx[index] = index_sq * spacingx2;
		sqtab_scaledy[index] = index_sq * spacingy2;
		sqtab_scaledz[index] = index_sq * spacingz2;
  }
#endif

  if(!(dx=(short *)malloc((xdim*ydim*zdim)*sizeof(short)))) {
    free(sq);
    return(-2);
  }
  if(!(dy=(short *)malloc((xdim*ydim*zdim)*sizeof(short)))) {
    free(dx);
    free(sq);
    return(-3);
  }
  if(!(dz=(short *)malloc((xdim*ydim*zdim)*sizeof(short)))) {
    free(dy);
    free(dx);
    free(sq);
    return(-4);
  }
  printf("START: ");

  printf("%d ", pass++);
  xydim = xdim * ydim;
  for(x=0;x<xdim;x++) {
    for(y=0;y<ydim;y++) {
      for(z=0;z<zdim;z++) {
        index = x + xdim * y + xydim * z;
        if(bdryIm[index]) {
          dx[index] = dy[index] = dz[index] = 0;
        }
        else {
          dx[index] = dy[index] = dz[index] = maxdim;
        }
      }
    }
  }

  printf("%d ", pass++);
  for(z=1;z<zdim;z++) {
    //part = (int) (20.0 + 25.0 * ((float) z - 1.0) / ((float) zdim - 1.0));
    upy = ydim - 1;
    stop = xdim - 1;
    for(y=0;y<ydim;y++) {
      index = z * xydim + y * xdim;
      for(x=0;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - xydim;
        if((testDist=DIS(dx[indexc],dy[indexc],dz[indexc]-1)) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc];
          dz[index] = dz[indexc] - 1;
          dist = testDist;
        }
        if(metric != 6) {
          if(y > 0) {
            indexc = index - xdim - xydim;
            if((testDist=DIS(dx[indexc],dy[indexc]-1,dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc] - 1;
              dist = testDist;
            }
          }
          if(y < upy) {
            indexc = index + xdim - xydim;
            if((testDist=DIS(dx[indexc],dy[indexc]+1,dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc] - 1;
              dist = testDist;
            }
          }
          if(x > 0) {
            indexc = index - 1 - xydim;
            if((testDist=DIS(dx[indexc]-1,dy[indexc],dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] - 1;
              dist = testDist;
            }
          }
          if(x < stop) {
            indexc = index + 1 - xydim;
            if((testDist=DIS(dx[indexc]+1,dy[indexc],dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] - 1;
              dist = testDist;
            }
          }
          if(metric == 26) {
            if(y > 0) {
              if(x > 0) {
                indexc = index - xydim - xdim - 1;
                if((testDist=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] - 1;
                  dist = testDist;
                }
              }
              if(x < stop) {
                indexc = index - xydim - xdim + 1;    
                if((testDist=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] - 1;
                  dist = testDist;
                }
              }
            }
            if(y < upy) {
              if(x > 0) {
                indexc = index - xydim + xdim - 1;
                if((testDist=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] - 1;
                  dist = testDist;
                }
              }
              if(x < stop) {
                indexc = index - xydim + xdim + 1;
                if((testDist=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] - 1;
                  dist = testDist;
                }
              }
            }
          }
        }
      }
    }
    for(y=1;y<ydim;y++) {
      index = z * xydim + y * xdim;
      for(x=0;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - xdim;
        if((testDist=DIS(dx[indexc],dy[indexc]-1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] - 1;
          dz[index] = dz[indexc];
          dist = testDist;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index - xdim - 1;
            if((testDist=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
          if(x < stop) {
            indexc = index - xdim + 1;
            if((testDist=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
        }
      }
      index = z * xydim + y * xdim + 1;
      for(x=1;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((testDist=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
      index = z * xydim + (y + 1) * xdim - 2;
      for(x=(xdim-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((testDist=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
    }
    for(y=(ydim-2);y>=0;y--) {
      index = z * xydim + y * xdim;
      for(x=0;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + xdim;
        if((testDist=DIS(dx[indexc],dy[indexc]+1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] + 1;
          dz[index] = dz[indexc];
          dist = testDist;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index + xdim - 1;
            if((testDist=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
          if(x < stop) {
            indexc = index + xdim + 1;
            if((testDist=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
        }
      }
      index = z * xydim + y * xdim + 1;
      for(x=1;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((testDist=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
      index = z * xydim + (y + 1) * xdim - 2;
      for(x=(xdim-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((testDist=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
    }
  }

  printf("%d ", pass++);
  for(z=zdim-2;z>=0;z--) {
    //part = (int) (45.0 + 25.0 * ((float) (zdim - 2 - z) / ((float) zdim - 1.0)));
    upy = ydim - 1;
    stop = xdim - 1;
    for(y=0;y<ydim;y++) {
      index = z * xydim + y * xdim;
      for(x=0;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + xydim;
        if((testDist=DIS(dx[indexc],dy[indexc],dz[indexc]+1)) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc];
          dz[index] = dz[indexc] + 1;
          dist = testDist;
        }
        if(metric != 6) {
          if(y > 0) {
            indexc = index - xdim + xydim;
            if((testDist=DIS(dx[indexc],dy[indexc]-1,dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc] + 1;
              dist = testDist;
            }
          }
          if(y < upy) {
            indexc = index + xdim + xydim;
            if((testDist=DIS(dx[indexc],dy[indexc]+1,dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc] + 1;
              dist = testDist;
            }
          }
          if(x > 0) {
            indexc = index - 1 + xydim;
            if((testDist=DIS(dx[indexc]-1,dy[indexc],dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] + 1;
              dist = testDist;
            }
          }
          if(x < stop) {
            indexc = index + 1 + xydim;
            if((testDist=DIS(dx[indexc]+1,dy[indexc],dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] + 1;
              dist = testDist;
            }
          }
          if(metric == 26) {
            if(y > 0) {
              if(x > 0) {
                indexc = index + xydim - xdim - 1;
                if((testDist=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] + 1;
                  dist = testDist;
                }
              }
              if(x < stop) {
                indexc = index + xydim - xdim + 1;    
                if((testDist=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] + 1;
                  dist = testDist;
                }
              }
            }
            if(y < upy) {
              if(x > 0) {
                indexc = index + xydim + xdim - 1;
                if((testDist=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] + 1;
                  dist = testDist;
                }
              }
              if(x < stop) {
                indexc = index + xydim + xdim + 1;
                if((testDist=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] + 1;
                  dist = testDist;
                }
              }
            }
          }
        }
      }
    }
    for(y=1;y<ydim;y++) {
      index = z * xydim + y * xdim;
      for(x=0;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - xdim;
        if((testDist=DIS(dx[indexc],dy[indexc]-1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] - 1;
          dz[index] = dz[indexc];
          dist = testDist;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index - xdim - 1;
            if((testDist=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
          if(x < stop) {
            indexc = index - xdim + 1;
            if((testDist=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
        }
      }
      index = z * xydim + y * xdim + 1;
      for(x=1;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((testDist=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
      index = z * xydim + (y + 1) * xdim - 2;
      for(x=(xdim-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((testDist=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
    }
    for(y=(ydim-2);y>=0;y--) {
      index = z * xydim + y * xdim;
      for(x=0;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + xdim;
        if((testDist=DIS(dx[indexc],dy[indexc]+1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] + 1;
          dz[index] = dz[indexc];
          dist = testDist;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index + xdim - 1;
            if((testDist=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
          if(x < stop) {
            indexc = index + xdim + 1;
            if((testDist=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = testDist;
            }
          }
        }
      }
      index = z * xydim + y * xdim + 1;
      for(x=1;x<xdim;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((testDist=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
      index = z * xydim + (y + 1) * xdim - 2;
      for(x=(xdim-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((testDist=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = testDist;
        }
      }
    }
  }

  printf("%d ", pass++);
  for(z=0;z<zdim;z++) {
    for(y=0;y<ydim;y++) {
      index = z * xydim + y * xdim;
      for(x=0;x<xdim;x++,index++) {
        if(((x+dx[index]) < 0) || ((x+dx[index]) > (xdim - 1)) ||
           ((y+dy[index]) < 0) || ((y+dy[index]) > (ydim - 1)) ||
           ((z+dz[index]) < 0) || ((z+dz[index]) > (zdim - 1))) {
          printf("Error: Unclassified point [%d,%d,%d]\n",x,y,z);
          return(-111);
        }
        if(!bdryIm[index]) {
          bdryIm[index] = bdryIm[ (z+dz[index])*xydim
			                     +(y+dy[index])*xdim
                                 +(x+dx[index])];
        }
      }
    }
  }

  printf("done\n");
  *odx = dx;
  *ody = dy;
  *odz = dz;
  free(sq);
  return(0);
}

