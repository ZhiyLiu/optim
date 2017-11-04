#ifndef _TSP_H_
#define _TSP_H_

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <assert.h>


class TSPNhd;

class SymmTable {
protected:
   double *data;

   int getIndex(int x,int y) {
      if(x<y)
         return getIndex(y,x);
      else if(x>y)
         return (y*size + x - (y*(y+3))/2 - 1);
      return -1;
   };
public:
   int size;

   SymmTable(int sideLength) {
      data = new double[sideLength*(sideLength-1)/2];
      size = sideLength;
   };

   double get(int x,int y) {
      return data[getIndex(x,y)];
   };
   void set(double value,int x,int y) {
      data[getIndex(x,y)] = value;
   }
   ~SymmTable() {
      delete data;
   };
};

class TSP {
protected:
   SymmTable *dist;

   void getDistIndex(int x,int y);
public:
   TSP(FILE *loadFile);
   double getDistance(int city1,int city2) {
	   return (city1==city2) ? 0 :dist->get(city1,city2);
   }
   int getSize() {
      return dist->size;
   }
   ~TSP() {
      delete dist;
   }
};

class TSPSolution {
protected:
   TSP *problem;
   int *cities;
   int size;

   double cost;
   void calcCost();
public:
   TSPSolution(TSP *problem);
   TSPSolution(TSP *problem,FILE *loadFile);
   TSPSolution(TSPSolution *soln);
   TSPSolution(TSP *problem,int cities[]);
   ~TSPSolution() {
      delete cities;
   }
   int getPrev(int city) {
      return city ? city-1 : size-1;
   }
   int getNext(int city) {
      return (city + 1)%size;
   }

   int getCityAt(int i) {
	   assert(i>=0&&i<size);
	   return cities[i];
   }

   int getSize() {
	   return size;
   }

   double getValue() {
      return cost;
   };

   double getTwoChangeDelta(int c1,int c2);
   TSPSolution *getBestTwoChange();
   TSPSolution *getRandomTwoChange();
   void applyTwoChange(int c1,int c2);
   void applyTwoChange(int c1,int c2,double bestDelta);
   TSPSolution *getLocalOptimum();
   TSPSolution *getNhdOptimum();
   TSPNhd * getNhd();

   void print(FILE *f) {
      fprintf(f,"NAME : PAUL'S SOLUTION.\n");
      fprintf(f,"COMMENT : COSTS %lg\n",cost);
      fprintf(f,"TYPE : TOUR\n");
      fprintf(f,"DIMENSION : %d\n",size);
      fprintf(f,"TOUR_SECTION\n");
      for(int i=0;i<size;i++) {
         fprintf(f,"%d\n",cities[i]+1);
      }
      fprintf(f,"-1\nEOF\n");
   }

   void citySwap(int c1,int c2);
   TSPSolution *getRandomCitySwap();


   int equals(TSPSolution *cmp);

   int nbindex;
};

#define UNDEFINED    0
#define IMMFAMILY       1
#define NONFAMILY    2

int tspSolutionSort(const void *a,const void *b);

class TSPNhd {
protected:
   TSPSolution *source;
   TSPSolution **nhd;
   int *familyStatus;
   int *familyMembers;
   int size;
   int family;
public:
   TSPNhd(TSPSolution *inSource,int inSize) {
      source = inSource;
      size = inSize;
      nhd = new TSPSolution*[size];
      memset(nhd,0,size*sizeof(TSPSolution *));

      familyStatus = new int[size];
      familyMembers=NULL;
      family=-1;
      memset(familyStatus,0,size*sizeof(int));
   }

   void setSolution(TSPSolution *soln,int index) {
      nhd[index] = soln;
   }

   TSPSolution *getSolution(int index) {
      return nhd[index];
   }

   ~TSPNhd() {
      for(int i=0;i<size;i++) {
         if(nhd[i])
            delete nhd[i];
      }
      delete nhd;
      delete familyStatus;
      if(familyMembers)
         delete familyMembers;
   }

   int getFamilyStatus(int index) {
      return familyStatus[index];
   }

   void calcFamilyStatus();
   int getSize() {
      return size;
   }

   int getFamilySize() {
      return family;
   }

   TSPSolution *getRandomFamilyMember() {
      return nhd[familyMembers[rand()%family]];
   }

   void sort() {
      qsort(nhd,size,sizeof(TSPSolution*),tspSolutionSort);
   }

};

#endif	/* _TSP_H_ */

