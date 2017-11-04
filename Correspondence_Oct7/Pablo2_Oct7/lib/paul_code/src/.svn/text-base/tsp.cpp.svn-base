#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "tsp.h"




/**
 * Class TSP
 */
TSP::TSP(FILE *file) {
   char buffer[200];
   int dim = 0;

   while(fgets(buffer,200,file)) {
      if(!strncmp(buffer,"DIMENSION",strlen("DIMENSION"))) {
         char *colon = strchr(buffer,':');
         dim = atoi(colon + strspn(colon,": \t"));
      }
      else if (!strncmp(buffer,"NODE_COORD_SECTION",strlen("NODE_COORD_SECTION"))) {
         break;
      }
   }

   if(!dim) {
      printf("No dimension in problem\n");
      fclose(file);
      return;
   }


   int *x = new int[dim];
   int *y = new int[dim];

   for(int i=0;i<dim;i++) {
      fgets(buffer,200,file);
      // First value is ignored
      strtok(buffer," \n\t");

      // Second value is X
      x[i] = atoi(strtok(NULL," \n\t"));

      // Third value is y
      y[i] = atoi(strtok(NULL," \n\t"));
   }

   // Create dimension array
   dist = new SymmTable(dim);
   for(int xi=0;xi<dim-1;xi++) {
      for(int yi=xi+1;yi<dim;yi++) {
	     double temp = (x[xi]-x[yi])*(x[xi]-x[yi])+
                            (y[xi]-y[yi])*(y[xi]-y[yi]);
         double dst = sqrt(temp);
         dist->set(dst,xi,yi);
      }
   }

   delete x;
   delete y;

   fclose(file);
}

/**
 * Class TSP Solution
 */
void TSPSolution::calcCost() {
   cost = 0;
   for(int i=0;i<size;i++) {
      cost+=problem->getDistance(cities[i],cities[(i+1)%size]);
   }
   //if(globalOpt && (cost - globalOpt->getValue() < -0.1))
   //   print();
}

TSPSolution::TSPSolution(TSP *inProblem) {
   problem = inProblem;
   cities = new int[problem->getSize()];
   size = problem->getSize();

   int *taken = new int[size];
   for (int i=0;i<size;i++) {
      taken[i] = 0;
   }

   // Create cities
   for (int i=0;i<size;i++) {
      int rnd = rand()%size;
      while(taken[rnd]) {
         rnd = (rnd + 1) % size;
      }
      taken[rnd] = 1;
      cities[i] = rnd;
   }

   delete taken;
   calcCost();
}

TSPSolution::TSPSolution(TSPSolution *copy) {
   cost = copy->cost;
   size = copy->size;
   problem = copy->problem;

   cities = new int[size];
   memcpy(cities,copy->cities,sizeof(int)*size);
}

TSPSolution::TSPSolution(TSP *inProblem,int inCities[]) {
   problem = inProblem;
   size = problem->getSize();
   cities = new int[size];
   memcpy(cities,inCities,sizeof(int)*size);
   calcCost();
}

TSPSolution::TSPSolution(TSP *inProblem, FILE *f) {
   char buffer[200];
   while(fgets(buffer,200,f)) {
      if (!strncmp(buffer,"TOUR_SECTION",strlen("TOUR_SECTION"))) {
         break;
      }
   }

   if(feof(f)) {
      printf("Not a valid problem.\n");
      return;
   }

   problem = inProblem;
   cities = new int[problem->getSize()];
   size = problem->getSize();

   for(int i=0;i<size;i++) {
      fgets(buffer,200,f);
      cities[i] = atoi(buffer)-1;
   }

   calcCost();
   fclose(f);
};

double TSPSolution::getTwoChangeDelta(int c1,int c2) {
   // Illegal case: c1-c2 is an edge, can't swap edge w/self
   if(c2==getNext(c1))
      return 0;

   // This is nice & easy for symmetric
   double delta = 0;

   delta -= problem->getDistance(cities[c1],cities[getNext(c1)]);
   delta -= problem->getDistance(cities[getPrev(c2)],cities[c2]);
   delta += problem->getDistance(cities[c1],cities[getPrev(c2)]);
   delta += problem->getDistance(cities[getNext(c1)],cities[c2]);

   return delta;
}

// Swaps edge beginning with c1 with edge ending w. c2
// (reverses path between c1&c2
void TSPSolution::applyTwoChange(int c1,int c2,double bestDelta) {

   int x=c1,y=c2;
   int oldx,temp;

   while(1) {
      oldx = x;
      x = getNext(x);
      y = getPrev(y);
      if((x==y) || (oldx==y))
         break;

      // Swap cities x & y
      temp=cities[x];
      cities[x]=cities[y];
      cities[y]=temp;
   }

   // Adjust the cost
   cost+=bestDelta;

   //if(globalOpt && (cost - globalOpt->getValue() < -0.1))
   //   print();
}

// Poor man's variation
void TSPSolution::applyTwoChange(int c1,int c2) {
   applyTwoChange(c1,c2,getTwoChangeDelta(c1,c2));
}

TSPSolution *TSPSolution::getBestTwoChange() {

   double bestDelta = 0;
   int best1=0,best2=0;

   // The diagram below shows which combinations are needed for 2-ch
   // neigbourhood checking
   // x - illegal
   // o - worth 0 for symmetric
   // * - used
   // = - matches another '*' in symmetric mode (across diag. of 'x')
   //   ABCDEFGH
   //   --------
   // a|o*****ox      (considering pairs (X,x))
   // b|xo=====o
   // c|oxo=====
   // d|*oxo====
   // e|**oxo===
   // f|***oxo==
   // g|****oxo=
   // h|*****oxo
   // thus, in symmetric modes we take pairs:
   // [a,NNN(a)] to [a,P(a)]                                  (1)
   // [y,NNN(y)] to [y,a], where y is in [N(a),PPP(A)]        (2)
   //
   // Resulting nhd size is n(n-3)/2


   int x,y;

   // (1)
   for(x=3;x<size;x++) {
      double delta = getTwoChangeDelta(0,x);
      if(delta < bestDelta) {
         bestDelta=delta;
         best1=0;
         best2=x;
      }
   }

   // (2)
   for(y=1;y<size-2;y++) {

      for(x=getNext(getNext(getNext(y)));
          x != 1;  // same as getNext(0)
          x = getNext(x))
         {
            double delta = getTwoChangeDelta(y,x);
            if(delta < bestDelta) {
               bestDelta=delta;
               best1=y;
               best2=x;
         }
      }
   }

   if(bestDelta < 0) {
      TSPSolution *newSolution = new TSPSolution(this);
      newSolution->applyTwoChange(best1,best2,bestDelta);
      return newSolution;
   }

   else
      return this;
}

TSPSolution *TSPSolution::getLocalOptimum() {
   TSPSolution *last,*curr;
   last = this;
   while(1) {
      curr = last->getBestTwoChange();
      if(curr==last) {
         break;
      }
      if(last!=this)
         delete last;
      last = curr;
   }

   return curr;
}

TSPSolution *TSPSolution::getNhdOptimum() {
   TSPNhd *nhd = getNhd();
   TSPSolution *best = this;

   nhd->calcFamilyStatus();
   for(int i=0;i<nhd->getSize();i++) {
      if(nhd->getFamilyStatus(i)==NONFAMILY) {
         TSPSolution *nbr = nhd->getSolution(i);
         TSPSolution *opt = nbr->getLocalOptimum();
         if(opt->getValue() - best->getValue() < 0.00001) {
            best = opt;
         }
         else {
            delete opt;
         }
      }
   }

   delete nhd;

   return best;
}

TSPSolution *TSPSolution::getRandomTwoChange() {
   int c1 = rand() % size;
   int c2 = rand() % (size-1);

   // Eliminate illegal combination.
   if(c2==getNext(c1))
      c2 = size-1;

   TSPSolution *newSolution = new TSPSolution(this);
   newSolution->applyTwoChange(c1,c2);
   return newSolution;
}

void TSPSolution::citySwap(int c1,int c2) {
   int tmp = cities[c1];
   cities[c1] = cities[c2];
   cities[c2] = tmp;
   calcCost();
}

TSPSolution *TSPSolution::getRandomCitySwap() {
   int c1 = rand()%size;
   int c2 = rand()%(size-1);
   if(c2==c1)
      c2 = size-1;

   TSPSolution *rtn = new TSPSolution(this);
   rtn->citySwap(c1,c2);
   return rtn;
}

void calcFamilyStatus();

TSPNhd *TSPSolution::getNhd() {
   // The diagram below shows which combinations are needed for 2-ch
   // neigbourhood checking
   // x - illegal
   // o - worth 0 for symmetric
   // * - used
   // = - matches another '*' in symmetric mode (across diag. of 'x')
   //   ABCDEFGH
   //   --------
   // a|o*****ox      (considering pairs (X,x))
   // b|xo=====o
   // c|oxo=====
   // d|*oxo====
   // e|**oxo===
   // f|***oxo==
   // g|****oxo=
   // h|*****oxo
   // thus, in symmetric modes we take pairs:
   // [a,NNN(a)] to [a,P(a)]                                  (1)
   // [y,NNN(y)] to [y,a], where y is in [N(a),PPP(A)]        (2)
   //
   // Resulting nhd size is n(n-3)/2


   int x,y;

   TSPNhd *nhd = new TSPNhd(this,size*(size-3)/2);
   int index = 0;

   // (1)
   for(x=3;x<size;x++) {
      TSPSolution *nbr = new TSPSolution(this);
      nbr->applyTwoChange(0,x);
      nhd->setSolution(nbr,index);
      index++;
   }

   // (2)
   for(y=1;y<size-2;y++) {
      for(x=getNext(getNext(getNext(y)));
          x != 1;  // same as getNext(0)
          x = getNext(x))
      {
         TSPSolution *nbr = new TSPSolution(this);
         nbr->applyTwoChange(y,x);
         nhd->setSolution(nbr,index);
         index++;
      }
   }

   return nhd;
}

int TSPSolution::equals(TSPSolution *cmp) {
   if(fabs(cost - cmp->cost) > 0.0001)
      return 0;
   int start1,start2;
   for(int i=0;i<size;i++) {
      if(cities[i]==0)
         start1=i;
      if(cmp->cities[i]==0)
         start2=i;
   }

   for(int i = 0;i < size;i ++) {
      if(cities[(start1+i)%size] != cmp->cities[(start2+i)%size])
         return 0;
   }

   return 1;
}

/**
 * Class TSP Nhd
 */
void TSPNhd::calcFamilyStatus() {
   family=0;
   for(int i=0;i<size;i++) {
      TSPSolution *opt = nhd[i]->getBestTwoChange();
      if(opt->equals(source)) {
         familyStatus[i] = IMMFAMILY;
         family++;
      } else {
         familyStatus[i] = NONFAMILY;
      }
      if(opt!=nhd[i]) {
         delete opt;
      }
   }

   familyMembers = new int[family];
   int index = 0;
   for(int i=0;i<size;i++) {
      if(familyStatus[i] == IMMFAMILY) {
         familyMembers[index] = i;
         index++;
      }
   }
}


int tspSolutionSort(const void *a,const void *b) {
   double ac = ((TSPSolution *)a)->getValue();
   double bc = ((TSPSolution *)b)->getValue();

   if (ac - bc < -0.0001) {
      return -1;
   }
   else if (ac - bc > 0.0001) {
      return 1;
   }
   return 0;
}


