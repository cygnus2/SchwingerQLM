#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#include <bits/stdc++.h>
#include<string>
#include <boost/dynamic_bitset.hpp>
#include "define.h"
// local hilbert space of the links 2S+1
int repS=2;
char gg = '2';   // gg is the char of repS
char gh = '1';   // gh is the char of repS-1
double spin=(repS-1.0)/2.0;
long long int phyS;

void conststates(){
  extern int checkGL(int);
  long long int NST,n;
  bool check;
  // count total Hilbert space from Gauge fields
  NST  = pow(2,LX);
  std::cout<<"Total Hilbert space with "<<LX<<" spin-1/2 links is="<<NST<<std::endl;
  
  // assign the value to sign
  for(size_t i=0;i<LX;i++){
      if(((i+1)%2)==0) sign[i]= 0;
      else             sign[i]=-1;
  }
  
  n=0; phyS=0;
  while(n<NST){
    check = checkGL(n);
    phyS += (int)(check);
    if(check==true) basisPhys.push_back(n);
    n++;
  }
  std::cout<<"States satisfying Gauss Law = "<<phyS<<std::endl;
}

// Gauss Law at site x: Gx = E_(i,i+1) - E_(i-1,i) - (n_x + sign)
// Imposing Gx = 0 ==> nx = [E(i,i+1) - E(i-1,i)] - sign
// this function passes a string to check Gauss Law
int checkGL(int val){
   int i, lf, lb;
   double LinkF, LinkB, OccF;
   boost::dynamic_bitset<> gC(LX, val);
   
   for(i=0;i<LX;i++){
      // lf = link between (i,i+x) and lb = link between (i-x,i)
      lf = (int)(gC[i]); lb = (int)(gC[next[DIM-1][i]]); 
      // compute the fermion number by using Gauss Law    
      LinkF = lf - spin; LinkB = lb - spin; 
      OccF  = LinkF - LinkB - sign[i] ;
      //std::cout<<"site i="<<i<<"; lf="<<LinkF<<"; lb="<<LinkB<<"; nf="<<OccF<<std::endl;
      if(OccF != 0.0 && OccF != 1.0) return false;
   }
   return true;
}



