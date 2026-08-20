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


void conststates1(){
  extern int checkGL2(std::string);
  extern int checkGL3(boost::dynamic_bitset, int);
  extern void storeConf(std::string);
  long long int NST,n;
  std::string str,gaugeConf; 
  long phyS;
  int check;
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
    // maximal 32 sites can be dealt with here; this is without using bitset
    // str = std::bitset<32> (n).to_string();
    // gaugeConf = str.substr(str.size() - LX);
    boost::dynamic_bitset<> gaugeConf(LX, n);
    
    //std::cout<<gaugeConf<<" "<<n<<std::endl;
    //check = checkGL2(gaugeConf);
    //if(check) storeConf(gaugeConf); 
    phyS += check;
    n++;
  }
  std::cout<<"States satisfying Gauss Law = "<<phyS<<std::endl;
  std::cout<<basisPhys.size()<<std::endl;
  for(size_t i=0; i<phyS; i++){
    for(size_t j=0; j<LX; j++){
    //   if(basisPhys[i][j]==true) printf("%d",1);
    //   if(basisPhys[i][j]==false) printf("%d",0);
    }
    printf("\n"); 
  }
}

// Gauss Law at site x: Gx = E_(i,i+1) - E_(i-1,i) - (n_x + sign)
// Imposing Gx = 0 ==> nx = [E(i,i+1) - E(i-1,i)] - sign
// this function passes a string to check Gauss Law
int checkGL2(std::string gC){
   int i, lf, lb;
   double LinkF, LinkB, OccF;
   for(i=0;i<LX;i++){
      // lf = link between (i,i+x) and lb = link between (i-x,i)
      lf = (int)(gC[i]-'0'); lb = (int)(gC[next[DIM-1][i]]-'0'); 
      // compute the fermion number by using Gauss Law    
      LinkF = lf - spin; LinkB = lb - spin; 
      OccF  = LinkF - LinkB - sign[i] ;
      //std::cout<<"site i="<<i<<"; lf="<<LinkF<<"; lb="<<LinkB<<"; nf="<<OccF<<std::endl;
      if(OccF != 0.0 && OccF != 1.0) return 0;
   }
   return 1;
}

// A similar function to check Gauss Law using boost
int checkGL3(boost::dynamic_bitset gC, int L){
}

void storeConf(std::string gC){
   std::vector<bool> conf(LX);
   int num=std::stoi(gC,nullptr,2);
   //conf = (bool)num;  
   std::cout<<gC<<" "<<num<<std::endl;
   //for(size_t i=0; i<LX; i++){
   //   conf[i] = (bool)std::stoi(gC[i]);
   //   std::cout<<i<<" "<<gC[i]<<" "<<conf[i];
   //}
   //conf = (bool)std::stoi(gC);
   std::cout<<std::endl;
   basisPhys.push_back(conf);
}


void conststates(){
  extern void addF(char *str);
  extern void addG(char *str);
  extern int checkGL(char *str1,char *str2);
  long long int NST,NF,NG,n,n1,n2;
  char fconf[LX],gconf[LX];
  long phyS;
  int flag;
  int i;

  // total Hilbert space; naive counting
  NF  = pow(2,LX);
  NG  = pow(2,LX);
  NST = NF*NG;
  std::cout<<"Total Hilbert space with "<<LX<<" sites and spin-1/2 links is="<<NST<<std::endl;

  //initialize fermionic and gauge configuration
  for(i=0;i<LX;i++){
     fconf[i]='0'; gconf[i]='0';
   }
   fconf[LX]='\0'; gconf[LX]='\0';
   // check print fermionic configuration
   //n=0;
   //while(n < NF){
   //    std::cout<<"Fermion config ="<<fconf<<std::endl;
   //    addF(fconf);
   //    n++;
   //}

   // in this notation, the fermion occupation at a site goes from 0 to 1
   // the gauge occuptation goes from 0 to repS = 2*S+1. The E value is
   // obtained as E = -spin+x.
   // For S=1/2: x=(0,1) ==> E=(-0.5, 0.5)
   // For S=1:   x=(0,1,2) ==> E=(-1,0,1)
   // For S=3/2; x=(0,1,2,3) ==> E=(-1.5,-0.5,0.5,1.5)
   // check print gauge configuration
   //std::cout<<gg<<" "<<gh<<std::endl;
   //n=0;
   //while(n < NG){
   //    std::cout<<"Gauge config ="<<gconf<<std::endl;
   //  addG(gconf);
   //    n++;
   //}
   //std::cout<<"Total config = "<<n<<std::endl;

   // actual state counting
   phyS=0; n1=0;
   while(n1 < NF){
     n2=0;
     while(n2 < NG){
       flag=0;
       flag=checkGL(fconf,gconf);
       if(flag == 1){
          phyS++;
          //std::cout<<"Fermion config ="<<fconf<<std::endl;
          //std::cout<<"Gauge   config ="<<gconf<<std::endl;
        }
       addG(gconf);
       n2++;
     }
     addF(fconf);
     n1++;
   }
   std::cout<<"Total number of physical states = "<<phyS<<std::endl;

}

void addF(char *str){
  int k,carry;
  carry=0; k=LX-1;
  if(str[k]=='0') { str[k]='1'; return; }
  else { str[k]='0'; carry=1; }
  while(k>0){
    k--;
    if((carry==1)&&(str[k]=='1')) str[k]='0';
    else if(str[k]<'1') { str[k]=str[k]+1; carry=0; return; }
  }
}

void addG(char *str){
  int k,carry;
  carry=0; k=LX-1;
  if(str[k]<gh) { str[k]=str[k]+1; return; }
  else { str[k]='0'; carry=1; }
  while(k>0){
    k--;
    if((carry==1)&&(str[k]==gh)) str[k]='0';
    else if(str[k]<gh) { str[k]=str[k]+1; carry=0; return; }
  }
}

int checkGL(char *fC,char *gC){
   int i,sign;
   int lf, lb, nf;
   double LinkF, LinkB, OccF;
   double Q;
   for(i=0;i<LX;i++){
      if(((i+1)%2)==0) sign=0;
      else sign=-1;
      lf = (int)(gC[i]-'0'); lb = (int)(gC[next[DIM-1][i]]-'0'); nf = (int)(fC[i]-'0');
      LinkF = lf - spin; LinkB = lb - spin; OccF = (double)(nf);
      Q  = LinkF - LinkB - ( OccF + sign );
      //std::cout<<"site i="<<i<<"; lf="<<LinkF<<"; lb="<<LinkB<<"; nf="<<OccF<<
      //"; Q="<<Q<<std::endl;
      if(Q != 0.0 ) return 0;
   }
   return 1;
}
