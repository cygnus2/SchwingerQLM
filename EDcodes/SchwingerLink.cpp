#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

/* according to the rules of cpp, the variables are declared here
 * and also in the header file as extern such that they are avl to
 * the other functions.
 */
int *next[2*DIM+1];
int *sign;
int LX,VOL;
double fMass, g2, V;
int NTOT,NH;
std::vector<int> basisPhys;
//std::vector<fermiN> nfBasis;
fermiN nfBasis;
int CHKDIAG;

int main(){
  FILE *fptr;
  char string[50];
  int i,d,p,q;
  int x,y;

  extern void initneighbor(void);
  extern void conststates(void);
  extern void FermiNo_decompose(void);
  extern void printbasis(void);
  extern int** allocateint2d(int, int);
  extern void deallocateint2d(int**,int,int);


  fptr = fopen("QUEUE","r");
  if(fptr == NULL){
      printf("could not open QUEUE FILE to open\n");
      exit(1);
  }
  fscanf(fptr,"%s %d\n",string,&LX);
  fscanf(fptr,"%s %lf\n",string,&fMass);
  fscanf(fptr,"%s %lf\n",string,&g2);
  fscanf(fptr,"%s %lf\n",string,&V);
  fclose(fptr);
  VOL = LX;

  // decide whether to check the results of the diagonalization
  CHKDIAG=0;

  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++){
    next[i] = (int *)malloc(VOL*sizeof(int));
  }
  sign = (int *)malloc(LX*sizeof(int));
  initneighbor();

  /* build basis states satisfying Gauss' Law */
  conststates();

  /* decomposing states using fermion number */
  FermiNo_decompose();

  /* construct the Hamiltonian */ 
  constH();

  /* breakup into translation sectors */
  //trans_decompose(sector);

  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); }
  free(sign);

  return 0;
}
