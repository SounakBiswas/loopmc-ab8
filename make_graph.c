#include "global.h"
int *nbr;
int *nnbr;
void addvtx(double ax, double ay){
  int i;
  int flag=0;
  double dist;

  for(i=0;i< n_vtx;i++){
    dist=sqrt((ax-xpos[i])*(ax-xpos[i]) + (ay-ypos[i])*(ay-ypos[i]) );
    if(dist<tol){
      flag=1;
      break;
    }
  }
  if(flag==0){
    xpos[n_vtx]=ax;
    ypos[n_vtx]=ay;
    n_vtx+=1;
  }


}
void make_nbrs(){
  int i,j;
  double dist;
  int *nbr=(int *)malloc(8*n_vtx*sizeof(int));
  int *nnbr=(int *)calloc(n_vtx,sizeof(int));
  for(i=0;i<n_vtx;i++){
    for(j=0;j<i;j++){
      dist=sqrt((xpos[i]-xpos[j])*(xpos[i]-xpos[j]) +(ypos[i]-ypos[j])*(ypos[i]-ypos[j]));
      if(fabs( dist-lspace)<(lspace/100) ){
        //printf("%d %d \n",i,j);
        nbr[8*i+nnbr[i]]=j;
        nbr[8*j+nnbr[j]]=i;
        nnbr[i]++;
        nnbr[j]++;
        nedges+=1;

      }

    }
  }
  printf("nedges %d \n",nedges);
  vptr=(int*)malloc(2*nedges*sizeof(int));
  fnbr=(int*)malloc((n_vtx+1)*sizeof(int));
  int nidx=0;
  fnbr[0]=nidx;
  for(i=0;i<n_vtx;i++){
    //printf("%d %d \n",i,nnbr[i]);
    for(j=0; j<nnbr[i]; j++){
      vptr[nidx]=nbr[8*i+j];
      nidx++;
    }
    fnbr[i+1]=nidx;

  }
  free(nbr);
  free(nnbr);
}
void ntiles(int ninfl, int *nt, int *nr){
  int tempr,tempt;
  while(ninfl>0){
    tempt=4*nr[0]+3*nt[0];
    tempr=2*nt[0]+3*nr[0];
    nr[0]=tempr;
    nt[0]=tempt;
    ninfl--;
  }

}


