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
void init_parallel_edges(){
  int i,nbr_i1,nbr_i2,j;
  int e1,e2,e3;//edges
  pedges_i=(int*)malloc(nedges*4*sizeof(int));
  pedges_j=(int*)malloc(nedges*4*sizeof(int));
  int ctr;
  double tol=1e-4;
  double xdiff,ydiff;
  double dist;
  tot_plaqs=0;
  for(i=0; i<n_vtx; i+=1){
    for(e1=fnbr[i]; e1<fnbr[i+1]; e1++){
      ctr=0;
      nbr_i1=vptr[e1];
      for(e2=fnbr[i]; e2<fnbr[i+1]; e2++){
        nbr_i2=vptr[e2];
        if(nbr_i2!=nbr_i1){
          for(e3=fnbr[nbr_i2]; e3<fnbr[nbr_i2+1]; e3+=1){
            j=vptr[e3];
            dist=(xpos[j]-xpos[nbr_i1])*(xpos[j]-xpos[nbr_i1])+(ypos[j]-ypos[nbr_i1])*(ypos[j]-ypos[nbr_i1]);
            if((j!=i)&&(fabs(dist-1.0)<tol)){
              pedges_i[2*e1+ctr]=nbr_i2;;
              pedges_j[2*e1+ctr]=j;;
              ctr+=1;
              assert(ctr<=2);
            }
          }

        }
      }
      if(ctr==1){
        pedges_i[2*e1+1]=-1;
        pedges_j[2*e1+1]=-1;

      }
      tot_plaqs+=ctr;
    }
  }
  tot_plaqs=tot_plaqs/8;

}

void make_orien(){
  int edge=0;
  int site,site2;
  double phi;
  int i,j;
  for(i=0; i<n_vtx; i++){
      for(j=fnbr[i];j<fnbr[i+1];j++){
	site=i;
	site2=vptr[j];
	phi=atan2(ypos[site2]-ypos[site],xpos[site2]-xpos[site])/(2*M_PI);
	orien[edge]=round(phi*8);

	edge++;
      }
  }

}


