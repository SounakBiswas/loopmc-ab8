#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "global.h"
int make_orien();
void ntiles(int,int*,int*);
void make_nbrs();
void find_charge();
int find_maxmatch();
void find_sitetypes();
void make_rtype();
void addvtx(double ax,double ay);
void init_parallel_edges();


int main(){
  FILE *fp=fopen("tiles_inf2.dat","r");
  double ax,ay,bx,by,cx,cy,dx,dy;
  int nt=16;
  int nr=16;
  ntiles(4,&nt,&nr);
  xpos=(double*)malloc((3*nt+4*nr)*sizeof(double)); 
  ypos=(double*)malloc((3*nt+4*nr)*sizeof(double));
  nedges=0;
  int first=1;
  int line=0;
  int typ;
  while(fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n", &typ,&ax, &ay, &bx, &by, &cx, &cy, &dx, &dy)!=EOF){
    if(first){
      lspace= sqrt((cx-bx)*(cx-bx)+(cy-by)*(cy-by));
      tol=lspace/100;
      first=0;
      printf("lspace=%f\n",lspace);
    }
    line++;
    addvtx(ax,ay);
    addvtx(bx,by);
    addvtx(cx,cy);
    if(typ==1)
      addvtx(dx,dy);
    //printf("line %d, ax %f \n",line,ax);
    //getchar();
    if(line%1000==0)
      printf("line=%d \n",line);
  }
  printf("nvertices,h=%d \n",n_vtx);
  fclose(fp);

  xpos=(double *)realloc(xpos, n_vtx*sizeof(double));
  ypos=(double *)realloc(ypos, n_vtx*sizeof(double));
  match=(int *)malloc(n_vtx*sizeof(double));
  charge=(int *)malloc(n_vtx*sizeof(double));
  make_nbrs();
  int i;
  for(i=0;i<n_vtx;i++){
    match[i]=-1;
    charge[i]=-1;
  }
  printf("edges: %d \n",nedges);
  find_charge();
  
  int mmatch=find_maxmatch();
  printf("max matching: %d \n",mmatch);
  find_sitetypes();
  make_rtype();
  fp=fopen("nodedata.dat","w");
  for(i=0; i<n_vtx; i++){
    fprintf(fp,"%d %f %f\n",i,xpos[i],ypos[i]);
  }
  fclose(fp);
  int j;
  fp=fopen("edgedata.dat","w");
  for(i=0; i<n_vtx; i++){
    for(j=fnbr[i];j<fnbr[i+1];j++)
     fprintf(fp,"%d %d\n",i,vptr[j]);
  }
  fclose(fp);
  fp=fopen("matched_edges.dat","w");
  for(i=0; i<n_vtx; i++){
    for(j=fnbr[i];j<fnbr[i+1];j++){
      if(match[i]==vptr[j])
     fprintf(fp,"%d %d\n",i,vptr[j]);
    }
  }
  fclose(fp);




  free(charge);
  free(xpos);
  free(ypos);
  free(match);
  free(vptr);
  free(fnbr);
  free(pedges_i);
  free(pedges_j);

   


}
