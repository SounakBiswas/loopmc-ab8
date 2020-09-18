#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "global.h"
int nedges;
int ntriangles(int,int,int);
void make_nbrs();
void make_orien();
int find_charge();
int find_maxmatch();
void find_sitetypes();
void make_rtype();
void addvtx(double ax,double ay);
void basic_loop();
void measure();
void monomer_move();
void init_genrand64(unsigned long long);
void init_parallel_edges();
void construct_probtabs();


int main(){
  n_vtx=769;
  rk_v=RK_V;
  temperature=TEMP;
  int cno=0;
  char fname[200];
  sprintf(fname,"./regions/clust%d.dat",cno);
  FILE *fp=fopen(fname,"r");
  xpos=(double*)malloc(n_vtx*sizeof(double)); 
  ypos=(double*)malloc(n_vtx*sizeof(double));
  int first=1;
  int line=0;
  int temp;
  lspace=1.0;
  while(fscanf(fp,"%d %lf %lf\n", &temp, &xpos[line], &ypos[line])!=EOF){
    line++;
  }
  printf("nvertices=%d \n",line);
  fclose(fp);

  match=(int *)malloc(n_vtx*sizeof(double));
  charge=(int *)malloc(n_vtx*sizeof(double));
  make_nbrs();
  int i;
  for(i=0;i<n_vtx;i++){
    match[i]=-1;
    charge[i]=-1;
  }
  printf("edges: %d \n",nedges);
  orien=(int *)malloc(nedges*sizeof(double));
  net_charge=find_charge();
  init_parallel_edges();
  //make_orien();
  construct_probtabs();
  
  int mmatch=find_maxmatch();
  printf("max matching: %d \n",mmatch);
  mdensity=(double *)malloc(n_vtx*sizeof(double));
  ddensity=(double *)malloc(nedges*sizeof(double));
  binno=0;
  nloops=0;
  binsize=100;
  corr1=corr2=0;
  looplen=0;
  sprintf(binmfname,"./outfiles/binfile_mon_v%.2f.dat",rk_v);
  sprintf(bindfname,"./outfiles/binfile_dim_v%.2f.dat",rk_v);
  sprintf(binplaqfname,"./outfiles/binfile_plaq_v%.2f.dat",rk_v);
  sprintf(bincorrfname,"./outfiles/bincorrfile_v%.2f.dat",rk_v);
  sprintf(bindorienfname,"./outfiles/bindorienfile_v%.2f.dat",rk_v);
  init_genrand64(4);
  //printf("charge10 %d charge11 %d netcharge %d",charge[10],charge[11],net_charge);
  //getchar();
  int j;
//  mmcorr=(double *)calloc(n_vtx*n_vtx,sizeof(double));
//  tmmcorr=(double *)calloc(n_vtx*n_vtx,sizeof(double));
//  ddcorr=(double *)calloc(nedges*nedges,sizeof(double));
  //ccorr_re=(double *)calloc(n_vtx*n_vtx,sizeof(double));
  //ccorr_im=(double *)calloc(n_vtx*n_vtx,sizeof(double));
//  lcounts=(double *)calloc(n_vtx*n_vtx,sizeof(double));
  for(i=0; i<10000; i++){
      basic_loop();
  }
  double avg_looplen=looplen/(1.0*nloops);
  int loops_in_mcstep=(1.0*n_vtx/(avg_looplen));
  loops_in_mcstep=(loops_in_mcstep>2)?loops_in_mcstep:2;
  printf("avg loop len=%f, loops=%d \n",avg_looplen,loops_in_mcstep);


  for(i=0; i<100000; i++){
    for(j=0; j<loops_in_mcstep; j++)
      basic_loop();
  //  //for(j=0; j<20; j++)
  //  //  monomer_move();
    measure();
    if(i%1000==0)
      printf("%d \n",i);
  }
  free(probtab);


  sprintf(fname,"nodedata%d.dat",cno);
  fp=fopen(fname,"w");
  for(i=0; i<n_vtx; i++){
    fprintf(fp,"%d %f %f\n",i,xpos[i],ypos[i]);
  }
  fclose(fp);
  printf("nedges=%d \n",nedges);
  sprintf(fname,"edgedata%d.dat",cno);
  fp=fopen(fname,"w");
  for(i=0; i<n_vtx; i++){
    if(charge[i]==0)
    for(j=fnbr[i];j<fnbr[i+1];j++)
      fprintf(fp,"%d %d\n",i,vptr[j]);
  }
  fclose(fp);
  //printf("%.16f %.16f %.16f\n",n_vtx*corr1/(1.0*nloops),n_vtx*corr2/(1.0*nloops),n_vtx*looplen/(1.0*nloops));
  //write_mmcors(cno);
  //write_tmmcors(cno);
  //write_ddcors(cno);



  free(xpos);
  free(match);
  free(ypos);
  free(charge);
  free(vptr);
  free(fnbr);
  free(mdensity);
//  free(mmcorr);
//  free(tmmcorr);
//  free(ddcorr);
  free(lcounts);
  free(ddensity);
  free(orien);
  return 1;



}
