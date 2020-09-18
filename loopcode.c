#include "global.h"
#include "assert.h"

double genrand64_real2();

void measure(){
  int i,j,k,l;
  if(nmeasure%binsize==0){
    binno=0;
    for(i=0;i<n_vtx;i++)
      mdensity[i]=0;
    for(i=0;i<nedges;i++)
      ddensity[i]=0;
    for(i=0;i<5;i++)
      dorien[i]=0;
    n_flippable=0;
    n_flippable2=0;
  }
  for(i=0;i<n_vtx;i++){
    mdensity[i]+=(match[i]==-1);
    //for(j=0;j<n_vtx;j++){
    //  tmmcorr[i+j*n_vtx]+= ((match[i]==-1)&&(match[j]==-1));
    //}
  }
  int edge=0;
  int edge2=0;
  edge=0;
  for(i=0; i<n_vtx; i++){
    if(charge[i]==0){
      for(j=fnbr[i];j<fnbr[i+1];j++){
        //dorien[orien[edge]]+=1.0;
        ddensity[edge]+=(match[i]==vptr[j]);
        edge++;
      }

    }
  }
  int s1,s2,temp_nf;
  temp_nf=0;
  for(i=0; i<n_vtx; i++){
    for(edge=fnbr[i];edge<fnbr[i+1];edge++){
      if(match[i]==vptr[edge]){
        s1=pedges_i[2*edge];
        s2=pedges_j[2*edge];
        if(match[s1]==s2)
          temp_nf+=1;
        if(pedges_i[2*edge+1]!=-1){
          s1=pedges_i[2*edge+1];
          s2=pedges_j[2*edge+1];
          if(match[s1]==s2)
            temp_nf+=1;

        }
      }
    }
  }
  n_flippable+=(temp_nf/4);
  n_flippable2+=(temp_nf/4)*(temp_nf/4);

  if(nmeasure%binsize==(binsize-1)){
    FILE *binfp;
    binfp=fopen(binmfname,"a");
    fprintf(binfp,"%d %d ",binno,nmeasure);
    for(i=0;i<n_vtx;i++)
      fprintf(binfp,"%.16f ",mdensity[i]/((double)binsize));
    fprintf(binfp,"\n");
    fclose(binfp);

    binfp=fopen(bindfname,"a");
    fprintf(binfp,"%d %d ",binno,nmeasure);
    for(i=0;i<nedges;i++)
      fprintf(binfp,"%.16f ",ddensity[i]/((double)binsize));
    fprintf(binfp,"\n");
    fclose(binfp);

    binfp=fopen(binplaqfname,"a");
    fprintf(binfp,"%d %d ",binno,nmeasure);
    fprintf(binfp,"%.16f %.16f",n_flippable/((double)binsize),n_flippable2/((double)binsize));
    fprintf(binfp,"\n");
    fclose(binfp);
    //binfp=fopen(bindorienfname,"a");
    //fprintf(binfp,"%d %d ",binno,nmeasure);
    //for(i=0;i<5;i++)
    //  fprintf(binfp,"%.16f ",dorien[i]/((double)binsize));
    //fprintf(binfp,"\n");
    //fclose(binfp);
  }
  nmeasure++;
  binno++;

}
int draw_exit_leg(int entry_leg, double *ptab,int dim){
  double z=genrand64_real2();
  int exit_leg=0;
  while(z>ptab[entry_leg+exit_leg*dim])
    exit_leg++;
  return exit_leg;

}
void basic_loop(){

  int site;int site3;
  int n_nbrs;
  int startsite2;
  int temp;
  nloops++;
  int site2=(int)(genrand64_real2()*n_vtx);
  int loopflag=0;
  int config;
  int j;
  int s1,s2,edge;
  int nparallel;
  double *ptab_l;//local prob table;
  int entry_leg,exit_leg;

  if(match[site2]!=-1){
    site=match[site2];
    assert(match[site2]!=-1);
    match[site2]=-1;

    startsite2=site2;
    n_nbrs=fnbr[site+1]-fnbr[site];
    while(1){
      n_nbrs=fnbr[site+1]-fnbr[site];
      config=0;
      for(edge=fnbr[site]; edge<fnbr[site+1]; edge++){
        nparallel=0;
        if(vptr[edge]==site2)
          entry_leg=edge-fnbr[site];
        s1=pedges_i[2*edge];
        s2=pedges_j[2*edge];
        if(match[s1]==s2)
          nparallel+=1;
        if(pedges_i[2*edge+1]!=-1){
          s1=pedges_i[2*edge+1];
          s2=pedges_j[2*edge+1];
          if(match[s1]==s2)
            nparallel+=1;
        }
        config+=nparallel*round(pow(3,edge-fnbr[site]));
      }
      ptab_l=(probtab+ptab_start[n_nbrs]+config*n_nbrs*n_nbrs);
      exit_leg=draw_exit_leg(entry_leg,ptab_l,n_nbrs);
      site3=vptr[fnbr[site]+exit_leg];
      looplen+=2.0;
      if(match[site3]==-1){
        match[site3]=site;
        match[site]=site3;
        break;
      }
      else{
        match[site]=site3;
        temp=match[site3];
        n_nbrs=fnbr[temp+1]-fnbr[temp];
        match[site3]=site;
        site=temp;
        site2=site3;
      }

    }
  }

}
void monomer_move(){

  int site;int site3;
  site=(int)(genrand64_real2()*n_vtx);
  int n_nbrs;
  int startsite2;
  int temp;
  int loopflag=0;
  int site2=-1; //out of the lattice
  int rand;
  while(match[site=(int)(genrand64_real2()*n_vtx)]!=-1);
  while(1){
    //printf("site=%d , site2=%d \n",site,site2);
    n_nbrs=fnbr[site+1]-fnbr[site];
    rand=(int)(genrand64_real2()*(n_nbrs));
    site3 = vptr[fnbr[site]+rand];
    if(site3==site2){
      site3=-1;
      match[site]=-1;
      break;
    }


    match[site]=site3;
    temp=match[site3];
    match[site3]=site;
    site=temp;
    site2=site3;

  }

}
