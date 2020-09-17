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
	//if(match[i]==vptr[j]){
	  dorien[orien[edge]]+=1.0;
	//  edge2=0;
	  ddensity[edge]+=(match[i]==vptr[j]);
	//  for(k=0; k<n_vtx; k++){
	//    if(charge[k]==0){
	//      for(l=fnbr[k];l<fnbr[k+1];l++){
	//        if(match[k]==vptr[l])
	//          ddcorr[edge+edge2*nedges]+=1;
	//        edge2++;
	//      }
	//    }

	//  }
	//}
	edge++;
      }

    }
  }

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
    binfp=fopen(bindorienfname,"a");
    fprintf(binfp,"%d %d ",binno,nmeasure);
    for(i=0;i<5;i++)
      fprintf(binfp,"%.16f ",dorien[i]/((double)binsize));
    fprintf(binfp,"\n");
    fclose(binfp);
  }
  nmeasure++;
  binno++;

}
void basic_loop(){

  int site;int site3;
  int n_nbrs;
  int startsite2;
  int temp;
  nloops++;
  int site2=(int)(genrand64_real2()*n_vtx);
  int loopflag=0;
  if(match[site2]!=-1){
    site=match[site2];
    assert(match[site2]!=-1);
    match[site2]=-1;
    if((charge[site2]==1 && net_charge>0) || (charge[site2]==0 && net_charge<0))
      loopflag=0;
    else 
      loopflag=1;

    startsite2=site2;
    n_nbrs=fnbr[site+1]-fnbr[site];
    while(1){
      //printf("site=%d , site2=%d \n",site,site2);
      n_nbrs=fnbr[site+1]-fnbr[site];
      while( (site3 = vptr[fnbr[site]+(int)(genrand64_real2()*n_nbrs)] ) == site2);
      assert(site3!=site2);
      assert(site3!=site);
     // lcounts[startsite2+n_vtx*site]+=1.0;
      //if(loopflag==1){
      //  mmcorr[startsite2+n_vtx*site]+=1.0/(1.0*n_nbrs);
      //  //mmcorr[site+n_vtx*startsite2]+=1.0/(1.0*n_nbrs);
      //}
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
