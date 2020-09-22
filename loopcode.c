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
  }
  int edge=0;
  int edge2=0;
  edge=0;
  for(i=0; i<n_vtx; i++){
    if(charge[i]==0){
      for(j=fnbr[i];j<fnbr[i+1];j++){
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
  int i;
  //printf("el:%d\t",entry_leg);
  //for(i=0; i<dim;i++)
  //  printf("%f\t",ptab[entry_leg+i*dim]);
  //printf("\n");
  while(z>ptab[entry_leg+exit_leg*dim])
    exit_leg++;
  assert(exit_leg<dim);
  assert(entry_leg<dim);
  return exit_leg;

}
void cancel_loop(int loopstack_top, int start_site){
  int site=start_site;
  int site3,site2,temp;
  loopstack_top--;
  while(loopstack_top>=0){
    site3=match[site];
    temp=match[site3];
    match[site3]=site;
    site2=loopstack[loopstack_top--];
    match[temp]=site2;
    //printf("site2=%d\t",site2);
    if(match[site2]==-1){
      match[site2]=temp;
      break;
    }
    site=temp;
  }
  assert(loopstack_top==-1);

}
void cancel_sloop(int loopstack_top, int start_site){
  int site=start_site;
  int site3,site2,temp;
  loopstack_top--;
  while(loopstack_top>=0){
    site3=match[site];
    temp=match[site3];
    match[site3]=site;
    //printf("s=%d s3=%d temp=%d\n",site,site3,temp);
    site2=loopstack[loopstack_top--];
    match[temp]=site2;
    //printf("site2=%d\n",site2);
    if(match[site2]==-1){
      match[site2]=temp;
      break;
    }
    site=temp;
  }
  assert(loopstack_top==-1);

}
void trim_loop(int loopstack_top, int start_site){
  int site=start_site;
  int site3,site2,temp;
  //printf("startsite=%d\n",site);
  while(loopstack_top>=0){
    site2=loopstack[loopstack_top--];
    //printf("site=%d top=%d site2=%d ms2%d\n",site,loopstack_top+1,site2,match[site2]);
    match[site]=site2;
    if(match[site2]==-1){
      match[site2]=site;
      //printf("here now bitch \n");
      break;
    }
    temp=match[site2];
    match[site2]=site;
    site=temp;
  }
  assert(loopstack_top==-1);

}
//long loop
/* Nomencalture for variables   
 * Each step of the loop starts with 'site' which was matched
 * to 'site2', but we now choose 'site3' to match it. match[site3]
 * now becomes the new 'site' for the next iteration of the loop.*/

int basic_loop(){
  int site;int site3;
  int n_nbrs;
  int startsite2;
  int temp;
  int site2=(int)(genrand64_real2()*n_vtx);
  int loopflag=0;
  int config;
  int j;
  int s1,s2,edge;
  int nparallel;
  double *ptab_l;//local prob table;
  int entry_leg,exit_leg;
  int temp_looplen=0;
  if(match[site2]!=-1){
    site=match[site2];
    assert(match[site2]!=-1);
    match[site2]=-1;

    startsite2=site2;
    n_nbrs=fnbr[site+1]-fnbr[site];
    int loopstack_top=0;
    //printf("loop start\n");
    while(1){
      loopstack[loopstack_top++]=site2;
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
      temp_looplen+=2;
      //printf("s2=%d s1=%d s3=%d\n",site2,site,site3);
      if(match[site3]==-1){
        match[site3]=site;
        match[site]=site3;
        break;
      }
      else{
        match[site]=site3;
        temp=match[site3];
        match[site3]=site;
        site=temp;
        site2=site3;
      }
      if(temp_looplen>(n_vtx*CUTOFF_FAC)){
        //printf("cancelled\n");
        //getchar();
        cancel_loop(loopstack_top,site);//undo changes
        return 0;
      }

    }
  }
  else
    return 0;//monomer start
  nloops++;
  looplen+=temp_looplen;
  return 1;

}
void check_matching(){
  int i;
  int unmatched=0;
  for(i=0; i<n_vtx; i++){
    if(match[i]!=-1){
      if(match[match[i]]!=i)
        printf("i %d match %d matchmatch%d\n",i,match[i],match[match[i]]);
      assert(match[match[i]]==i);
    }
    else
      unmatched+=1;
  }
  assert(unmatched==1);
}

int short_loop(){
  check_matching();
  int site;int site3;
  int n_nbrs;
  int startsite2;
  int temp;
  int site2=(int)(genrand64_real2()*n_vtx);
  int loopflag=0;
  int config;
  int j;
  int s1,s2,edge;
  int nparallel;
  //printf("nsloop=%d\n",nsloops);
  double *ptab_l;//local prob table;
  int entry_leg,exit_leg;
  int temp_looplen=0;
  if(match[site2]!=-1){
    site=match[site2];
    assert(match[site2]!=-1);
    match[site2]=-1;

    startsite2=site2;
    n_nbrs=fnbr[site+1]-fnbr[site];
    int loopstack_top=0;
    int backtrack=0;
    //printf("loop start\n");
    //printf("nsloops=%d\n",nsloopctr);
    loopstack[loopstack_top++]=site2;
    while(1){
      n_nbrs=fnbr[site+1]-fnbr[site];
      assert(n_nbrs>=2);
      config=0;
      entry_leg=-1;
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
      if(entry_leg==-1){
        printf("site2=%d\n",site2);
        printf("site=%d\n",site);
        for(edge=fnbr[site]; edge<fnbr[site+1]; edge++)
          printf("e %d nbr1 %d\n",edge,vptr[edge]);
        for(edge=fnbr[site2]; edge<fnbr[site2+1]; edge++)
          printf("e %d nbr2 %d\n",edge,vptr[edge]);
      }

      assert(entry_leg!=-1);
      ptab_l=(probtab+ptab_start[n_nbrs]+config*n_nbrs*n_nbrs);
      exit_leg=draw_exit_leg(entry_leg,ptab_l,n_nbrs);
      site3=vptr[fnbr[site]+exit_leg];
      //test for bounce processes
      temp_looplen+=2;
      //printf("s2=%d s1=%d s3=%d\n",site2,site,site3);
      //printf("vs2=%d vs1=%d vs3=%d\n",lvisits[site2],lvisits[site],lvisits[site3]);
      //getchar();
      backtrack=0;
      if(loopstack_top>0)
        if(site3==loopstack[loopstack_top-1]){
          //printf("lstop %d val %d\n",loopstack_top-1,loopstack[loopstack_top-1]);
          loopstack_top--;
          lvisits[site3]=-1;
          backtrack=1;
         // printf("back_track\n");
          //getchar();
        }

      if(match[site3]==-1){
        match[site3]=site;
        match[site]=site3;
        break;
      }
      else{
        match[site]=site3;
        temp=match[site3];
        match[site3]=site;
        site=temp;
        site2=site3;
      }
      //check and correct for self-intersections
      if(!backtrack){
        if(lvisits[site3]==nsloopctr){
          temp_looplen=0;
          //printf("self intersection\n");
          //printf("s3=%d vs3=%d\n",site3,lvisits[site3]);
          //getchar();
          loopstack_top--;
          //printf("top=%d\n",loopstack_top);
          while(loopstack[loopstack_top--]!=site3){
            temp_looplen+=2;
            //printf ( "s2 %d\n",loopstack[loopstack_top] );
            // getchar();
          }
          //printf("top=%d\n",loopstack_top);
          //printf("top=%d, lsval%d site=%d\n",loopstack_top,loopstack[loopstack_top],site);
          //getchar();
          temp_looplen+=2;
          assert(loopstack[loopstack_top+1]==site3);
          //printf("call\n");
          trim_loop(loopstack_top,site);
          check_matching();
          break;

        }
        else{
          lvisits[site3]=nsloopctr;
          loopstack[loopstack_top++]=site2;
          //printf("top=%d, lsval%d \n",loopstack_top-1,site2);

        }
      }
      if((!backtrack)&& (temp_looplen>(n_vtx*CUTOFF_FAC))){
        //loopstack_top--;
        //loopstack_top--;
        //printf("top=%d, lsval%d site=%d s3 %d\n",loopstack_top-1,loopstack[loopstack_top-1],site,site3);
        //int st;
        //for(st=0;st<loopstack_top;st++)
        //  printf("%d %d\n",st,loopstack[st]);
        //getchar();
        loopstack_top--;
        cancel_sloop(loopstack_top,site);//undo changes
        nsloopctr++;
        return 0;
      }

    }
  }
  else{
    nsloopctr++;
    return 0;//monomer start
  }
  nsloops++;
  nsloopctr++;
  slooplen+=temp_looplen;
  return 1;

}
int loop_nocutoff(){
  int site;int site3;
  int n_nbrs;
  int startsite2;
  int temp;
  int site2=(int)(genrand64_real2()*n_vtx);
  int loopflag=0;
  int config;
  int j;
  int s1,s2,edge;
  int nparallel;
  double *ptab_l;//local prob table;
  int entry_leg,exit_leg;
  int temp_looplen=0;
  if(match[site2]!=-1){
    site=match[site2];
    assert(match[site2]!=-1);
    match[site2]=-1;

    startsite2=site2;
    n_nbrs=fnbr[site+1]-fnbr[site];
    //printf("loop start\n");
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
      temp_looplen+=2;
      //printf("s2=%d s1=%d s3=%d\n",site2,site,site3);
      if(match[site3]==-1){
        match[site3]=site;
        match[site]=site3;
        break;
      }
      else{
        match[site]=site3;
        temp=match[site3];
        match[site3]=site;
        site=temp;
        site2=site3;
      }
    }
  }
  else
    return 0;//monomer start
  nloops++;
  looplen+=temp_looplen;
  return 1;

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
