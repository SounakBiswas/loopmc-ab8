#include "global.h"
int *visited;
int *queue;
int *parent;
int qhead,qtail,qlen;
void find_charge(){
  visited=(int *)malloc(n_vtx*sizeof(int));
  int *stack=(int *)malloc(n_vtx*sizeof(int));
  int site,nbrsite;
  int i,nid;
  for(i=0; i<n_vtx; i++){
    visited[i]=-1;
  }

  //find bipartite charge
  int stop=-1;
  stack[++stop]=0;
  charge[0]=0;
  while(stop!=-1){
    site=stack[stop--];
    if(visited[site]==-1){
      for(nid=fnbr[site]; nid<fnbr[site+1]; nid++){
	nbrsite=vptr[nid];
	if(charge[nbrsite]==-1){
	  charge[nbrsite]= (! charge[site]);
	  //printf("%d %d %d %d\n",site,nbrsite,charge[site],charge[nbrsite]);
	  stack[++stop]=nbrsite;
	}
      }
    }
  }

  int sum=0;
  for(i=0;i<n_vtx;i++){
    sum+= (2*charge[i]-1);
    //printf("%d %d\n",i,charge[i]);
  }
  printf("net charge=%d\n",sum);
  free(stack);
  free(visited);
}
void enqueue(int k){ 
  queue[qtail]=k;
  qtail=(qtail!=0)?qtail-1:n_vtx-1;
  qlen++;
}
void dequeue(){
  qhead=(qhead!=0)?qhead-1:n_vtx-1;
  qlen--;
}
void prune(int init_head){
  int i;
  if(qhead<init_head)
    for(i=qhead+1;i<init_head;i++)
      visited[match[queue[i]]]=n_vtx;
  else{
    for(i=qhead+1;i<n_vtx;i++)
      visited[match[queue[i]]]=n_vtx;
    for(i=0;i<init_head;i++)
      visited[match[queue[i]]]=n_vtx;
  }
}
void augment(int vtx){
  int temp;
  int vtxp;
  while(1){
    vtxp=parent[vtx];
    temp=match[vtxp];
    match[vtxp]=vtx;
    match[vtx]=vtxp;
    vtx=temp;
    if(temp==-1)
      break;
  }
  qhead=qtail=n_vtx-1;
  qlen=0;
}

int bfs_phase(int vtx,int n_phase){
  enqueue(vtx);
  int nid,nbrsite;
  while(qlen!=0){
    vtx=queue[qhead];
    dequeue();
    for(nid=fnbr[vtx]; nid<fnbr[vtx+1]; nid++){
      nbrsite=vptr[nid];
      if(visited[nbrsite]<n_phase){
	visited[nbrsite]=n_phase;
	parent[nbrsite]=vtx;
	if(match[nbrsite]==-1){
	  augment(nbrsite);
	  return 1;
	}
	else
	  enqueue(match[nbrsite]);
      }
    }

  }
  return 0;  //no augmenting path
}
int find_maxmatch(){
  visited=(int *)malloc(n_vtx*sizeof(int));
  queue=(int *)malloc(n_vtx*sizeof(int));
  parent=(int *)malloc(n_vtx*sizeof(int));
  int nphase;
  qhead=qtail=n_vtx-1;
  qlen=nphase=0;
  int vtx;
  int matching=0;
  int init_head;
  int i;
  for(i=0;i<n_vtx;i++)
    visited[i]=-1;
  for(vtx=0;vtx<n_vtx;vtx++){
    if(charge[vtx]==0 && match[vtx]==-1){
      init_head=qhead;
      if(bfs_phase(vtx,nphase)) //start bfs matching;
      matching++;
      else
	prune(init_head);
    }
    nphase+=1;
  }
  free(visited);
  free(parent);
  free(queue);
  return matching;
}
void mark_sites(int vtx,int nphase){
  enqueue(vtx);
  int nid,nbrsite;
  while(qlen!=0){
    vtx=queue[qhead];
    dequeue();
    sitetype[vtx]=1;
    for(nid=fnbr[vtx]; nid<fnbr[vtx+1]; nid++){
      nbrsite=vptr[nid];
      if(visited[nbrsite]==-1){
	visited[nbrsite]=nphase;
	sitetype[nbrsite]=2;
	if(match[nbrsite]!=-1)
	  enqueue(match[nbrsite]);
      }
    }

  }
}

void find_sitetypes(){
  int i,nphase;
  visited=(int *)malloc(n_vtx*sizeof(int));
  queue=(int *)malloc(n_vtx*sizeof(int));
  sitetype=(int*)calloc(n_vtx,sizeof(int));
  qhead=qtail=n_vtx-1;
  qlen=nphase=0;
  for(i=0;i<n_vtx;i++){
    visited[i]=-1;
    sitetype[i]=0;
  }
  for(i=0;i<n_vtx;i ++){
    if(match[i]==-1){
      mark_sites(i,nphase);
    }
    nphase++;

  }
  free(visited);
  free(queue);

}
void make_rtype(){

  int i;
  int *burn=(int*) malloc(n_vtx*sizeof(int));
  int *pocket=(int*) malloc(n_vtx*sizeof(int));
  for(i=0; i<n_vtx; i++)
    burn[i]=-1;
  int cidx=0;
  int pck_ctr,nbr;
  int cctr;
  int ctype,actype,nid;
  FILE *gp=fopen("csize_list.dat","w");
  int ch;
  for(i=0;  i<n_vtx; i++){
    ch=0;
    
    if(burn[i]==-1){
      pck_ctr=-1;
      cctr=0;
      int site,site0;
      site0=i;
      pocket[++pck_ctr]=site0;
      burn[site0]=cidx;
      ch+=2*charge[site0]-1;
      ctype=sitetype[site0];
      if(ctype==0){
       	while(pck_ctr!=-1){
       	  site=pocket[pck_ctr--];
       	  cctr++;
       	  for(nid=fnbr[site]; nid<fnbr[site+1]; nid++) {
       	    nbr=vptr[nid];
       	    if(burn[nbr]==-1 && sitetype[nbr]==0){
       	      pocket[++pck_ctr]=nbr;
       	      burn[nbr]=cidx;
       	    }
       	  }
       	}
       ;

      }
      else{
	while(pck_ctr!=-1){
	  site=pocket[pck_ctr--];
	  cctr++;
	  ctype=sitetype[site];
	  actype=1+!(ctype-1);
	  for(nid=fnbr[site]; nid<fnbr[site+1]; nid++) {
	    nbr=vptr[nid];
	    if(burn[nbr]==-1 && sitetype[nbr]==actype){
	      pocket[++pck_ctr]=nbr;
	      burn[nbr]=cidx;
              ch+=2*charge[nbr]-1;
	    }

	  }
	}
      }
      if(cctr!=1)
	fprintf(gp,"cidx=%d, cctr=%d ctype=%d charge=%d\n",cidx,cctr,ctype,ch);
      cidx++;
    }
  }
  fclose(gp);

  printf("number of monomer regions: %d \n",cidx);
  //for(i=0;i<n_vtx;i++)
  //  printf("%d %d\n",burn[i],sitetype[i]);
  FILE **fp;
  fp=(FILE **)malloc(cidx*sizeof(FILE*));
  char fname[200];
  for(i=0; i<cidx; i++){
    sprintf(fname,"./regions/clust%d.dat",i);
    fp[i]=fopen(fname,"w");
  }
  for(i=0; i<n_vtx; i++){
    fprintf(fp[burn[i]],"%d %.16f %.16f\n", i,xpos[i],ypos[i]);
  }
  for(i=0; i<cidx; i++){
    fclose(fp[i]);

  }
  free(fp);
  printf("nvtx=%d\n",n_vtx);
  gp=fopen("./clust0p.dat","w");
  for(i=0; i<n_vtx; i++){
    if((burn[i]==0)&&(fnbr[i+1]-fnbr[i])!=8)
    fprintf(gp,"%d %.16f %.16f\n", i,xpos[i],ypos[i]);
  }
  fclose(gp);


  free(sitetype);
  free(burn);
  free(pocket);

}
