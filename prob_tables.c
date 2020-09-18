#include "global.h"
#include "lp_lib.h"
/* In the detailed-balance matrix, the variables are numbered as follows
 * 0
 * 1 2
 * 3 4 5
 * 6 7 8 9*/
void lpp_wrapper(int dim,double *w,double *a){ //dim by dim detailed balance eq
  lprec *lp;
  int i,Ncol,*colno=NULL,j,ret=0;
  double *row = NULL;
  Ncol=dim*(dim+1)/2;
  colno = (int *) malloc(Ncol * sizeof(int));
  row = (double *) malloc(Ncol * sizeof(double));
  lp=make_lp(0,Ncol);
  int rowno,varno,startvar,var,incr,temp;
  int vctr=0;
  char vn[5];
  for(rowno=0; rowno<dim; rowno ++){
    for(temp=0;temp<=rowno;temp++){
      sprintf(vn,"a%d%d",rowno,temp);
      set_col_name(lp,vctr+1,vn);
      vctr++;
    }
  }

  set_add_rowmode(lp,TRUE);
  for(rowno=0; rowno<dim; rowno++){
    var=rowno*(rowno+1)/2;
    incr=rowno+1;
    for(varno=0; varno<dim; varno++){
      if(varno<rowno){
        colno[varno]=var+1;
        var++;
      }
      else{
        colno[varno]=var+1;
        var+=incr;
        incr++;
      }
      row[varno]=1;

    }
    add_constraintex(lp, dim, row, colno, EQ, w[rowno]);
  }

  set_add_rowmode(lp, FALSE); /*  rowmode should be turned off again when done building the model */
  var=0;
  for(varno=0; varno<dim; varno++){
    colno[varno]=var+1;
    var+=(varno+2);
    row[varno]=1;
  }
  set_obj_fnex(lp,dim,row,colno);
  //write_LP(lp,stdout);
  set_verbose(lp,SEVERE);
  ret=solve(lp);
  assert(ret==0);
  get_variables(lp,row);
  //printf("Objective value: %f\n", get_objective(lp));
  //for(j = 0; j < Ncol; j++)
  //    printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
  vctr=0;
  for(rowno=0; rowno<dim; rowno ++){
    for(temp=0;temp<=rowno;temp++){
      a[rowno+temp*dim]=row[vctr];
      assert(row[vctr]>=-1e-8);
      if(temp!=rowno)
        a[temp+rowno*dim]=row[vctr];
      vctr++;
    }
  }
  free(row);
  free(colno);
  delete_lp(lp);
  double tol;
  //for(i=0; i<dim; i++){
  //  for(j=0;j<dim; j++)
  //   printf("a%d%d=%f\t",i,j,a[i+j*dim]);
  //  printf("\n");
  //}
  //normalize
  for(i=0; i<dim; i++){
    for(j=0;j<dim; j++)
      a[i+j*dim]/=w[i];
  }
  double psum;
  for(i=0; i<dim; i++){
    psum=0;
    for(j=0;j<dim; j++)
      psum+=a[i+j*dim];
    assert(fabs(psum-1)<1e-4);
  }
  //make cumulative probability tables
  for(i=0; i<dim; i++){
    for(j=1;j<dim; j++)
      a[i+j*dim]+=a[i+(j-1)*dim];
    assert(fabs(a[i+dim*(dim-1)]-1)<1e-4);
  }

}


void construct_probtabs(){
  int ncoord=2;
  int nconfigs;
  int config;
  fugacity=1;
  int edge;
  int nflippable;
  double *w,*a;
  int temp;
  //counting loop
  int n_ptab_entries=0;;
  for(ncoord=2; ncoord<=8; ncoord++){
    ptab_start[ncoord]=n_ptab_entries;
    nconfigs=round(pow(3,ncoord));
    //for(config=0; config<nconfigs; config+=1){
    n_ptab_entries+=nconfigs*ncoord*ncoord;
  }
  probtab=(double*)malloc(n_ptab_entries*sizeof(double));
  printf("ptab entries:%d size:%f\n",n_ptab_entries,8.0*n_ptab_entries/(pow(1024,3.0)));
  int startpos=0;
  for(ncoord=2; ncoord<=8; ncoord++){
    nconfigs=round(pow(3,ncoord));
    //for(config=0; config<nconfigs; config+=1){
    for(config=0; config<nconfigs; config+=1){
      w=malloc(ncoord*sizeof(double));
      a=malloc(ncoord*ncoord*sizeof(double));
      temp=config;
      for(edge=0; edge<ncoord; edge++){
        nflippable=temp%3;
        temp=temp/3;
        w[edge]=exp(-rk_v*nflippable);

      }
      lpp_wrapper(ncoord,w,probtab+startpos);
      startpos+=ncoord*ncoord;
      free(w);
      free(a);
    }

  }
}
