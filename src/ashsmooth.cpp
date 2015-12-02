#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List cxxtitable(SEXP sig){
  NumericVector signal=sig; 
  int n=(int) signal.size();
  int J=(int) log2((double)n);  

  NumericVector tempvec(2*n);
  NumericMatrix sumtable(J+1,n);
  NumericMatrix difftable(J+1,n);
  sumtable(0,_) = signal;
  difftable(0,_) = signal;

  for (int D=0; D<J; D++){
    int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
    for (int l=0; l<pD; l++){
      int a=l*nD+1;
      double d;
      for (int i=0; i<nD-1; i++){
        d=sumtable(D,a+i-1);
        tempvec(i)=d;
        tempvec(i+nD+1)=d;
      }
      //i=nD-1
      d=sumtable(D,a+nD-2);
      tempvec(nD-1)=d;
      tempvec(nD)=d;

      for (int i=0; i<nD; i++){
        sumtable(D+1,a+i-1)=tempvec(2*i)+tempvec(2*i+1);
        difftable(D+1,a+i-1)=tempvec(2*i)-tempvec(2*i+1);
      }
    }
  }
  return(List::create(Named("sumtable")=sumtable, Named("difftable")=difftable));
}



// [[Rcpp::export]]
NumericMatrix cxxtirtable(SEXP sig){
  NumericVector signal=sig; 
  int n=(int) signal.size();
  int J=(int) log2((double)n);

  NumericVector tempvec(2*n);
  NumericMatrix sumtable(J+1,n);
  NumericMatrix vrtable(J+1,n);
  sumtable(0,_) = signal;
  vrtable(0,_) = signal;

  for (int D=0; D<J; D++){
    int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
    for (int l=0; l<pD; l++){
      int a=l*nD+1;
      double d;
        for (int i=0; i<nD-1; i++){
        d=sumtable(D,a+i-1);
        tempvec(i)=d;
        tempvec(i+nD+1)=d;
      }
      //i=nD-1
      d=sumtable(D,a+nD-2);
      tempvec(nD-1)=d;
      tempvec(nD)=d;

      for (int i=0; i<nD; i++){
        sumtable(D+1,a+i-1)=tempvec(2*i)+tempvec(2*i+1);
        vrtable(D+1,a+i-1)=log((double)tempvec(2*i))-log((double)tempvec(2*i+1));
      }
    }
  }
  return(vrtable);
}


// [[Rcpp::export]]
NumericVector cxxreverse_gwave(SEXP estimate, SEXP pmat, SEXP qmat){
  NumericMatrix pp=pmat;
  NumericMatrix qq=qmat;
  NumericVector est1=estimate;
  int np=(int) pp.ncol();
  int J=(int) pp.nrow();
  NumericVector est(np,est1(0));
  //for(int D=J; D-->0;){
  for(int D=0; D<J; D++){
    //int nD=pow(2., (int) (J-D)), pD=pow(2., (int) D);
    int nD=pow(2., (int) (D+1)), pD=pow(2., (int) (J-1-D));
    int nDo2=nD/2;
    NumericVector tempvecl(nD), tempvecr(nD);
    for(int l=0; l<pD; l++){
      int a=l*nD+1; 
      double dep, deq, dp, dq;
      for (int i=0; i<nDo2; i++){
        dep=est(a+i-1)/2;
        //dp=pp(D,a+i-1);
        dp=pp(J-1-D,a+i-1);
        //dq=qq(D,a+i-1);
        dq=qq(J-1-D,a+i-1);
        tempvecl(2*i)=dep+dp;
        tempvecl(2*i+1)=dep+dq;
      }
      for (int i=nDo2; i<nD-1; i++){
        dep=est(a+i)/2;
        //dp=pp(D,a+i);
        dp=pp(J-1-D,a+i);
        deq=est(a+i-1)/2;
        //dq=qq(D,a+i-1);
        dq=qq(J-1-D,a+i-1);
        tempvecr(2*(i-nDo2))=deq+dq;
        tempvecr(2*(i-nDo2)+1)=dep+dp;
      }         
      //i=nD-1
      dep=est(a+nDo2-1)/2;
      //dp=pp(D,a+nDo2-1);
      dp=pp(J-1-D,a+nDo2-1);
      deq=est(a+nD-2)/2; 
      //dq=qq(D,a+nD-2);
      dq=qq(J-1-D,a+nD-2);
      tempvecr(nD-2)=deq+dq;
      tempvecr(nD-1)=dep+dp;
      for(int i=0; i<nD; i++)
        est(a+i-1)=0.5*(tempvecl(i)+tempvecr(i));
    }
  }
  return(est);
}



// [[Rcpp::export]]
NumericVector cxxreverse_gvwave(SEXP estimate, SEXP pmat, SEXP qmat){
  NumericMatrix pp=pmat;
  NumericMatrix qq=qmat;
  NumericVector est1=estimate;
  int np=(int) pp.ncol();
  int J=(int) pp.nrow();
  NumericVector est(np,est1(0));
  //for(int D=J; D-->0;){
  for(int D=0; D<J; D++){
    //int nD=pow(2., (int) (J-D)), pD=pow(2., (int) D);
    int nD=pow(2., (int) (D+1)), pD=pow(2., (int) (J-1-D));
    int nDo2=nD/2;
    NumericVector tempvecl(nD), tempvecr(nD);
    for(int l=0; l<pD; l++){
      int a=l*nD+1; 
      double dep, deq, dp, dq;
      for (int i=0; i<nDo2; i++){
        dep=est(a+i-1)/4;
        //dp=pp(D,a+i-1);
        dp=pp(J-1-D,a+i-1);
        //dq=qq(D,a+i-1);
        dq=qq(J-1-D,a+i-1);
        tempvecl(2*i)=dep+dp;
        tempvecl(2*i+1)=dep+dq;
      }
      for (int i=nDo2; i<nD-1; i++){
        dep=est(a+i)/4;
        //dp=pp(D,a+i);
        dp=pp(J-1-D,a+i);
        deq=est(a+i-1)/4;
        //dq=qq(D,a+i-1);
        dq=qq(J-1-D,a+i-1);
        tempvecr(2*(i-nDo2))=deq+dq;
        tempvecr(2*(i-nDo2)+1)=dep+dp;
      }         
      //i=nD-1
      dep=est(a+nDo2-1)/4;
      //dp=pp(D,a+nDo2-1);
      dp=pp(J-1-D,a+nDo2-1);
      deq=est(a+nD-2)/4; 
      //dq=qq(D,a+nD-2);
      dq=qq(J-1-D,a+nD-2);
      tempvecr(nD-2)=deq+dq;
      tempvecr(nD-1)=dep+dp;
      for(int i=0; i<nD; i++)
        est(a+i-1)=0.5*(tempvecl(i)+tempvecr(i));
    }
  }
  return(est);
}





// [[Rcpp::export]]
NumericVector cxxSParentTItable(SEXP sig){
  NumericVector signal=sig; 
  int n=(int) signal.size();
  int J=(int) log2((double)n);

  NumericVector parent(2*J*n);
  NumericMatrix TItable(J+1,n);
  TItable(0,_) = signal;
  for (int D=0; D<J; D++){
    int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
    for (int l=0; l<pD; l++){
      int a=l*nD+1, b=2*l*nD+1, c=2*D*n+b, d;
      for (int i=0; i<nD-1; i++){
        d=TItable(D,a+i-1);
        parent(c+i-1)=d;
        parent(c+i+nD)=d;
      }
      //i=nD-1
      d=TItable(D,a+nD-2);
      parent(c+nD-2)=d;
      parent(c+nD-1)=d;

	
      for (int i=0; i<nD; i++)
        TItable(D+1,a+i-1)=parent(c+2*i-1)+parent(c+2*i);
    }
  }
  return(parent);
}








// [[Rcpp::export]]
NumericVector cxxreverse_pwave(SEXP estimate, SEXP pmat, SEXP qmat){
  NumericMatrix pp=pmat;
  NumericMatrix qq=qmat;
  NumericVector est1=estimate;
  int np=(int) pp.ncol();
  int J=(int) pp.nrow();
  NumericVector est(np,est1(0));
  //for(int D=J; D-->0;){
  for(int D=0; D<J; D++){
    //int nD=pow(2., (int) (J-D)), pD=pow(2., (int) D);
    int nD=(int) pow(2., (int) (D+1)), pD=(int) pow(2., (int) (J-1-D));
    int nDo2=nD/2;
    NumericVector tempvecl(nD), tempvecr(nD);
    for(int l=0; l<pD; l++){
      int a=l*nD+1; 
      double dep, deq, dp, dq;
      for (int i=0; i<nDo2; i++){
	dep=est(a+i-1);
	//dp=pp(D,a+i-1);
	dp=pp(J-1-D,a+i-1);
	//dq=qq(D,a+i-1);
	dq=qq(J-1-D,a+i-1);
	tempvecl(2*i)=dep+dp;
	tempvecl(2*i+1)=dep+dq;
      }
      for (int i=nDo2; i<nD-1; i++){
	dep=est(a+i);
	//dp=pp(D,a+i);
	dp=pp(J-1-D,a+i);
	deq=est(a+i-1);
	//dq=qq(D,a+i-1);
	dq=qq(J-1-D,a+i-1);
	tempvecr(2*(i-nDo2))=deq+dq;
	tempvecr(2*(i-nDo2)+1)=dep+dp;
      }         
      //i=nD-1
      dep=est(a+nDo2-1);
      //dp=pp(D,a+nDo2-1);
      dp=pp(J-1-D,a+nDo2-1);
      deq=est(a+nD-2); 
      //dq=qq(D,a+nD-2);
      dq=qq(J-1-D,a+nD-2);
      tempvecr(nD-2)=deq+dq;
      tempvecr(nD-1)=dep+dp;
      for(int i=0; i<nD; i++)
               est(a+i-1)=0.5*(tempvecl(i)+tempvecr(i));
    }
  }
  return(est);
}






