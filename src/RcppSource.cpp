#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector PARAMETERS(int n,int d,int iteration){
  Function sampleC("sample.int");
  Function quantileC("quantile");
  double meansh=0;
  double meansl=0;
  for(int itr=0;itr<iteration;itr++){
    NumericVector distvec(n*(n-1)/2);
    for(int k=0;k<d;k++){
      int l=0;
      IntegerVector ind=sampleC(n);
      for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){
          int dif=ind[i]-ind[j];
          distvec[l]=distvec[l]+dif*dif;
          l++;
        }
      }
    }
    NumericVector sh=quantileC(distvec,_["probs"]=0.5);
    meansh=meansh+sqrt(sh[0]);
    NumericVector sl=quantileC(distvec,_["probs"]=0.05);
    meansl=meansl+sqrt(sl[0]);
  }
  meansh=meansh/(iteration+0.0)/n;
  meansl=meansl/(iteration+0.0)/n;
  double m=ceil(log(meansh/meansl)/log(2));
  NumericVector parameters(m+1);
  for(int i=0;i<=m;i++)
    parameters[i]=meansh/pow(2,i+0.5);
  return parameters;
}

// [[Rcpp::export]]
List CGKPRESET(int n,int d,NumericVector parameters){
  int l=parameters.length();
  List preset;
  for(int i=0;i<l;i++){
    List subpreset;
    NumericVector s1v(n);
    double np=1/(n*parameters[i]);
    for(int j=0;j<n;j++){
      double x=j*np;
      double y=-x*x/2;
      s1v[j]=exp(y);
    }
    subpreset.push_back(s1v);
    NumericVector s1vc=cumsum(s1v);
    NumericVector s2v(n);
    for(int j=0;j<n;j++)
      s2v[j]=(s1vc[j]+s1vc[n-1-j]-1)/n;
    subpreset.push_back(s2v);
    double v1=0;
    for(int j=0;j<n;j++){
      double s1vd=pow(s1v[j],d);
      v1=v1+(n-j)*s1vd;
    }
    v1=(2*v1-n)/(n*n);
    subpreset.push_back(v1);
    double v2=0;
    for(int j=0;j<n;j++){
      double s2vd=pow(s2v[j],d);
      v2=v2+s2vd;
    }
    v2=v2/n;
    subpreset.push_back(v2);
    double v3=0;
    for(int j=0;j<n;j++){
      v3=v3+(n-j)*s1v[j];
    }
    v3=(2*v3-n)/(n*n);
    v3=pow(v3,d);
    subpreset.push_back(v3);
    preset.push_back(subpreset);
  }
  return preset;
}

// [[Rcpp::export]]
NumericVector CGK(NumericMatrix X,int n,int d, List preset){
  Function rankC("rank");
  String str("random");
  int pl=preset.length();
  NumericVector VAL(pl);
  for(int xpl=0;xpl<pl;xpl++){
    List subpreset=preset[xpl];
    NumericVector s1v=subpreset[0];
    NumericVector s2v=subpreset[1];
    double v1=subpreset[2];
    double v2=subpreset[3];
    double v3=subpreset[4];
    NumericMatrix Y(n,d);
    for(int i=0;i<d;i++){
      NumericVector vec1=X(_,i);
      NumericVector vec2=rankC(vec1,_["ties.method"]=str);
      Y(_,i)=vec2;
    }
    IntegerVector distvec0=rep(1,n*(n-1)/2);
    NumericVector distvec=as<NumericVector>(distvec0);
    for(int k=0;k<d;k++){
      NumericVector vec=Y(_,k);
      int l=0;
      for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){
          int dif=vec[i]-vec[j];
          if(dif<0)
            dif=-dif;
          distvec[l]=distvec[l]*s1v[dif];
          l++;
        }
      }
    }
    double s1=sum(distvec);
    s1=(2*s1+n)/(n*n);
    double s2=0;
    for(int i=0;i<n;i++){
      double pd=1;
      for(int j=0;j<d;j++){
        pd=pd*s2v[Y(i,j)-1];
      }
      s2=s2+pd;
    }
    s2=s2/n;
    double cgksq=(s1-2*s2+v3)/(v1-2*v2+v3);
    if(cgksq>0){
      VAL[xpl]=sqrt(cgksq);
    }else{
      VAL[xpl]=0;
    }
  }
  return VAL;
}

// [[Rcpp::export]]
NumericVector MCGK(List D,int n,int d,List testpreset){
  NumericVector valsum(n);
  NumericVector valmax(n);
  for(int i=0;i<n;i++){
    NumericMatrix Y(n-1,d);
    for(int j=0;j<d;j++){
      NumericMatrix Dj=D[j];
      NumericVector dij=Dj(_,i);
      dij.erase(i);
      Y(_,j)=dij;
    }
    NumericVector val=CGK(Y,n-1,d,testpreset);
    valsum[i]=sum(val);
    valmax[i]=max(val);
  }
  NumericVector mcgk(2);
  mcgk[0]=sum(valsum);
  mcgk[1]=sum(valmax);
  return mcgk;
}

// [[Rcpp::export]]
List MCGKTEST(List D,int n,int d,List testpreset,int iteration){
  Function sampleC("sample.int");
  NumericVector tstat=MCGK(D,n,d,testpreset);
  NumericMatrix ND(2,iteration);
  for(int itr=0;itr<iteration;itr++){
    List E;
    for(int k=0;k<d;k++){
      NumericMatrix Dk=D[k];
      NumericMatrix Ek(n,n);
      IntegerVector ind=sampleC(n);
      ind=ind-1;
      for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){
          Ek(i,j)=Dk(ind(i),ind(j));
          Ek(j,i)=Ek(i,j);
        }
      }
      E.push_back(Ek);
    }
    NumericVector ndist=MCGK(E,n,d,testpreset);
    ND(_,itr)=ndist;
  }
  List result=List::create(tstat,ND);
  return result;
}
