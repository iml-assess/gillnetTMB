#include <TMB.hpp> 
#include <iostream>

template<class Type>
  Type objective_function<Type>::operator() ()
{

  // input data
  DATA_VECTOR(length)
  DATA_VECTOR(mesh)
  DATA_VECTOR(cpn)
  DATA_INTEGER(rtype)
  DATA_INTEGER(distr)
  DATA_IVECTOR(Nid)

  // parameters
  PARAMETER_VECTOR(par)
  PARAMETER_VECTOR(logN)

  // basics
  int n=cpn.size();
  Type nll=0;
  Type one=1;
  Type two=2;

  // selectivity
  vector<Type> logsel(n);
  Type k1=par(0);
  Type k2=par(1);
  Type m1=mesh.minCoeff();

  switch(rtype){
  	case 1:  //norm.loc
      logsel = -pow(length -k1 * mesh, two) / (two * pow(k2, two));
  	  break;

  	case 2:  //norm.sca
  	 logsel = -pow(length - k1 * mesh, two) / (two * k2 * pow(mesh, two));
  	 break;

  	case 3:  //gamma
  	 logsel = (k1 - one) * (log(length) - log((k1 - one) * k2 * mesh)) + k1 - one - length / (k2 * mesh);
  	 break;

  	case 4:  //lognorm
  	 logsel = log(one / length) + k1 + log(mesh / m1) - pow(k2, two) / two - (pow(log(length) - k1 - log(mesh / m1), two)) / (two * pow(k2, two));
  	 break;

  	default:
  	  error("Unknown rtype code");
      return(0);
      break;
  }
  vector<Type> sel = exp(logsel);

  // predictions
   vector<Type> logpred(n);
   int ix;
   int onei=1;

   for(int i=0; i<n; i++){
    ix = Nid(i) - onei;
  	logpred(i) = logN(ix) + logsel(i);
   }

  // likelihood
  vector<Type> pred = exp(logpred);
  Type theta = exp(par(par.size()-onei));

  switch(distr){
  	case 1: // poisson
		  nll -= sum(dpois(cpn, pred, true));
  	  break; 

  	case 2: // nbinom
      for(int i=0; i<n; i++){
       nll -= theta * (log(theta) - log(theta+pred(i))) + lgamma(theta+cpn(i)) - lgamma(cpn(i)+one)- lgamma(theta) + cpn(i) * (log(pred(i)) - log(theta+pred(i)));  
     }
      //nll -= sum(dnbinom(cpn, theta, theta / (theta + pred), true)); // size than prob, or inverse -> compilation error (arguments of function don't match function?)
  	  break;

  	default:
  	  error("Unknown error distribution code");
      return(0);
      break;
  }

  // return
  ADREPORT(logpred)
  ADREPORT(logsel)
    
  return nll;
}
