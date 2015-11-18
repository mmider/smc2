#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <stdio.h>

using namespace Rcpp;
using namespace std;

double runif()
{
  double U = rand()/((double)RAND_MAX);
  return U;
}

double rnormC(double mean = 0, double sd = 1)
{
  double U_1 = runif();
  double U_2 = runif();
  double Z = sqrt(-2 * log(U_1)) * cos(2 * PI * U_2);
  Z = sd * Z + mean;
  return Z;
}

double dnormC(double x, double mean = 0, double sd = 1)
{
  double denom = 1/(sqrt(2 * PI) * sd);
  double exponent = -(x-mean) * (x-mean)/(2*sd*sd);
  return denom * exp(exponent);
}

int sampleC(vector<double> probs)
{
  int n = probs.size();
  double total = probs[0];
  double realisation = runif();
  for (int i = 0; i < n; i++){
    if (realisation < total){
      return i;
    }
    total += probs[i+1];
  }
  return -1;
}

// [[Rcpp::export]]
int sample_wrap(NumericVector probs){
  int n = probs.size();
  vector<double> prob(n);
  for (int i = 0; i < n; i++){
    prob[i] = probs(i);
  }
  return sampleC(prob);
}

// [[Rcpp::export]]
NumericMatrix particle_filter(double sigma, double tau, NumericVector y, int Nx)
{
  int T = y.size();
  NumericMatrix history(Nx, 3 * T - 1);
  vector<double> x(Nx);
  vector<double> x_temp(Nx);
  vector<double> weights(Nx);
  vector<int> a(Nx);
  vector<double> W(Nx);
  double temp;
  
  for(int i = 0; i < Nx; i++){
    temp = rnormC(0,sigma);
    x[i] = temp;
    history(i,0) = temp;
  }
  double total_weights = 0;
  for (int i = 0; i < Nx; i++){
    temp = dnormC(y(1), x[i], tau);
    weights[i] = temp;
    history(i,T) = temp;
    total_weights += temp;
  }
  
  for (int t = 1; t < T; t++){
    for (int i = 0; i < Nx; i++){
      W[i] = weights[i]/total_weights;
    }
    total_weights = 0;
    for (int i = 0; i < Nx; i++){
      a[i] = sampleC(W);
      history(i, 2*T+t-1) = a[i];
      x_temp[i] = rnormC(x[a[i]],sigma);
    }
    for (int i = 0; i < Nx; i++){
      x[i] = x_temp[i];
      weights[i] = dnormC(y(t), x_temp[i], tau);
      if (weights[i] < 0){
        printf("ABORT!!!!      %f\n", weights[i]);
        printf("y_t: %f\nx_i: %f\ntau: %f\n", y(t), x_temp[i], tau);
      }
      history(i,t) = x_temp[i];
      history(i,T+t) = weights[i];
      total_weights += weights[i];
    }
  }
  return history;
}
