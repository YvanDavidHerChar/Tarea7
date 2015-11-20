#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define USAGE "./random_gauss.x n_points"
#define pi 3.14159265359

float *matriz(int n_puntos)
{
  float *array;
  int i;

  array = malloc(n_puntos*sizeof(float));

  for(i=0;i<n_puntos;i++)
  {
    array[i] = 0.0;
  }
  return array;
}

float likelyhood(float *y_obs, float *y_model, int n_puntos)
{
  float suma=0;
  float resultado;
  int i;


  for(i = 0 ; i< n_puntos ; i++)
    {
      suma = suma + (y_obs[i]-y_model[i])*(y_obs[i]-y_model[i]);
    }

  resultado = (-1.0)*suma/2.0;

  return resultado;
}

void *sacarLosDatos(float *x_obs, float *y_obs)
{
  FILE *in;
  int laUnidad=-1;
  char linea[4096];
  char filename[100]= "lotka_volterra_obs.dat";

  in = fopen(filename, "r");

  while(fgets(linea,sizeof(linea),in) != 0)
    {
      float tiempo ,presas, predadores;
 
      sscanf(linea, "%f %f %f", &tiempo, &presas, &predadores);

      x_obs[laUnidad] = presas;
      y_obs[laUnidad] = predadores;
     laUnidad++;
    }

  fclose(in);
}

float fPrima(float x, float y, float alpha, float beta)
{
  float resultado;

  resultado = x*(alpha-beta*y);

  return resultado;
}

float gPrima(float x, float y, float gamma, float delta)
{
  float resultado;

  resultado = -y*(gamma-delta*x);

  return resultado;
}

float *RungeKutta(float *x, float *y, float x_0, float y_0, float alpha, float beta, float gamma, float delta, int iteraciones, float paso)
{
  float k1,k2,k3,k4;
  float l1,l2,l3,l4;
  float promedioK, promedioL;
  int t;

  x[0] = x_0;
  y[0] = y_0;

  for(t=0; t<iteraciones-1 ; t++)
    {
      k1 = fPrima(x[t], y[t], alpha,beta);
      l1 = gPrima(x[t], y[t], gamma, delta);

      k2 = fPrima(x[t] + 1/2*paso*k1, y[t] + 1/2*paso*l1, alpha, beta); 
      l2 = gPrima(x[t] + 1/2*paso*k1, y[t] + 1/2*paso*l1, gamma, delta); 

      k3 = fPrima(x[t] + 1/2*paso*k2, y[t] + 1/2*paso*l2, alpha, beta); 
      l3 = gPrima(x[t] + 1/2*paso*k2, y[t] + 1/2*paso*l2, gamma, delta);
 
      k4 = fPrima(x[t] + paso*k3, y[t] + paso*l3, alpha, beta); 
      l4 = gPrima(x[t] + paso*k3, y[t] + paso*l3, gamma, delta);

      promedioK = paso*(k1+ 2*k2 + 2*k3 + k4)/6; 
      promedioL = paso*(l1+ 2*l2 + 2*l3 + l4)/6;

      x[t+1] = x[t] + promedioK;
      y[t+1] = y[t] + promedioL;

    }

  return x;
}


void estadoInicial(float *a, float *b, float *c, float *d, float *l, float *yInicial, float *x_obs, float *x, float *y, int n_puntos, float y_0, float x_0, float paso)
{
  int i;
  srand48(time(NULL));
  a[0] = drand48();
  b[0] = drand48();
  c[0] = drand48();
  d[0] = drand48();

  yInicial = RungeKutta(x, y, x_0, y_0, a[0], b[0], c[0], d[0], n_puntos, paso);

  l[0] =  likelyhood(x_obs, yInicial, n_puntos);


}

void *elCamino(float *x_obs, float *x_aComparar,float *y_ayudaAx, float *xDelMomento,float *xCandidato, float *a, float *b, float *c, float *d,float *l,float paso, int t, int iteraciones, int burn)
{
  int i;
  float aPrima;
  float bPrima;
  float cPrima;
  float dPrima;
  float lPresente;
  float lCandidato;
  float gamma;
  float beta;
  float alpha;
  float x_0 = 15 ;
  float y_0 = 13; 

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for(i=0; i < (iteraciones-1) ; i++)
    {
      aPrima = gsl_ran_gaussian(r, 0.1)+a[i];
      //bPrima = gsl_ran_gaussian(r, 0.1)+b[i];
      //cPrima = gsl_ran_gaussian(r, 0.1)+c[i];
      //dPrima = gsl_ran_gaussian(r, 0.1)+d[i];

      xDelMomento = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i], b[i], c[i], d[i], t, paso);
      xCandidato = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, aPrima, b[i], c[i], d[i], t, paso);

      lPresente = likelyhood(x_obs, xDelMomento, t);
      lCandidato = likelyhood(x_obs, xCandidato , t);
 
      gamma = (lCandidato - lPresente);
      if(gamma>=0.00)
	{
	  a[i+1]=aPrima;
	}
      else
	{
	  beta = drand48();
	  alpha = exp(gamma);
	  if(beta<=alpha)
	    {
	      a[i+1]=aPrima;
	    }
	  else
	    {
	      a[i+1]=a[i];
	    }
	}
      //aPrima = gsl_ran_gaussian(r, 0.1)+a[i];
      bPrima = gsl_ran_gaussian(r, 0.1)+b[i];
      //cPrima = gsl_ran_gaussian(r, 0.1)+c[i];
      //dPrima = gsl_ran_gaussian(r, 0.1)+d[i];

      xDelMomento = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i], b[i], c[i], d[i], t, paso);
      xCandidato = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i], bPrima, c[i], d[i], t, paso);

      lPresente = likelyhood(x_obs, xDelMomento, t);
      lCandidato = likelyhood(x_obs, xCandidato , t);

      gamma = (lCandidato - lPresente);
      if(gamma>=0.00)
	{
	  b[i+1]=bPrima;
	}
      else
	{
	  beta = drand48();
	  alpha = exp(gamma);
	  if(beta<=alpha)
	    {
	      b[i+1]=bPrima;
	    }
	  else
	    {
	      b[i+1]=b[i];
	    }
	}

      //aPrima = gsl_ran_gaussian(r, 0.1)+a[i];
      //bPrima = gsl_ran_gaussian(r, 0.1)+b[i];
      cPrima = gsl_ran_gaussian(r, 0.1)+c[i];
      //dPrima = gsl_ran_gaussian(r, 0.1)+d[i];


      xDelMomento = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i], b[i], c[i], d[i], t, paso);
      xCandidato = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i], b[i], cPrima, d[i], t, paso);

      lPresente = likelyhood(x_obs, xDelMomento, t);
      lCandidato = likelyhood(x_obs, xCandidato , t);

      gamma = (lCandidato - lPresente);
      if(gamma>=0.00)
	{
	  c[i+1]=cPrima;
	}
      else
	{
	  beta = drand48();
	  alpha = exp(gamma);
	  if(beta<=alpha)
	    {
	      c[i+1]=cPrima;
	    }
	  else
	    {
	      c[i+1]=c[i];
	    }
	}
      //aPrima = gsl_ran_gaussian(r, 0.1)+a[i];
      //bPrima = gsl_ran_gaussian(r, 0.1)+b[i];
      //cPrima = gsl_ran_gaussian(r, 0.1)+c[i];
      dPrima = gsl_ran_gaussian(r, 0.1)+d[i];


      xDelMomento = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i], b[i], c[i], d[i], t, paso);
      xCandidato = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i], b[i], c[i], dPrima, t, paso);

      lPresente = likelyhood(x_obs, xDelMomento, t);
      lCandidato = likelyhood(x_obs, xCandidato , t);

      gamma = (lCandidato - lPresente);
      if(gamma>=0.00)
	{
	  d[i+1]=dPrima;
	}
      else
	{
	  beta = drand48();
	  alpha = exp(gamma);
	  if(beta<=alpha)
	    {
	      d[i+1]=dPrima;
	    }
	  else
	    {
	      d[i+1]=d[i];
	    }
	}

      xCandidato = RungeKutta(x_aComparar, y_ayudaAx, x_0, y_0, a[i+1], b[i+1], c[i+1], d[i+1], t, paso);

      l[i+1] = likelyhood(x_obs,xCandidato,t);

      if(i+1>burn)
	{
	 printf("%f %f %f %f %f\n", a[i+1],b[i+1],c[i+1],d[i+1],l[i+1]);
	}
    }
}

int main(int argc, char **argv)
{
  float *a_walk;
  float *b_walk;
  float *c_walk;
  float *d_walk;
  float *l_walk;
  float *yInicial;
  float *xDelMomento;
  float *xCandidato;
  float *x_aComparar;
  float *y_ayudaAx;
  float *x_obs;
  float *y_obs;
  float *x;
  float *y;
  int n_puntos;
  int iteraciones;
  int i;
  int burn;
  float x_0,y_0;
  float paso;

  iteraciones = atoi(argv[1]);
  burn = atoi(argv[2]);
  n_puntos = 96;
  paso = 0.8/n_puntos;

  x = matriz(n_puntos);
  y = matriz(n_puntos);
  x_obs = matriz(n_puntos);
  y_obs = matriz(n_puntos);
  a_walk = matriz(iteraciones);
  b_walk = matriz(iteraciones);
  c_walk = matriz(iteraciones);
  d_walk = matriz(iteraciones);
  yInicial = matriz(n_puntos);
  xDelMomento = matriz(n_puntos);
  xCandidato = matriz(n_puntos);
  l_walk = matriz(iteraciones);
  x_aComparar = matriz(n_puntos);
  y_ayudaAx = matriz(n_puntos);



  sacarLosDatos(x_obs, y_obs);

  x_0 = x_obs[0];
  y_0 = y_obs[0];

  estadoInicial(a_walk, b_walk, c_walk, d_walk, l_walk, yInicial, x_obs, x, y, n_puntos, y_0, x_0, paso);

  elCamino(x_obs, x_aComparar, y_ayudaAx, xDelMomento, xCandidato, a_walk, b_walk, c_walk, d_walk, l_walk, paso, n_puntos, iteraciones, burn);
}
