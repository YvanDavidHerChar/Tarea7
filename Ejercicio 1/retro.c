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
    array[i] = i+1;
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

float *modelo( float *x_obs, float *y_modelo, float a, float b, float c, float d, int n_puntos)
{
  int i;

  for(i=0; i < n_puntos; i++)
    {
      y_modelo[i] = a*cos((2*pi/d)*x_obs[i] + b) + c ; 
    }

  return y_modelo;
}
void *sacarLosDatos(float *manchasPorAnho)
{
  FILE *in;
  float sumaDelAnho=0;
  int i;
  int j;
  int laUnidad=0;
  char linea[4096];
  char filename[100]= "monthrg.dat";

  in = fopen(filename, "r");

  while(fgets(linea,sizeof(linea),in) != 0)
    {
      float manchas,anho, mes, dia, esteTrabajovaleun5;
      i++;

      if(i>2280)
	{
	 if(j<=10)
	    {
	      sscanf(linea, "%f %f %f %f %f", &anho, &mes, &dia, &manchas, &esteTrabajovaleun5);
	      sumaDelAnho = sumaDelAnho + manchas;
	      j++; 
	    }
	  else
	    {

	      sscanf(linea, "%f %f %f %f %f", &anho, &mes, &dia, &manchas, &esteTrabajovaleun5);
	      sumaDelAnho = sumaDelAnho + manchas;
	      j++; 

	      j=0;
	      manchasPorAnho[laUnidad] = sumaDelAnho;
	      sumaDelAnho = 0;
	      laUnidad++;

	    }
	}
    }

  fclose(in);
}

void estadoInicial(float *a, float *b, float *c, float *d, float *l, float *yInicial,float *x_obs, float *y_obs, int n_puntos)
{
  srand48(time(NULL));
  a[0] = drand48();
  b[0] = drand48();
  c[0] = drand48();
  d[0] = drand48();

  yInicial = modelo(x_obs, yInicial,a[0],b[0],c[0],d[0],n_puntos);

  l[0] =  likelyhood(y_obs, yInicial, n_puntos);

}

void *elCamino(float *x_obs, float *y_obs, float *yDelMomento,float *yCandidato, float *a, float *b, float *c, float *d,float *l, int n_puntos, int iteraciones, int burn)
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

      yDelMomento = modelo(x_obs, yDelMomento ,a[i],b[i],c[i],d[i],n_puntos);
      yCandidato = modelo(x_obs, yCandidato ,aPrima,b[i],c[i],d[i],n_puntos);

      lPresente = likelyhood(y_obs, yDelMomento, n_puntos);
      lCandidato = likelyhood(y_obs, yCandidato , n_puntos);
 
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

      yDelMomento = modelo(x_obs, yDelMomento ,a[i],b[i],c[i],d[i],n_puntos);
      yCandidato = modelo(x_obs, yCandidato ,a[i],bPrima,c[i],d[i],n_puntos);

      lPresente = likelyhood(y_obs, yDelMomento, n_puntos);
      lCandidato = likelyhood(y_obs, yCandidato , n_puntos);

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

      yDelMomento = modelo(x_obs, yDelMomento ,a[i],b[i],c[i],d[i],n_puntos);
      yCandidato = modelo(x_obs, yCandidato ,a[i],b[i],cPrima,d[i],n_puntos);

      lPresente = likelyhood(y_obs, yDelMomento, n_puntos);
      lCandidato = likelyhood(y_obs, yCandidato , n_puntos);

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

      yDelMomento = modelo(x_obs, yDelMomento ,a[i],b[i],c[i],d[i],n_puntos);
      yCandidato = modelo(x_obs, yCandidato ,a[i],b[i],c[i],dPrima,n_puntos);

      lPresente = likelyhood(y_obs, yDelMomento, n_puntos);
      lCandidato = likelyhood(y_obs, yCandidato , n_puntos);

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

      yCandidato = modelo(x_obs, yCandidato, a[i+1],b[i+1],c[i+1],d[i+1],n_puntos);

      l[i+1] = likelyhood(y_obs,yCandidato,n_puntos);

      if(i+1>burn)
	{
	  printf("%f %f %f %f %f\n", a[i+1],b[i+1],c[i+1],d[i+1], l[i+1]);
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
  float *yDelMomento;
  float *yCandidato;
  float *x_obs;
  float *y_obs;
  int n_puntos;
  int iteraciones;
  int i;
  int burn;

  iteraciones = atoi(argv[1]);
  burn = atoi(argv[2]);
  n_puntos = 1995-1800+1;

  x_obs = matriz(n_puntos);
  y_obs = matriz(n_puntos);
  a_walk = matriz(iteraciones);
  b_walk = matriz(iteraciones);
  c_walk = matriz(iteraciones);
  d_walk = matriz(iteraciones);
  yInicial = matriz(n_puntos);
  yDelMomento = matriz(n_puntos);
  yCandidato = matriz(n_puntos);
  l_walk = matriz(iteraciones);



  sacarLosDatos(y_obs);

  estadoInicial(a_walk, b_walk, c_walk, d_walk, l_walk, yInicial, x_obs, y_obs, n_puntos);

  elCamino(x_obs, y_obs, yDelMomento,yCandidato, a_walk, b_walk, c_walk, d_walk, l_walk, n_puntos, iteraciones, burn);
}
