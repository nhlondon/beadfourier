
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <fftw3.h>
#include <complex.h>

#include "system.h"
using namespace std;

double getPot(double* rad, double rEq, int numPart, int potType)
{
	double arg,arg2,pot;

	pot = 0.0;
	for(int i=0; i<numPart; i++)
	{
	  arg = rad[i] - rEq;
    arg2 = arg*arg;
    switch(potType)
    {
      case 1:
        pot += 0.5*arg2;
        break;

      case 2:
        pot += 0.5*arg2 + 0.1*arg2*arg + 0.01*arg2*arg2;
        break;

      case 3:
        pot += 0.25*arg2*arg2;
        break;

      case 4:
        arg = 1.0 - exp(-1.1605*arg);
        pot += 0.18748*arg*arg;
        break;
      
      case 5:
        pot += 0.5*arg2 + 1.21/2.0*arg2*arg + 7.0*1.4641/24.0*arg2*arg2;
        break;
    }
  }

	return pot;
}

double getPotSingle(double rad, double rEq, int potType)
{
	double arg,arg2,pot;

	pot = 0.0;
  arg = rad - rEq;
  arg2 = arg*arg;
  switch(potType)
  {
    case 1:
      pot += 0.5*arg2;
      break;

    case 2:
      pot += 0.5*arg2 + 0.1*arg2*arg + 0.01*arg2*arg2;
      break;

    case 3:
      pot += 0.25*arg2*arg2;
      break;

    case 4:
      arg = 1.0 - exp(-1.1605*arg);
      pot += 0.18748*arg*arg;
      break;
    
    case 5:
      pot += 0.5*arg2 + 1.21/2.0*arg2*arg + 7.0*1.4641/24.0*arg2*arg2;
      break;
  }

	return pot;
}
void getForce(double* rad, double** q, double** force, double* fScalar, double rEq, 
    int numPart, int dim, int potType)
{
	double arg,arg2;

	for(int i=0; i<numPart; i++)
	{
    arg = rad[i]-rEq;
    arg2 = arg*arg;

    switch(potType)
    {
      case 1:
        fScalar[i] = -1.0*arg;
        break;

      case 2:
        fScalar[i] = -1.0*arg - 0.3*arg2 - 0.04*arg2*arg;
        break;

      case 3:
        fScalar[i] = -1.0*arg2*arg;
        break;
    
      case 4:
        arg = exp(-1.1605*arg);
        fScalar[i] = -2.0*0.18748*1.1605*(1.0-arg)*arg;
        break;
      
      case 5:
        fScalar[i] = -1.0*arg - 1.21*3.0/2.0*arg2 - 7.0*1.2641/6.0*arg2*arg;
        break;
    }
    
    if(dim > 1)
    {
      for(int j=0; j<dim; j++)
      {
        force[i][j] += q[i][j]/rad[i]*fScalar[i];
      }
    }
    else
      force[i][0] += fScalar[i];
	}	
	return;
}

void getQCentForce(double qCentR, double* qCentQ, double* qCentForce, double fScalar,
		int numPart, int dim)
{

		for(int i=0; i<dim; i++)
		{
			qCentForce[i] = qCentQ[i]/qCentR * fScalar;
		}
		
		return;
}
