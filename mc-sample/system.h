
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <fftw3.h>
#include <complex.h>

using namespace std;

double getPot(double* rad, double rEq,
	 	int numPart, int potType);

double getPotSingle(double rad, double rEq, int potType);

void getForce(double* rad, double** q, double** force, double* fScalar, double rEq, 
    int numPart, int dim, int potType);

void getQCentForce(double qCentRad, double* qCentQ, double* qCentForce, double fScalar,
		int numPart, int dim);


