
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>

#include "system.h"
#include "utility.h"

using namespace std;

double dotProd(double** a, double** b, int part, int dim);

void getCOM(double** q, double com[], double* m, double mTot, 
		int numPart, int dim);

void getRadius(double** q, double* rad, int numPart, int dim);

void getTheta(double** q, double* theta, int numPart, int dim);

void matrixTranspose(double** a, double** b, int dim);

void nmTransform(double** a, double** b, double** nmTransform, 
		int numPart, int dim, int dir);

void nmSetup(double* nmFreq, double** nmTrans, double beta, int numPart, int dim);

void freeRP(double** qNM, double** pNM, double* nmFreq, double* m, 
		double mass, double tstep, int numPart, int dim);

void pileSetup(double* pileC1, double* pileC2, double* nmFreq, 
		double tauT, double tstep, int numPart);

void pile(double** pNM, double* pileC1, double* pileC2, double* m, double beta,
		int numPart, int dim);

void getqCent(double* rad, double* theta, double* qCent, int numPart, int dim);

double getDist(double* a, double* b, double* dif, int dim);

double getBeadPot(double** q, double k, int numPart, int dim);

void getBeadForce(double** q, double** force, double kBead, 
		int numPart, int dim);

double getg(double* r, double* qCent, double k, double rCon, 
		int conCase);

void getGradg(double** q, double* com, double* rad,double*** grad, 
		double * m, double mTot, double k, 
		int numPart, int dim, int conType, int conCase);

void rattleR(int maxiter, double** q, double** p, double* lambdaR, 
		double k, double* rCon, double tol, double* m, double mass, double mTot, double tstep, 
		int numPart, int dim, int* constraints, int numCon);

void rattleV(int maxiter, double** q, double** p, double* lambdaV,
		double k, double tol, double* m, double mTot, double tstep,
		int numPart, int dim, int* constraints,int numCon);

void splineSetup(double *xv, double *yv, double *y2, double yp1, double ypn, int n);

double interp(double *xx, double *yy, double *y2, double x, int j);

double locate(double *xx, double x,int n, int m);

double window(double t, double tHalf, double tau);

double trapezoid(double* x, double* y, int n);

double getFourierPot(double*** a, double kbead, int numPart, int kmax, int dim);

double getXiPot(double** q, double*** a, double*** qxi, double** xiPot, double rEq, 
    double dxi, int numPart, int kmax, int dim, int nxi, int potType);

void histogram(double val, double** hist, double xMin, double xMax, int numBins);

//void forceHistogram(double val, double force,  double** hist, double xMin, double xMax, 
//    int numBins, int dim);

void forceHistogram(double val, double* pos, double* force,  double** hist, 
    double xMin, double xMax, int numBins, int dim);

