
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <fftw3.h>

#include "functions.h"
#include "utility.h"

using namespace std;

double dotProd(double** a, double** b, int part, int dim)
{
	double c=0.0;

	for(int i=0;i<dim;i++)
	{
		c+=a[part][i]*b[part][i];
	}

	return c;
}

void getCOM(double** q, double* com, double* m, double mTot, 
		int numPart, int dim)
{
	//static double com[dim];
	for(int i=0;i<dim;i++)
	{
		com[i]=0.0;
		for(int j=0;j<numPart;j++)
		{
			com[i] += m[j]*q[j][i];
		}
		com[i] = com[i]/mTot;
	}

	return;
}

void getRadius(double** q, double* rad, int numPart, int dim)
{
	for(int i=0; i<numPart; i++)
	{
		rad[i]=0.0;
    if(dim > 1)
    {
      for(int j=0; j<dim; j++)
        rad[i] += q[i][j]*q[i][j];
	    rad[i] = sqrt(rad[i]);
		}
    else
    {
      rad[i]=q[i][0];
    }
	}

	return;
}

void getTheta(double** q, double* theta, int numPart, int dim)
{
	int dir;
	double pi,dif;
	bool wrap;
	pi = 4.0*atan(1.0);
	wrap = false;

	for(int i=0; i<numPart; i++)
	{
		theta[i] = atan2(q[i][1],q[i][0]);
		//cout << "theta: " << i << " " << theta[i] << endl;	
	}

	/*theta[0] = theta[0]/(2.0*pi);
	theta[0] = theta[0] - round(theta[0]);
	theta[0] = theta[0]*pi*2.0;
	*/
	for(int i=1; i<numPart; i++)
	{
//		if(!first)
//		{
			dif = theta[i] - theta[0];
			//dif = theta[i] - thetaOld[i];
			dif = dif/(2.0*pi);
			dif = dif - round(dif);
			dif = dif*pi*2.0;
			theta[i] = theta[0] + dif;
		//cout << "theta: " << i << " " << theta[i] << endl;	
			//theta[i] = thetaOld[i] + dif;
//		}
	}

		/*	if(!first)
		{
			if(abs(theta[i]-thetaOld[i]) > pi)
			{
				if(theta[i] > 0.0) theta[i] -= 2*pi;
				else theta[i] += 2*pi;
			}
		}*/
//		thetaOld[i] = theta[i];
//	}
/*
	for(int i=0; i<numPart; i++)
	{
		if(abs(theta[i]) > 2.0*pi)
		{
			wrap = true;
			cout << "wrap" << endl;
			if(theta[i] > 0.0) dir = -1;
			else dir = 1;
			break;
		}
	}

	if(wrap)
	{
		for(int i=0; i<numPart; i++)
		{
			theta[i] += dir*2*pi;
			thetaOld[i] = theta[i];
		}
	}
*/
	return;
}

void getqCent(double* rad, double* theta, double* qCent, int numPart, int dim)
{
	double rAvg, thetaAvg;

	rAvg = 0.0;
	thetaAvg = 0.0;

	for(int i=0; i<numPart; i++)
	{
		rAvg += rad[i];
		if(dim > 1)
      thetaAvg += theta[i];
	//	cout << "angle: " << theta << endl;
	}

	rAvg = rAvg/numPart;
	thetaAvg = thetaAvg/numPart;

//	cout << "rAvg: " << rAvg << endl;
//	cout << "thetaAvg: " << thetaAvg << endl;
	if(dim > 1)
  {
    qCent[0] = rAvg*cos(thetaAvg);
    qCent[1] = rAvg*sin(thetaAvg);
    qCent[2] = rAvg;
	  qCent[3] = thetaAvg;
  }
  else
  {
    qCent[0] = rAvg;
  }

	return;
}

void nmTransform(double** a, double** b, double** nmTrans, int numPart, int dim, int dir)
{

	double nmTransT[numPart], aTmp[numPart];

	//writeOutArray(a,"a",numPart*dim);
	//Forward transformation
	if(dir == 1)
	{
		for(int i=0; i<dim; i++)
		{
			for(int j=0; j<numPart; j++)
			{
				b[j][i] = 0.0;
				for(int k=0; k<numPart; k++)
				{
					b[j][i] += nmTrans[k][j]*a[k][i];
				//	nmTransT[k] = nmTrans[k][j];
				//	aTmp[k] = a[k*dim+i];
				}
				//b[i*numPart+j] = dotProd( nmTransT, aTmp, numPart);
			}
		}
	}
	//Backwards transformation
	else if(dir == -1)
	{
		for(int i=0; i<dim; i++)
		{
			for(int j=0; j<numPart; j++)
			{
				//b[j*dim+i] = dotProd( nmTrans[j], &a[i*numPart], numPart);
				b[j][i] = 0.0;
				for(int k=0; k<numPart; k++)
				{
					b[j][i] += nmTrans[j][k]*a[k][i];
				}
			}
		}
	}
	//writeOutArray(b,"b",numPart*dim);

	return;
}	
void nmSetup(double* nmFreq, double** nmTrans, double beta, int numPart, int dim)
{
	double freq,pi;

	pi = 4.0*atan(1.0);
	freq = ((double)numPart)/beta;

	for(int i=1;i<numPart;i++)
	{
		nmFreq[i] = 2.0*freq*sin(i*pi/numPart);
	}
	nmFreq[0] = 0.0;

	//writeOutArray(nmFreq,"nmFreq",numPart);

	for(int j=0;j<numPart;j++) 
	{
		nmTrans[j][0] = 1.0/sqrt(numPart);
		
		for(int k=1;k<=numPart/2;k++) 
		{
			nmTrans[j][k] = sqrt(2.0/numPart)*cos(2*pi*(j+1)*k/numPart);
		}
		for(int k=(numPart/2)+1;k<numPart;k++) 
		{
			nmTrans[j][k] = sqrt(2.0/numPart)*sin(2*pi*(j+1)*k/numPart);
		}

		if(numPart%2==0)
		{
			nmTrans[j][numPart/2] = 1.0/sqrt(numPart)*pow(-1.0,j+1);
		}
	}

//	writeOutMatrix(nmTrans,"nmTrans",numPart,numPart);

	return;
}

void freeRP(double** qNM, double** pNM, double* nmFreq, double* m, 
		double mass, double tstep, int numPart, int dim)
{
	double omegaPrime, qTemp;

	for(int i=0; i<dim; i++)
	{
		qNM[0][i] = qNM[0][i] + pNM[0][i]*tstep/m[0];
	}

	if(numPart > 1)
	{
		for(int i=1; i<numPart; i++)
		{
			omegaPrime = sqrt(mass/m[i])*nmFreq[i];
			//writeOutDbl(nmFreq[i],"nmFreq");
			//writeOutDbl(omegaPrime,"omegaPrime");
			for(int j=0; j<dim; j++)
			{
				qTemp = qNM[i][j]*cos(omegaPrime*tstep) + 
					pNM[i][j]/(m[i]*omegaPrime)*sin(omegaPrime*tstep);
				pNM[i][j] = -m[i]*omegaPrime*qNM[i][j]*sin(omegaPrime*tstep) +
					pNM[i][j]*cos(omegaPrime*tstep);
				qNM[i][j] = qTemp;
			}
		}
	}

	return;
}

void pileSetup(double* pileC1, double* pileC2, double* nmFreq, 
		double tauT, double tstep, int numPart)
{
	double gamma[numPart];

	gamma[0] = 1.0/tauT;

	for(int i=1; i<numPart; i++)
	{
		gamma[i] = 2.0*nmFreq[i];
	}

	for(int i=0; i<numPart; i++)
	{
		pileC1[i] = exp(-(tstep/2.0)*gamma[i]);
		pileC2[i] = sqrt(1.0-(pileC1[i]*pileC1[i]));
	}

	//writeOutArray(gamma,"gamma",numPart);
	//writeOutArray(pileC1,"pileC1",numPart);
	//writeOutArray(pileC2,"pileC2",numPart);

	return;
}


void pile(double** pNM, double* pileC1, double* pileC2, double* m, double beta, 
		int numPart, int dim)
{
	double factor;

	//factor = sqrt(m[0]/beta);
	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<numPart; j++)
		{
			pNM[j][i] = pNM[j][i]*pileC1[j] + 
				sqrt(m[j]/beta)*pileC2[j]*randGaus(1.0,0.0);
		}
	}

	return;
}
/*
double getDist(double* a, double* b, double* dif, int dim)
{
	double dist=0.0;

	for(int i=0;i<dim;i++)
	{
		dif[i] = a[i] - b[i];
		dist += dif[i]*dif[i];
	}

	dist = sqrt(dist);

	return dist;
}
*/

double getBeadPot(double** q, double kBead, int numPart, int dim)
{
	int j;
	double dif, pot;
	
	pot = 0.0;
  if(numPart == 2)
  {
    for(int k=0; k<dim; k++)
    {
      dif = q[0][k] - q[1][k];
      pot += kBead*dif*dif;
    }
  }
  else
  {
    for(int i=0; i<numPart; i++)
    {
    //	cout << i << endl;
      if(i==numPart-1) j=0;
      else j=i+1;

    //	cout << j << endl;
      for(int k=0; k<dim; k++)
      {
    //		cout << i*dim+k << " " << j*dim+k << endl;
        dif = q[i][k] - q[j][k];
  //			cout << dif << endl;
        pot += 0.5*kBead*dif*dif;
    //cout << "pot: " << pot << endl;
      }
    }
  }

	//cout << "pot: " << pot << endl;
	//cout << endl;
	return pot;
}

void getBeadForce(double** q, double** force, double kBead, int numPart, int dim)
{
	int j,k;
	double dif;

	for(int i=0; i<numPart; i++)
	{
		//cout << i << endl;
		if(numPart > 2)
		{
			if(i==0)
			{
				j=i+1;
				k=numPart-1;
			}
			else if(i==numPart-1)
			{
				j=0;
				k=i-1;
			}
			else
			{
				j=i+1;
				k=i-1;
			}
		}
		else
		{
			if(i==0)
			{
				j=i+1;
				k=i+1;
			}
			else
			{
				j=i-1;
				k=i-1;
			}
		}
		//cout << j << " " << k << endl;

		for(int l=0; l<dim; l++)
		{
			force[i][l] += -kBead*(2.0*q[i][l] - q[j][l] - q[k][l]);
		}
	}

	return;
}


double getg(double* q, double* qCent, double k, double rCon, int conCase)
{
	double g,pi;

	switch(conCase)
	{
		case 1:
//			cout << " g-COM: ";
			g = k*q[0]*q[0] + q[1]*q[1] - rCon*rCon;
//			cout << g << endl;
			break;
		case 2:
//			cout << " g-qcent: ";
			g = qCent[2] - rCon;
//			cout << g << endl;
			break;
		case 3:
//			cout << " g-qcent: ";
			pi = 4.0*atan(1.0);	
			g = qCent[3] - rCon;
			g = g/(2.0*pi);
			g = g - round(g);
			g = g*2.0*pi;
			//writeOutDbl(g,"g");
			//			cout << g << endl;
			break;
	}
//	cout << g << endl;
	return g;
}

void getGradg(double** q, double* com, double* rad, double*** grad, double* m, 
		double mTot, double k, int numPart, int dim, int conType, int conCase)
{
	//static double grad[numPart*dim];

	//cout << "in getGradg" << endl;
	//writeOutInt(conType,"conType");
	
	//conType = 1;

	switch(conType)
	{
		case 1:
		//	cout << "grad-COM: " << endl;
			for(int i=0;i<numPart;i++)
			{
				grad[conCase][i][0] = 2.0*k*com[0]*m[i]/mTot;
				grad[conCase][i][1] = 2.0*com[1]*m[i]/mTot;
			}
			break;
		case 2:
		//	cout << "grad-qCentR: " << endl;
			for(int i=0;i<numPart;i++)
			{
				grad[conCase][i][0] = q[i][0]/rad[i]/numPart;
				grad[conCase][i][1] = q[i][1]/rad[i]/numPart;
			}
				//writeOutMatrix(grad[conCase],"grad",numPart,dim);
			break;
		case 3:
			//cout << "grad-qCentTheta: " << endl;
			for(int i=0;i<numPart;i++)
			{
				//writeOutInt(i,"i");	
				double x2 = q[i][0]*q[i][0];
				double y2 = q[i][1]*q[i][1];
				double denom = x2 + y2;
				//writeOutDbl(denom,"denom");

				grad[conCase][i][0] =	- q[i][1]/denom/numPart;
				grad[conCase][i][1] =  q[i][0]/denom/numPart;	

		//		double denom = 1.0 + q[i*dim+1]*q[i*dim+1]/x2;

		//		grad[i*dim] = -q[i*dim+1]/(x2*denom)/numPart;
			//	grad[i*dim+1] = 1.0/(q[i*dim]*denom)/numPart;
				//writeOutMatrix(grad[conCase],"grad",numPart,dim);
			}
	}
	
	return;
}

void rattleR(int maxiter, double** q, double** p, double* lambdaR, 
		double k, double* rCon, double tol, double* m, double mass,
		double mTot, double tstep, 
		int numPart, int dim, int* constraints, int numCon)
{
	double g,update,conVal;
	double rPrime[numPart*dim];
	double com[dim],rad[numPart],theta[numPart],qCent[2*dim];
	double comOld[dim];
	double ***grad;
	double ***gradOld;
	//double gradOld[numPart*dim*numCon];

	grad = new double**[numCon];
	gradOld = new double**[numCon];
	for(int i=0; i<numCon; i++)
	{
		grad[i] = new double*[numPart];
		gradOld[i] = new double*[numPart];
	
		for(int j=0; j<numPart; j++)
		{
			grad[i][j] = new double[dim];
			gradOld[i][j] = new double[dim];
		}
	}

	int iter=0;
	int numCvgd=0,cvgd[numCon];
	//p[0] = pOld[0];
	//p[1] = pOld[1];	
	//lambdaR = 0.0;

	int conType;

	getCOM(q,comOld,m,mTot,numPart,dim);
	getRadius(q,rad,numPart,dim);
	getTheta(q,theta,numPart,dim);
	getqCent(rad,theta,qCent,numPart,dim);

	//cout << "RattleR" << endl;
	//writeOutArray(rCon,"rCon",numCon);
	for(int i=0; i<numCon; i++)
	{
		lambdaR[i] = 0.0;
		cvgd[i] = 0;
		conType = constraints[i];
		//writeOutInt(conType,"conType");	
		//cout << "get gradG" << endl;
		getGradg(q,comOld,rad,gradOld,m,mTot,k,numPart,dim,conType,i);
	//cout << "got gradG" << endl;
		//getGradg(r,comOld,rad,&gradOld[i*numPart*dim],m,mTot,k,numPart,dim,conCase);
	}
	
	while(numCvgd < numCon)
	{
		iter++;
		//writeOutInt(iter, "iter");
		for(int j=0; j<numCon; j++)
		{
			//writeOutInt(j,"j");
			if(cvgd[j] != 1)
			{
				conVal = rCon[j];
				conType = constraints[j];
				getCOM(q,com,m,mTot,numPart,dim);
				getRadius(q,rad,numPart,dim);
				getTheta(q,theta,numPart,dim);
				getqCent(rad,theta,qCent,numPart,dim);
				g = getg(com,qCent,k,conVal,conType);
			
				//writeOutDbl(g,"g");	
				if(conType == 2)
				{
					//writeOutArray(rad, "rad",numPart);
					//writeOutDbl(rCon[0],"radCon");
				}
				else if(conType == 3)
				{
					//writeOutArray(theta, "theta",numPart);
					//writeOutDbl(rCon[1],"radTheta");
				}

				if(abs(g) < tol)
				{
					//cout << "convergerd" << endl;
					cvgd[j] = 1;
					numCvgd++;
					continue;
				}
				getGradg(q,com,rad,grad,m,mTot,k,numPart,dim,conType,j);
				//getGradg(r,com,rad,&grad[j*numPart*dim],m,mTot,k,numPart,dim,conCase);

				//writeOutMatrix(gradOld[j],"gradOld",numPart,dim);
				//writeOutMatrix(grad[j],"grad",numPart,dim);
				update = 0.0;
				for(int k=0;k<numPart;k++)
				{
					//update -= (1.0/m[k]) * dotProd(&grad[k*dim+j*numPart*dim],
					//		&gradOld[k*dim+j*numPart*dim],dim);
				
					update -= (1.0/mass) * dotProd(grad[j],gradOld[j],k,dim);
					//writeOutDbl(update,"update");	
				}
				update = g/update;
				//writeOutDbl(update,"update");	
				for(int k=0;k<numPart;k++)
				{
					for(int l=0; l<dim; l++)
					{
						q[k][l] = q[k][l] + update*gradOld[j][k][l]/mass;
						p[k][l] = p[k][l] + update*gradOld[j][k][l]/tstep;
					}
				}

			//writeOutMatrix(q,"q",numPart,dim);
				lambdaR[j] -= update;
				//writeOutArray(lambdaR,"lambdaR",numCon);
			}
	
		}
		if(iter == maxiter) 
		{
			cout << "RattleR: Too many iterations" << endl;
			for(int j=0; j<numCon; j++)
			{
				if(cvgd[j] != 1)
					cout << "Constriant " << j << " not met" << endl;
			}
			writeOutMatrix(q,"q",numPart,dim);
			writeOutArray(rad,"rad",numPart);
			writeOutArray(theta,"theta",numPart);
			writeOutArray(rCon,"rCon",dim);
			exit(1);
		}
	}

	for(int j=0; j<numCon; j++)
	{
		for(int k=0; k<numPart; k++)
		{
			delete [] grad[j][k];
			delete [] gradOld[j][k];
		}
		delete grad[j];
		delete gradOld[j];
	}
	delete [] grad;
	delete [] gradOld;

	return;
	//return lambdaR;
}

void rattleV(int maxiter, double** q, double** p, double* lambdaV,
		double k, double tol, double* m, double mTot, double tstep,
		int numPart, int dim, int* constraints,int numCon)
{
	double update,dotSum;
	double rad[numPart],theta[numPart];
	double com[dim],qCent[2*dim];
	double dot[numPart];
	//double gradDot[numPart*numCon];
	double ***grad, **gradDot;

	grad = new double**[numCon];
	gradDot = new double*[numCon];
	for(int i=0; i<numCon; i++)
	{
		grad[i] = new double*[numPart];
		gradDot[i] = new double[numPart];
		for(int j=0; j<numPart; j++)
		{
			grad[i][j] = new double[dim];
		}
	}

	int iter=0;
	int conType;
	int numCvgd=0,cvgd[numCon];
	getCOM(q,com,m,mTot,numPart,dim);
	getRadius(q,rad,numPart,dim);
	getTheta(q,theta,numPart,dim);
	getqCent(rad,theta,qCent,numPart,dim);
	
	//cout << "RattleV" << endl;
	for(int i=0; i<numCon; i++)
	{
		lambdaV[i] = 0.0;
		cvgd[i] = 0;
		conType = constraints[i];
		getGradg(q,com,rad,grad,m,mTot,k,numPart,dim,conType,i);
	}
	
	//writeOutMatrix(p,"p",numPart,dim);
	//writeOutArray(grad,"grad",numPart*dim*numCon);
	for(int i=0;i<numCon;i++)
	{
		for(int j=0;j<numPart;j++)
		{
			gradDot[i][j] = dotProd(grad[i],grad[i],j,dim);
		}
	}

	//for(int i=0;i<maxiter;i++)
	while(numCvgd < numCon)
	{
		iter++;
		//cout << "iter: " << iter << endl;
		for(int j=0; j<numCon; j++)
		{
			//cout << "j: " << j << endl;
			if(cvgd[j] != 1)
			{
				conType = constraints[j];
				dotSum = 0.0;
				for(int k=0;k<numPart;k++)
				{
					dot[k] = dotProd(grad[j],p,k,dim);
					dotSum += dot[k];
				}
				//cout << "dotSum: " << dotSum << endl;
				if(abs(dotSum) < tol)
				{
					cvgd[j] = 1;
					numCvgd++;
					continue;
				}

				update = 0.0;
				for(int k=0; k<numPart; k++)
				{
					update -= gradDot[j][k];
				}
				update = dotSum/update;

				//cout << "update: " << update << endl;
				for(int k=0;k<numPart;k++)
				{
					for(int l=0; l<dim; l++)
					{
						p[k][l] = p[k][l] + update*grad[j][k][l];
					}
				}
				//writeOutArray(p,"p",numPart*dim);
				lambdaV[j] += update;
			}
		}
		//writeOutInt(numCvgd,"numCvgD");
		//writeOutInt(numCon,"numCon");
		//writeOutArray(p,"p",numPart*dim);
		//if(numCvgd == numCon) break;
		if(iter == maxiter) 
		{
			cout << "RattleV: Too many iterations" << endl;
			for(int j=0; j<numCon; j++)
			{
				if(cvgd[j] != 1)
					cout << "Constriant " << j << " not met" << endl;
			}
			//writeOutDbl(dotSum,"dotSum");
			//writeOutMatrix(q,numPart,dim);
			//writeOutArray(rCon,dim);
			exit(1);
		}
	}
	
	//cout << "lambdaV: " << lambdaV << endl;
//		cout << "iter: " << iter << endl;

//	return lambdaV;
	for(int j=0; j<numCon; j++)
	{
		for(int k=0; k<numPart; k++)
		{
			delete [] grad[j][k];
		}
		delete grad[j];
		delete gradDot[j];
	}
	delete[] grad;
	delete[] gradDot;

	return;	
}

void splineSetup(double *xv, double *yv, double *y2, double yp1, double ypn, int n)
{
	int i,k;
	double p,qn,sig,un;
	double *u;

	u = new double[n-1];

	if(yp1 > 0.99e99)
		y2[0] = u[0] = 0.0;
	else
	{
		y2[0] = -0.5;
		u[0] = (3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0]-yp1));
	}

	for(i=1; i<n-1; i++)
	{
		sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
		u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
	}

	if(ypn > 0.99e99)
		qn=un=0.0;
	else
	{
		qn=0.5;
		un=(3.0*u[i]/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
	}

	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];

	delete [] u;

	return;
}

double interp(double *xx, double *yy, double *y2, double x, int j)
{
	int klo,khi;
	double y,h,a,b;

	klo = j;
	khi = klo+1;
	//writeOutDbl(x,"x");
	//writeOutInt(klo,"klo");
	//writeOutDbl(xx[klo],"xx[klo]");

	h=xx[khi]-xx[klo];

	a=(xx[khi]-x)/h;
	b=(x-xx[klo])/h;

	y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0;

	return y;
}

double locate(double *xx, double x,int n, int m)
{
	int ju,jm,jl;
	bool ascnd=(xx[n-1] >= xx[0]);
	jl=0;
	ju=n-1;

	while (ju-jl > 1)
	{
		jm = (ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	return max(0,min(n-m,jl-((m-2)>>1)));
}

double window(double t, double tHalf, double tau)
{
	double win;

	win = exp((t-tHalf)/tau);
	win += 1.0;
	win = 1.0/win;

	return win;
}

double trapezoid(double* x, double* y, int n)
{
	double integrand;

	integrand = 0.0;
	for(int i=1; i<n; i++)
		integrand += ((y[i]+y[i-1])*abs(x[i]-x[i-1]))/2.0;

	return integrand;
}

double heaviside(double val, double ref)
{
  double heavi;

  heavi = 0.0;  
  if(val - ref >= 0)
    heavi = 1.0;

  return heavi;
}

double step(double val, double ref1, double ref2)
{
  double step;
  
  step = heaviside(val,ref1) - heaviside(val,ref2);

  return step;
}

void histogram(double val, double** hist, double xMin, double xMax, int numBins)
{
  int bin;
  double dBins, stepVal;

  dBins = (xMax-xMin)/numBins;
  
  bin = (int) ((val-xMin)/dBins);
  if(bin >=0 && bin < numBins)
    hist[1][bin] += 1.0/dBins;
  /*
  if(val - xMin > 0.0)
  {
    bin = (int) ((val-xMin)/dBins);
    if(bin >= 0 && bin < numBins)
      hist[1][bin] +=  1.0/dBins;
  }
  */
  return;
}

double getFourierPot(double*** a, double kbead, int numPart, int kmax, int dim)
{
  double pot,pi;

  pi = 4.0*atan(1.0);
  pot = 0.0;
  for(int i=0; i<numPart; i++)
  {
    for(int k=0; k<kmax; k++)
    {
      for(int j=0; j<dim; j++)
        pot += 0.25*kbead*((double)(k+1)*pi)*((double)(k+1)*pi)*a[i][k][j]*a[i][k][j];
    }
  }

  return pot; 
}
double getXiPot(double** q, double*** a, double*** qxi, double** xiPot, double rEq, 
    double dxi, int numPart, int kmax, int dim, int nxi, int potType)
{
  int iup;
  double pi,pot;

  double **rxi,*xi;

  pi = 4.0*atan(1.0);

  rxi = new double*[numPart];
  xi = new double[nxi];
 
  for(int i=0; i<numPart; i++)
  {
    rxi[i] = new double[nxi];
  }
  for(int i=0; i<nxi; i++)
    xi[i] = i*dxi;

  //writeOutArray(xi,"xi",nxi);
  for(int i=0; i<numPart; i++)
  {
    for(int l=0; l<nxi; l++)
    {
      rxi[i][l] = 0.0;
      for(int j=0; j<dim; j++)
      {
        if(numPart>1)
        {
          if(i==(numPart-1))
          {
            iup = 0;
          }
          else
          {
            iup = i+1;
          }
          qxi[i][j][l] = q[i][j] + (q[i][j]+q[iup][j])*xi[l];
        }
        else
        {
          qxi[i][j][l] = q[i][j];
        }
        for(int k=0; k<kmax; k++)
          qxi[i][j][l] += a[i][k][j]*sin((double)(k+1)*pi*xi[l]);
        rxi[i][l] += qxi[i][j][l]*qxi[i][j][l];
      }
      rxi[i][l] = sqrt(rxi[i][l]);
    }
  }
  //writeOutTensor(qxi,"qxi",numPart,dim,nxi);
  for(int i=0; i<numPart; i++)
  {
    for(int j=0; j<nxi; j++)
      xiPot[i][j] = getPotSingle(rxi[i][j],rEq,potType);
  }
  //writeOutMatrix(xiPot,"xiPot",numPart,nxi);
  pot = 0.0;
  for(int i=0; i<numPart; i++)
    pot += trapezoid(xi,xiPot[i],nxi);
  pot = pot/((double)numPart);

  //writeOutDbl(pot,"trapPot");
  for(int i=0; i<numPart; i++)
    delete [] rxi[i];
  delete [] rxi;
  delete [] xi;
  return pot;
}
void forceHistogram(double val, double* pos, double* force,  double** hist, 
    double xMin, double xMax, int numBins, int dim)
{
  int startBin,bin;
  double dBins, stepVal, dot;

  dBins = (hist[0][1]-hist[0][0]);
  
  dot = 0.0;
  for(int i=0; i<dim; i++)
    dot += pos[i]*force[i];

  bin = (int) ((val-xMin)/dBins);
  if(bin >=0 && bin < numBins)
    for(int i=bin; i<numBins; i++)
    hist[1][i] += dot/dBins/pow(val,dim);
    //hist[1][bin] += force/dBins;


  /*
  for(int i=startBin; i<numBins; i++)
  {
    stepVal = step(val,hist[0][i]-dBins/2.0,hist[0][i]+dBins/2.0);
    if(stepVal == 1.0)
    {
      hist[1][i] += force/dBins;
      break;
    }
  }
  */
  return;
}
/*
double forceTrapezoid(double* x, double** y, int n, int dim)
{
	double integrand,ref;

	integrand = 0.0;
	for(int i=1; i<n; i++)
		integrand += ((y[i]+y[i-1])*abs(x[i]-x[i-1]))/2.0;

	return integrand;
}
*/

