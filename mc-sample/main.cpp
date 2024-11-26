
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <mpi.h>
#include <filesystem>

#include "system.h"
#include "functions.h"
#include "utility.h"

using namespace std;

int main()
{
  clock_t timer;
  timer = clock();
	bool first,thermostat,constrain;
	int		numPart,dim,kmax,nxi;
	double spread;
  string holdS;

//	srand(1989);

	double g,k=1.0,tol,boltz,mass,dxi,tConv,dConv;
	double kBead,rEq,potPar[2],tstep,tauT;
	double mTot,energy,pot,kinetic,dist, energy0, temp, beta, sigma, tscale;
	double engAcc, eng2Acc, engAvg, engStd, energyCon;
	double *m,*mk,**q,**p,***a,***pk,***qxi;
  double **force,***forcek,**potxi;
  double ***posHist, ***momHist,*posHistNorm,*momHistNorm;

	double pi = 4.0*atan(1.0);
	int maxiter;
	int maxstep;
	int numSample,numCor,numFreq,corStep;
	double conStart,conEnd,delCon,cortStep;
	int printStep,progStep,statStep;
	int* constraints;
	int numCon,numTraj,numConfig,numBins,deCor;
  double rStart, rEnd,dBins;

  int potType;
 
  int idnode,mxnode;

  MPI_Init(NULL,NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &idnode);
  MPI_Comm_size(MPI_COMM_WORLD, &mxnode);

  srand(1989);
	//srand(time(NULL)+idnode);

  string potName, inFile="sample.txt";
	ifstream ip;
	ofstream fRDF;
  ofstream fConfig;
  ip.open(inFile,ios::in);
  timer = clock_t();
	if(ip.fail()) {cout << "Can't open file: "<< inFile << endl; exit(0); }

  ip >> numPart >> dim;  
	ip >> kmax >> nxi;
  ip >> temp;
	ip >> maxstep;
  ip >> deCor;
  ip >> tstep;
  ip >> numBins;
  ip >> rStart >> rEnd;
  ip >> rEq;
  ip >> potName;
  ip.close();

  if( potName.compare("harmonic") == 0)
  {
    potType = 1;
    mass = 1.0;
    boltz = 1.0;
    beta = temp;
    cout << "Potential Type: Harmonic" << endl;
  }
  else if( potName.compare("mildanharmonic") == 0)
  {
    potType = 2;
    mass = 1.0;
    boltz = 1.0;
    beta = temp;
    cout << "Potential Type: Mildly Anharmonic" << endl;
  }
  else if( potName.compare("quartic") == 0)
  {
    potType = 3;
    mass = 1.0;
    boltz = 1.0;
    beta = temp;
    cout << "Potential Type: Quartic" << endl;
  }
  else if( potName.compare("ohmodel") == 0)
  {
    potType = 4;
    mass = 1741.1;
    boltz = 3.1668115363e-6;
    beta = 1.0/(temp*boltz);
    cout << "Potential Type: OH-model" << endl;
  }
  else if( potName.compare("qtip4pf") == 0)
  {
    potType = 5;
    mass = 1.0;
    boltz = 1.0;
    beta = temp;
    cout << "Potential Type: Mildly Anharmonic" << endl;
  }
  else
  {
    cout << "Unrecognized potential type" << endl;
    exit(1);
  }

  if(dim > 2)
  {
    cout << "Dimensions greater than 2" << endl;
    exit(1);
  }
  
  //numPart = 1;
	m = new double[numPart];
	mk = new double[kmax];
  //rad = new double[numPart];
	//theta = new double[numPart];

	//dif = new double[dim];
	//com = new double[dim];

	q = new double*[numPart];
	p = new double*[numPart];
  qxi = new double**[numPart];
  a = new double**[numPart];
  pk = new double**[numPart];
	force = new double*[numPart];
	forcek = new double**[numPart];
  potxi = new double*[numPart];

	for(int i=0; i<numPart; i++)
	{
    q[i] = new double[dim];  
    p[i] = new double[dim];  
	  force[i] = new double[dim];
    a[i] = new double*[kmax];
    pk[i] = new double*[kmax];
	  forcek[i] = new double*[kmax];
    potxi[i] = new double[nxi];
    qxi[i] = new double*[dim];
    for(int j=0; j<kmax; j++)
    {
      a[i][j] = new double[dim];
      pk[i][j] = new double[dim];
      forcek[i][j] = new double[dim];
    }
    for(int j=0; j<dim; j++)
      qxi[i][j] = new double[nxi];
	}

  tConv = 2.4188843265857e-2; //convert a.u. to fs
  dConv =  0.529177210903;    //convert a.u. to angstrom

  if(maxstep < 1000) printStep=1;
  else printStep=maxstep/100;

  //statStep = 1000;
  statStep = 0;
  thermostat=true;

  writeOutDbl(beta,"beta");
  writeOutDbl(mass,"mass");
  writeOutInt(dim,"dim");
  writeOutInt(numPart,"numPart");
  writeOutInt(kmax,"kmax");
  writeOutInt(nxi,"nxi");
  mass = mass /((double) numPart);
  sigma = ((double) numPart * (double) dim)/(2.0*beta);
	spread = pow(1e-3,sqrt(beta/mass))/numPart;
  kBead = mass * ((double)numPart*(double)numPart)/(beta*beta);
  //rEq = 0.0;
  potPar[0] = 0.18748;
  potPar[1] = 1.1605;
  tauT = 8.0e1;
  dxi = 1.0/((double)nxi-1);

  //nmSetup(nmFreq,nmTrans,beta,numPart,dim);
  //writeOutMatrix(nmTrans,"nmTrans",numPart,numPart);

  mTot = 0.0;
  for(int i=0;i<numPart;i++)
  {
    m[i] = mass;
    mTot += m[i];
  }
  for(int k=0;k<kmax;k++)
    mk[k] = 0.5/((double)(k+1));

  writeOutArray(m,"m",numPart);
  writeOutDbl(kBead,"kbead");
  writeOutDbl(dxi,"dxi");
  /*
  for(int i=0; i<numBins; i++)
  {
    radHist[0][i] = (double)(i+0.5)*dBins + rStart;
    radHist[1][i] = 0.0;
    radHist0[i] = 0.0;
    forceHist[0][i] = (double)(i+0.5)*dBins + rStart;
    forceHist[1][i] = 0.0;
    forceHist0[i] = 0.0;
  }
*/
  for(int i=0; i<numPart; i++)
  {  
    for(int j=0; j<dim; j++)
    {
      q[i][j] = spread*(2.0*randGen() - 1.0);
      p[i][j] = randGaus(sqrt(m[i]/beta),0.0);
      for(int k=0; k<kmax; k++)
      {
        a[i][k][j] = spread*(2.0*randGen() - 1.0);
        pk[i][k][j] = randGaus(sqrt(mk[k]/beta),0.0);
      }
    }
  }

  //writeOutMatrix(q,"q",numPart,dim);
  //getCOM(q,com,m,mTot,numPart,dim);
  //writeOutArray(com,"com",dim);

  for(int i=0; i<numPart; i++)
  {
    for(int j=0; j<dim; j++)
    {
      q[i][j] = q[i][j] + randGaus(0.05,rEq)/sqrt(2.0);///sqrt(2);
    }
    //q[i][1] -= com[1] + sqrt(rEq*rEq - q[i][0]*q[i][0]);//-rCon/sqrt(2);
  }
  writeOutMatrix(q,"q",numPart,dim);
  writeOutTensor(a,"a",numPart,kmax,dim);
  //getCOM(q,com,m,mTot,numPart,dim);
  //writeOutArray(com,"com",dim);
  //getRadius(q,rad,numPart,dim);
  //writeOutArray(rad,"rad",dim);
  
  pot = 0.0;
  if(numPart > 1) pot += getBeadPot(q, kBead, numPart, dim);
  writeOutDbl(pot,"pot");
  pot += getFourierPot(a, kBead, numPart, kmax, dim);
  writeOutDbl(pot,"pot");
  pot += getXiPot(q,a,qxi,potxi,rEq,dxi,numPart,kmax,dim,nxi,potType);
  writeOutDbl(pot,"pot");
  //pot += getPot(rad, rEq, numPart,potType)/((double)numPart); 
  

  writeOutMatrix(q,"q",numPart,dim);
  
/*  for(int i=0; i<numPart; i++)
  {
    writeOutInt(i,"bead");
    writeOutMatrix(a[i],"ajk",kmax,dim);
  }
*/
  fConfig.open("beadconfig.xvg",ios::out);
  fConfig << "Configs" << endl;
  fConfig.close();
  fConfig.open("fourierconfig.xvg",ios::out);
  fConfig << "Configs" << endl;
  fConfig.close();

  int dipCount=0;
 
  //double step = tstep;
  double *stepSize,*accepRatio;
  int *attempts,*success;
  //bool stepGood = false;
  int goodSteps; 

  int nAccep,testStep,totDim;
  int part,partk,mover;
  double *pos0, pot0,rad0;
  double  accepLow, accepHigh, dStep;

  testStep = 50000;
  accepLow = 0.48;
  accepHigh = 0.52;
  dStep = 0.01;

  goodSteps=0;
  totDim = kmax+1;
  pos0 = new double[dim];
  stepSize = new double[totDim];
  success = new int[totDim];
  attempts = new int[totDim];
  accepRatio = new double[totDim];
  posHist = new double**[dim*totDim];
  momHist = new double**[dim*totDim];
  posHistNorm = new double[dim*totDim];
  momHistNorm = new double[dim*totDim];

  for(int i=0; i<totDim; i++)
  {
    stepSize[i] = tstep;
    attempts[i] = 0;
    success[i] = 0;
  } 
  for(int i=0; i<dim*totDim; i++)
  {
    posHist[i] = new double*[2];
    momHist[i] = new double*[2];
    for(int j=0; j<2; j++)
    {
      posHist[i][j] = new double[numBins];
      momHist[i][j] = new double[numBins];
    }
  }

  dBins = (rEnd-rStart)/(double)numBins;
  for(int i=0; i<numBins; i++)
  {
    posHist[0][0][i] = (double)(i+0.5)*dBins + rStart;
    momHist[0][0][i] = (double)(i+0.5)*dBins + rStart;
  }
  while(goodSteps!=totDim)
  {
    for(int i=0; i<totDim; i++)
    {
      attempts[i] = 0;
      success[i] = 0;
    } 
    for(int i=0; i<testStep; i++)
    {
      if(randGen() < 0.5) //move bead
      {
        for(int j=0; j<numPart; j++)
        {
          part = (int)(numPart*randGen());
          if(part >= numPart)
            part = numPart-1;
          
          for(int k=0; k<dim; k++)
            pos0[k] = q[part][k];
          pot0 = pot;
          //rad0 = rad[part];

          for(int k=0; k<dim; k++)
            q[part][k] = q[part][k] + stepSize[0]*(randGen() - 0.5);
          
          //getRadius(q,rad,numPart,dim);
          pot = 0.0;
          if(numPart > 1) pot += getBeadPot(q, kBead, numPart, dim);
          pot += getFourierPot(a, kBead, numPart, kmax, dim);
          pot += getXiPot(q,a,qxi,potxi,rEq,dxi,numPart,kmax,dim,nxi,potType);
          //pot += getPot(rad, rEq, numPart,potType)/((double)numPart); 
   
          attempts[0]++;
          if(randGen() < exp(-beta*(pot-pot0)))
          {
            success[0]++;
          }
          else
          {
            pot = pot0;
            //rad[part] = rad0;
            for(int k=0; k<dim; k++)
              q[part][k] = pos0[k];
          }
        }
      }
      else //move fourier amplitude
      {
        partk = (int)(kmax*randGen()); // choose fourier component
        mover = partk+1;
        for(int j=0; j<numPart; j++)
        {
          part = (int)(numPart*randGen());
          if(part >= numPart)
            part = numPart-1;
          
          for(int k=0; k<dim; k++)
            pos0[k] = a[part][partk][k];
          pot0 = pot;
          //rad0 = rad[part];

          for(int k=0; k<dim; k++)
            a[part][partk][k] = a[part][partk][k] + stepSize[mover]*(randGen() - 0.5);
          
          //getRadius(q,rad,numPart,dim);
          pot = 0.0;
          if(numPart > 1) pot += getBeadPot(q, kBead, numPart, dim);
          pot += getFourierPot(a, kBead, numPart, kmax, dim);
          pot += getXiPot(q,a,qxi,potxi,rEq,dxi,numPart,kmax,dim,nxi,potType);
          //pot += getPot(rad, rEq, numPart,potType)/((double)numPart); 
   
          attempts[mover]++;
          if(randGen() < exp(-beta*(pot-pot0)))
          {
            success[mover]++;
          }
          else
          {
            pot = pot0;
            //rad[part] = rad0;
            for(int k=0; k<dim; k++)
              a[part][partk][k] = pos0[k];
          }
        }
      }
    }

    goodSteps=0;
    for(int i=0; i<totDim; i++)
    {
      accepRatio[i] = (double)success[i]/(double)attempts[i];
      
      if(accepRatio[i] < accepLow)
      {
        stepSize[i] -= dStep;
      }
      else if (accepRatio[i] > accepHigh)
      {
        stepSize[i] += dStep;
      }
      else if (accepRatio[i] > accepLow && accepRatio[i] < accepHigh)
        goodSteps++;
    }
    writeOutArray(accepRatio,"accepRatio",totDim);
    writeOutArray(stepSize,"stepSize",totDim);
    
    /*
    accepRatio = (double)nAccep/(double)(testStep*numPart);
    writeOutDbl(accepRatio,"accepRatio");

    if(accepRatio < accepLow)
    {
      step = step - dStep;
      writeOutDbl(step,"step");
    }
    else if(accepRatio > accepHigh)
    {
      step = step + dStep;
      writeOutDbl(step,"step");
    }
    else if(accepRatio > accepLow && accepRatio < accepHigh)
      stepGood = true;
  */
  }

    
  //writeOutDbl(step,"step");
  cout << "Starting config generation" << endl;
  nAccep = 0;
  for(int iter=1;iter<=maxstep+statStep;iter++)
  {
    for(int i=0; i<totDim; i++)
    {
      attempts[i] = 0;
      success[i] = 0;
    } 
    for(int i=0; i<deCor; i++)
    {
      if(randGen() < 0.5) //move bead
      {
        for(int j=0; j<numPart; j++)
        {
          part = (int)(numPart*randGen());
          if(part >= numPart)
            part = numPart-1;
          
          for(int k=0; k<dim; k++)
            pos0[k] = q[part][k];
          pot0 = pot;
          //rad0 = rad[part];

          for(int k=0; k<dim; k++)
          {
            q[part][k] = q[part][k] + stepSize[0]*(randGen() - 0.5);
            p[part][k] = randGaus(sqrt(m[part]/beta),0.0);
          }

          //getRadius(q,rad,numPart,dim);
          pot = 0.0;
          if(numPart > 1) pot += getBeadPot(q, kBead, numPart, dim);
          pot += getFourierPot(a, kBead, numPart, kmax, dim);
          pot += getXiPot(q,a,qxi,potxi,rEq,dxi,numPart,kmax,dim,nxi,potType);
          //pot += getPot(rad, rEq, numPart,potType)/((double)numPart); 
   
          attempts[0]++;
          if(randGen() < exp(-beta*(pot-pot0)))
          {
            success[0]++;
          }
          else
          {
            pot = pot0;
            //rad[part] = rad0;
            for(int k=0; k<dim; k++)
              q[part][k] = pos0[k];
          }
        }
      }
      else //move fourier amplitude
      {
        partk = (int)(kmax*randGen()); // choose fourier component
        mover = partk+1;
        for(int j=0; j<numPart; j++)
        {
          part = (int)(numPart*randGen());
          if(part >= numPart)
            part = numPart-1;
          
          for(int k=0; k<dim; k++)
            pos0[k] = a[part][partk][k];
          pot0 = pot;
          //rad0 = rad[part];

          for(int k=0; k<dim; k++)
          {
            a[part][partk][k] = a[part][partk][k] + stepSize[mover]*(randGen() - 0.5);
            pk[part][partk][k] = randGaus(sqrt(mk[k]/beta),0.0);
          }

          //getRadius(q,rad,numPart,dim);
          pot = 0.0;
          if(numPart > 1) pot += getBeadPot(q, kBead, numPart, dim);
          pot += getFourierPot(a, kBead, numPart, kmax, dim);
          pot += getXiPot(q,a,qxi,potxi,rEq,dxi,numPart,kmax,dim,nxi,potType);
          //pot += getPot(rad, rEq, numPart,potType)/((double)numPart); 
   
          attempts[mover]++;
          if(randGen() < exp(-beta*(pot-pot0)))
          {
            success[mover]++;
          }
          else
          {
            pot = pot0;
            //rad[part] = rad0;
            for(int k=0; k<dim; k++)
              a[part][partk][k] = pos0[k];
          }
        }
      }
    }

    if(iter > statStep)
    { 
        //histogram(qCent[0],radHist,rStart,rEnd,numBins);
        //forceHistogram(qCent[0],forceBAvg,forceHist,rStart,rEnd,numBins,0);
      for(int i=0; i<numPart; i++)
      {
        for(int j=0; j<totDim; j++)
        {
          for(int k=0; k<dim; k++)
          {
            if(j==0)
            {
              histogram(q[i][k],posHist[k*totDim],rStart,rEnd,numBins);
              histogram(p[i][k],momHist[k*totDim],rStart,rEnd,numBins);
            }
            else
            {
              histogram(a[i][j-1][k],posHist[k*totDim+j],rStart,rEnd,numBins);
              histogram(pk[i][j-1][k],momHist[k*totDim+j],rStart,rEnd,numBins);
            }
          }
        }
      }
    }
   
  // Write out trajectory info		

    //writeOutMatrix(p,"p",numPart,dim);
    //writeOutTensor(pk,"pk",numPart,kmax,dim);
    if(idnode == 0 && iter > statStep)
    {
      fConfig.open("beadconfig.xvg",ios::app);
      fConfig << setw(15) << iter-statStep;
      for(int j=0; j<numPart; j++)
      {
        for(int k=0; k<dim; k++)
        {
          fConfig << setw(15) << q[j][k];
        }
      }
      for(int j=0; j<numPart; j++)
      {
        for(int k=0; k<dim; k++)
        {
          fConfig << setw(15) << p[j][k];
        }
      }
      fConfig << endl;
      fConfig.close();
      fConfig.open("fourierconfig.xvg",ios::app);
      fConfig << setw(15) << iter-statStep;
      for(int j=0; j<numPart; j++)
      {
        for(int k=0; k<kmax; k++)
        {
          for(int l=0; l<dim; l++)
            fConfig << setw(15) << a[j][k][l];
        }
      }
      for(int j=0; j<numPart; j++)
      {
        for(int k=0; k<kmax; k++)
        {
          for(int l=0; l<dim; l++)
            fConfig << setw(15) << pk[j][k][l];
        }
      }
      fConfig << endl;
      fConfig.close();
    }
  }
  
  //}
  //accepRatio = (double)nAccep/(double)(maxstep*numPart*deCor);
  //writeOutDbl(accepRatio,"accepRatio");

  for(int i=0; i<totDim; i++)
    accepRatio[i] = (double)success[i]/(double)attempts[i];
  writeOutArray(accepRatio,"accepRatio",totDim);
  MPI_Barrier(MPI_COMM_WORLD);
/*
  for(int j=0; j<numBins; j++)
  {
    forceHist[1][j] = beta*forceHist[1][j]/((double)(maxstep-statStep));
      //forceHist[0][j]/pow(forceHist[0][j],dim);
      //1.0/pow(forceHist[0][j],dim);
    //`mapHist[1][j] = mapHist[1][j]/((double) numTraj*maxstep);
  }
  */
  /*
  MPI_Reduce(radHist[1],radHist0,numBins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(forceHist[1],forceHist0,numBins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  double normRad, normTheta, normForce, shift;
*/
  if(idnode == 0)
  {
    //for(int i=0; i<numBins; i++)
      //forceRDF[i] = trapezoid(forceHist[0],forceHist0,i);
   //   forceRDF[i] = (forceHist0[i] > 0.0) ? forceHist0[i] : 0.0;

   // normRad = trapezoid(radHist[0],radHist0,numBins);
    //normForce = trapezoid(forceHist[0],forceRDF,numBins);
    for(int i=0; i<dim*totDim; i++)
    {
      posHistNorm[i] = trapezoid(posHist[0][0],posHist[i][1],numBins);
      momHistNorm[i] = trapezoid(posHist[0][0],momHist[i][1],numBins);
    } 

    fRDF.open("rdf-mc.xvg",ios::out); 
    fRDF << "RDF" << endl;
    for(int i=0; i<numBins; i++)
    {
      fRDF << setw(15) << posHist[0][0][i];
      for(int j=0; j<dim*totDim; j++)
        fRDF << setw(15) << posHist[j][1][i]/posHistNorm[j];
        //if(dim > 1)
        //  fHist << setw(15) << thetaHist[0][i]
        //  << setw(15) << thetaHist[1][i]/normTheta;
        fRDF << endl;
    }
    fRDF.close();
  }

  MPI_Barrier(MPI_COMM_WORLD);


  if(idnode==0)
  {
    timer = clock();
    double simTime = (double)timer/CLOCKS_PER_SEC;
    cout << "Calculation took " << simTime << "s" << endl;
  }
    

  MPI_Finalize();
  return 0;
}
