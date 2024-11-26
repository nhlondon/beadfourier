
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <fftw3.h>

#include "utility.h"

using namespace std;

void writeOutInt(int val, string name)
{
	cout << name << ": " << val << endl;
	return;
}

void writeOutDbl(double val, string name)
{
	cout << name << ": " << val << endl;
	return;
}

void writeOutArray(double* array, string name, int length)
{
	cout << name << ": ";
	for(int i=0; i<length; i++)
	{
		cout << array[i] << " ";
	}
	cout << endl;
}

void writeOutMatrix(double** mat, string name, int lengthA, int lengthB)
{

	for(int i=0; i<lengthA; i++)
	{
		cout << name+to_string(i) << ": ";
		for(int j=0; j<lengthB; j++)
		{
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
}

void writeOutTensor(double*** mat, string name, int lengthA, int lengthB,
    int lengthC)
{

	for(int i=0; i<lengthA; i++)
	{
		for(int j=0; j<lengthB; j++)
		{
		cout << name+to_string(i)+to_string(j) << ": ";
    for(int k=0; k<lengthC; k++)
			cout << mat[i][j][k] << " ";
		}
		cout << endl;
	}
}

double randGen()
{
	return ((double)(rand())+1.0)/((double)(RAND_MAX)+1.0);
}

double randGaus(double sigma, double mu)
{
	double v1,v2,x;

	v1 = randGen();
	v2 = randGen();

	x = (cos(2.0*3.14*v2)*sqrt(-2.0*log(v1)))*sigma + mu;

	return x;
}
