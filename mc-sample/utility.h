
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>

using namespace std;

void writeOutInt(int val, string name);

void writeOutDbl(double val, string name);

void writeOutArray(double* val, string name, int length);

void writeOutMatrix(double** mat, string name, int lengthA, int lengthB);

void writeOutTensor(double*** mat, string name, int lengthA, int lengthB,
    int lengthC);

double randGen();

double randGaus(double sigma, double mu);
