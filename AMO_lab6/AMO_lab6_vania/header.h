#ifndef HEADER_H
#define HEADER_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

double function(double);
double *Seidel(double **, int, double);
double* vector_A(int);
double **Matrix_C(int);
double* vector_D(double*, int);
double* vector_B(double*, double*, int);
void Spline(int);

#endif //header