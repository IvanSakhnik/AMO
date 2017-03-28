#ifndef HEADER_H
#define HEADER_H

#include "stdio.h"
#include "stdlib.h"
#include "iostream"
#include "cmath"

#define n 4 //Кількість коренів у СЛАР

void Gauss(double Matrix[n][n + 1], double Root[n]);
void Seidel(double Matrix[n][n + 1], double Root[n], double &eps);

#endif //MYHEADER_H