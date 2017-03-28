#include "header.h"

double  a = 2;
double b = 9;

double function(double x){
	return 0.2*pow(M_E, pow(x, (1.0/3)))*log(x)*sin(3 * x);
}

double *Seidel(double **matrix, int N, double eps){
	double *PrevEquation = new double[N];
	double *SRoot = new double[N];
	double tmp = 0, q = 0, max = 0, M_norma;
	int k;

	for (int i = 0; i < N; i++){
		tmp = matrix[i][i];
		for (int j = 0; j <N + 1; j++)
			matrix[i][j] = matrix[i][j] / tmp;
	}
	for (int i = 0; i < N; i++)
		SRoot[i] = matrix[i][N];

	q = 0;
	for (int i = 0; i < N; i++){
		max = 0;
		for (int j = 0; j < N; j++)
			if (i != j)
				max += fabs(matrix[i][j]);
		if (max > q)
			q = max;
	}
	k = 1;
	do{
		for (int i = 0; i < N; i++)
			PrevEquation[i] = SRoot[i];
		for (int i = 0; i < N; i++){
			SRoot[i] = matrix[i][N];
			for (int j = 0; j < N; j++)
				if (i != j)
					SRoot[i] -= matrix[i][j] * SRoot[j];
		}

		for (int i = 0; i < N; i++)
			PrevEquation[i] = SRoot[i] - PrevEquation[i];

		M_norma = fabs(PrevEquation[0]);
		for (int i = 1; i < N; i++)
			if (fabs(PrevEquation[i])>M_norma)
				M_norma = fabs(PrevEquation[i]);

		for (int i = 0; i < N; i++)
			PrevEquation[i] = SRoot[i];
		k++;
	} while (M_norma> eps*((1 - q) / q));

	return SRoot;
}

double* vector_A(int N){
	double* A = new double[N];
	double x = a, h = (b - a) / N;
	for (int i = 0; i < N; i++)	{
		A[i] = function(x);
		x += h;
	}
	return A;
}

double **Matrix_C(int N){
	double **matrix = new double*[N];
	for (int i = 0; i < N; i++)
		matrix[i] = new double[N+1];

	for (int i = 0; i < N ; i++)
		for (int j = 0; j < N + 1; j++)
			matrix[i][j] = 0;

	double h = (b - a) / N;

	for (int i = 0; i < N; i++){
		if (i>0)
			matrix[i][i - 1] = h;
		if (i < N-1)
			matrix[i][i + 1] = h;
		matrix[i][i] = 4 * h;
		matrix[i][N] = 6 * (function(a + h * (i - 1)) - 2 * function(a + h * i) + function(a + h * (i + 1))) / h;
	}
	return matrix;
}

double* vector_D(double* C, int N){
	double* D = new double[N];
	double h = (b - a) / N;

	D[0] = C[0];
	for (int i = 1; i < N; ++i)
		D[i] = (C[i] - C[i - 1]) / h;
	return D;
}

double* vector_B(double* C, double* D, int N){
	double* B = new double[N];
	double h = (b - a) / N;
	B[0] = 0;

	for (int i = 1; i < N; ++i)
		B[i] = (h * C[i] / 2) - (h * h * D[i] / 2) + (function(a + h*i) - function(a + h*(i - 1))) / h;
	return B;
}

void Spline(int N){

	double* A;
	double* C;
	double* D;
	double* B;
	double **matrix = new double*[N];
	for (int i = 0; i < N; i++)
		matrix[i] = new double[N + 1];

	A = vector_A(N);
	matrix = Matrix_C(N);
	C = Seidel(matrix, N, 1e-2);
	D = vector_D(C, N);
	B = vector_B(C, D, N);

	double h = (b - a) / N;
	int k = (b - a) * 10;
	double step = (b - a) /(double)k;
	double x = a, y, xj;

	ofstream fout("res.txt");
	for (int i = 0; i <= k; i++){
		int j = 0;
		while (x >= (a + j * h))
			j++;
		j--;
		xj = a + h * j;
		y = A[j] + B[j] * (x - xj) + C[j] * (x - xj) * (x - xj) / 2 + D[j] * (x - xj) * (x - xj) * (x - xj) / 6;
		fout << x;
		fout << " ; ";
		fout << y;
		fout << endl;
		x += step;
	}
	fout.close();
	return;
}