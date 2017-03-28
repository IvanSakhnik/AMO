#include "header.h"

void Gauss(double Matrix[n][n + 1], double GRoot[n]){
	double main1 = 0, main2 = 0, main3 = 0, helper[n + 1];
	//під головною діагоналлю "0", в головній діагоналі "1"
	for (int i = 0; i < n; i++){
		if (Matrix[i][i] != 0){
			main1 = Matrix[i][i];
			for (int j = 0; j <= n; j++){
				Matrix[i][j] /= main1;
			}
			for (int k = i; k < n - 1; k++){
				main2 = Matrix[k + 1][i];
				for (int j = 0; j <= n; j++){
					helper[j] = Matrix[i][j] * main2;
					Matrix[k + 1][j] -= helper[j];
				}
			}
		}
		else{
			printf("Divide by Zero, Coefficient %d %d ", i, i);
			printf(" = Zero; \n");
		}
	}
	//над головною діагоналлю "0"
	for (int i = n - 1; i >= 0; i--){
		for (int ki = i; ki > 0; ki--){
			main3 = Matrix[ki - 1][i];
			for (int j = 0; j <= n; j++){
				helper[j] = Matrix[i][j] * main3;
				Matrix[ki - 1][j] -= helper[j];
			}
		}
	}
	//корені рівняння
	for (int i = 0; i < n; i++){
		GRoot[i] = Matrix[i][n];
	}
}

void Seidel(double Modern_matrix[n][n + 1], double SRoot[n], double &eps){
	double PrevEquation[n], tmp = 0, q = 0, max = 0, M_norma;
	int k;

	for (int i = 0; i < n; i++){
		tmp = Modern_matrix[i][i];
		for (int j = 0; j <n+1; j++)
			Modern_matrix[i][j] = Modern_matrix[i][j] / tmp;
	}
	for (int i = 0; i < n; i++)
		SRoot[i] = Modern_matrix[i][n];

	q = 0;
	for (int i = 0; i < n; i++){
		max = 0;
		for (int j = 0; j < n; j++)
			if (i != j)
				max += fabs(Modern_matrix[i][j]);
		if (max > q)
			q = max;
	}
	k = 1;
	do{
		for (int i = 0; i < n; i++)
			PrevEquation[i] = SRoot[i];
		for (int i = 0; i < n; i++){
			SRoot[i] = Modern_matrix[i][n];
			for (int j = 0; j < n; j++)
				if (i != j)
					SRoot[i] -= Modern_matrix[i][j] * SRoot[j];
		}

		for (int i = 0; i < n; i++)
			PrevEquation[i] = SRoot[i] - PrevEquation[i];

		M_norma = fabs(PrevEquation[0]);
		for (int i = 1; i < n; i++)
			if (fabs(PrevEquation[i])>M_norma)
				M_norma = fabs(PrevEquation[i]);

		for (int i = 0; i < n; i++)
			PrevEquation[i] = SRoot[i];
		k++;
	} while (M_norma> eps*((1 - q) / q));
	printf("Iteration= %d\n", k);
}
