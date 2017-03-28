/*Lab3_AMO 2.11.2016p.
		made by Ivan Sakhnik KB-42*/

#include "header.h"
int main(){

	double M[n][n + 1] = { { 3, 2, 8, 3, 35 }, 
						 { 17, 51, 13, 20, 407 }, 
						 { 0, 3, 7, 3, 28 }, 
						 { 16, 4, 20, 16, 124 } },
		ModernM[n][n + 1] = { { 3, -1, 1, 0, 7 }, 
							{ 17, 51, 13, 20, 407 }, 
							{ 0, 3, 7, 3, 28 }, 
							{ 1, 3, 1, 10, 36 } },
		GRoot[n], SRoot[n], eps = 1e-6;

	Gauss(M, GRoot);
	printf("Gauss Method: \n");
	for (int i = 0; i < 4; i++){
		printf("x%d = %.12f\n", i, GRoot[i]);
	}
	printf("\n");

	printf("Seidel Method: \n");
	printf("eps = %.13f \n", eps);
	Seidel(ModernM, SRoot, eps);
	
	for (int i = 0; i < 4; i++){
		printf("x%d = %.12f\n", i, SRoot[i]);
	}
	return 0;
}