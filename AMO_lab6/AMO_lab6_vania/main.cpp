#include "header.h"

int main(){
	int N = 80;
	printf("N = %d\n", N);
	Spline(N);
	printf("Spline built, results in file 'res.txt' \n");

	return 0;
}