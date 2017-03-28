#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double P = M_PI;

double myfunc(double x) {
	return 0.1*x*x - sin(x - P / 4) - 1.5;
}

double first_myfunc(double x) {
	return 0.2*x - cos(x - P / 4);
}

double second_myfunc(double x) {
	return 0.2 + sin(x - P / 4);
}

int find_roots(double mezha1, double mezha2, double step, double *mas) {

	int i = -1;
	double a;

	for (a = mezha1; a <= mezha2; a += step)
		if (myfunc(a) * myfunc(a + step) < 0)
			mas[++i] = a;
	return i;
}


double* method_iteration(double a, double b, double eps) {
	double m, M;
	double alfa, q;
	double x1, x;
	int n;
	double *res= (double*)malloc(sizeof(double) * 3);

	if (fabs(first_myfunc(a)) > fabs(first_myfunc(b))) {
		M = first_myfunc(a);
		m = first_myfunc(b);
	}
	else {
		M = first_myfunc(b);
		m = first_myfunc(a);
	}
	alfa = 1 / M;
	q = 1 - fabs(m / M);

	x = (b + a) / 2;
	x1 = x;
	n = 0;
	do {
		n++;
		x = x1;
		x1 = x - alfa*myfunc(x);
	} while (fabs(x1 - x) > (1 - q) / q*eps);
	res[0] = x1;
	res[1] = fabs(fabs(x1) - fabs(x))*q / (1 - q);
	res[2] = (double)n;


	return res;
}

double*  method_hord(double a, double b, double eps) {
	double c, m, x, x1;
	int n;
	double *res = (double*)malloc(sizeof(double) * 3);

	m = fabs(first_myfunc(a));
	for (x1 = a; x1 < b; x1 += 0.001) {
		if (fabs(first_myfunc(x1)) < m) 
			m = fabs(first_myfunc(x1));
	}

	if (myfunc(a) * second_myfunc(a) > 0){
		c = a;
		x1 = b;
	}
	else{
		c = b;
		x1 = a;
	}

	x = x1;
	n = 0;
	while (fabs(myfunc(x)) / m > eps) {
		n++;
		x -= myfunc(x)*(x - c) / (myfunc(x) - myfunc(c));
	}

	res[0] = x;
	res[1] = fabs(myfunc(x)) / m;
	res[2] = (double)n;

	return res;
}

main() {
	double mas[4];
	double *res1 = (double*)malloc(sizeof(double) * 3);
	double *res2 = (double*)malloc(sizeof(double) * 3);
	int i, k;
	double eps, step = 0.3;
	double a, b;

	k = find_roots(-10, 10, step, mas);

	printf("ITERATION:\n\n");
	for (i = 0; i <= k; i++){
		a = mas[i];
		b = a + step;
		printf(" EPS            ROOT_%d              ACCURACY\n\n", i + 1);
		for (eps=1e-2; eps >= 1e-14; eps*= 1e-3){	
			res1 = method_iteration(a, b, eps);
			printf("%0.E  %20.15f  %20.15f\n", eps, res1[0], res1[1]);
		}
		printf("\n");
	}

	printf("-----------------------------------------------------------\n");
	printf("HORD:\n\n");
	for (i = 0; i <= k; i++){
		a = mas[i];
		b = a + step;
		printf(" EPS            ROOT_%d              ACCURACY\n\n", i + 1);
		for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3){
			res2= method_hord(a, b, eps);
			printf("%0.E  %20.15f  %20.15f\n", eps, res2[0], res2[1]);
		}
		printf("\n");
	}

	printf("-----------------------------------------------------------\n");
	printf(" EPS      ITERATION      HORD \n\n");
	a = mas[0];
	b = a + step;
	for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3){
		res1 = method_iteration(a, b, eps);
		res2 = method_hord(a, b, eps);
		printf("%0.E   %6.0d     %6.0d\n", eps, (int)res1[2], (int)res2[2]);
	}


	system("PAUSE");

	return;
}

