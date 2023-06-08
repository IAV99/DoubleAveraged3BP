#include "Header.h"

double int1d_simps(vector <double> f, vector <double> x) {
	int n = (int)x.size();
	double sum = 0;
	for (int k = 0; k < n/2 ; k++) {
		sum += (x[2*k + 2] - x[2*k]) / 6 * (f[2*k+2] + 4*f[2*k + 1] + f[2*k]);
	}
	return sum;

}

double int2d(vector <vector <double>> f, vector <double> x, vector <double> y) {
	int n = (int)x.size(), m = (int)y.size();
	vector <double> sum;
	for (int i = 0; i < n; i++) {
		sum.push_back(int1d_simps(f[i], y));
	}
	return int1d_simps(sum, x);
}