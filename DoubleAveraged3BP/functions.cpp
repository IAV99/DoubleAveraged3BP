#include "Header.h"

vect_2d r1(double a, double omega_d, double E1, double e1) {
	vect_2d r;
	r.x = a * (cos(E1) - e1);
	r.y = a * sqrt(1 - e1 * e1) * sin(E1);
	vect_2d r1;
	r1.x = cos(omega_d) * r.x - sin(omega_d) * r.y;
	r1.y = sin(omega_d) * r.x + cos(omega_d) * r.y;
	return r1;
}

vect_2d dr1_domegad(double a, double omega_d, double E1, double e1) {
	vect_2d r;
	r.x = a * (cos(E1) - e1);
	r.y = a * sqrt(1 - e1 * e1) * sin(E1);
	vect_2d r1;
	r1.x = -sin(omega_d) * r.x - cos(omega_d) * r.y;
	r1.y = cos(omega_d) * r.x - sin(omega_d) * r.y;
	return r1;
}

vect_2d dr1_de1(double a, double omega_d, double E1, double e1) {
	vect_2d r;
	r.x = a * (-1);
	r.y = a * (-e1) / sqrt(1 - e1 * e1) * sin(E1);
	vect_2d r1;
	r1.x = cos(omega_d) * r.x - sin(omega_d) * r.y;
	r1.y = sin(omega_d) * r.x + cos(omega_d) * r.y;
	return r1;
}


vect_2d r2(double E2, double e2) {
	vect_2d r;
	r.x = (cos(E2) - e2);
	r.y = sqrt(1 - e2 * e2) * sin(E2);
	return r;
}

vect_2d dr2_de2(double E2, double e2) {
	vect_2d r;
	r.x = -1;
	r.y = - e2 / sqrt(1 - e2 * e2) * sin(E2);
	return r;
}

double W0(double a, double omega_d, double e1, double e2, double E1, double E2) {
	vect_2d r_1 = r1(a, omega_d, E1, e1);
	vect_2d r_2 = r2(E2, e2);
	double mod_dr = sqrt((r_1.x - r_2.x) * (r_1.x - r_2.x) + (r_1.y - r_2.y) * (r_1.y - r_2.y));
	return (1 - e1 * cos(E1)) * (1 - e2 * cos(E2)) / mod_dr;
}

double dW_de1(double a, double omega_d, double e1, double e2, double E1, double E2) {
	vect_2d r_1 = r1(a, omega_d, E1, e1);
	vect_2d dr = dr1_de1(a, omega_d, E1, e1);
	vect_2d r_2 = r2(E2, e2);
	double mod_dr = sqrt((r_1.x - r_2.x) * (r_1.x - r_2.x) + (r_1.y - r_2.y) * (r_1.y - r_2.y));
	double first = - ((r_1.x - r_2.x) * dr.x + (r_1.y - r_2.y) * dr.y) * (1 - e1 * cos(E1)) * (1 - e2 * cos(E2));
	double second = - cos(E1) * (1 - e2 * cos(E2));
 	return first / (mod_dr * mod_dr * mod_dr) + second / mod_dr;
}

double dW_de2(double a, double omega_d, double e1, double e2, double E1, double E2) {
	vect_2d r_1 = r1(a, omega_d, E1, e1);
	vect_2d r_2 = r2(E2, e2);
	vect_2d dr = dr2_de2(E2, e2);
	double mod_dr = sqrt((r_1.x - r_2.x) * (r_1.x - r_2.x) + (r_1.y - r_2.y) * (r_1.y - r_2.y));
	double first = ((r_1.x - r_2.x) * dr.x + (r_1.y - r_2.y) * dr.y) * (1 - e1 * cos(E1)) * (1 - e2 * cos(E2));
	double second = - cos(E2) * (1 - e1 * cos(E1));
	return first / (mod_dr * mod_dr * mod_dr) + second / mod_dr;
}

double dW_domegad(double a, double omega_d, double e1, double e2, double E1, double E2) {
	vect_2d r_1 = r1(a, omega_d, E1, e1);
	vect_2d dr = dr1_domegad(a, omega_d, E1, e1);
	vect_2d r_2 = r2(E2, e2);
	double mod_dr = sqrt((r_1.x - r_2.x) * (r_1.x - r_2.x) + (r_1.y - r_2.y) * (r_1.y - r_2.y));
	double first = -((r_1.x - r_2.x) * dr.x + (r_1.y - r_2.y) * dr.y) * (1 - e1 * cos(E1)) * (1 - e2 * cos(E2));
	return first / (mod_dr * mod_dr * mod_dr);
}

/*
// dW/domega_d
vector <vector <double>> dWdo;
for (int i = 0; i < n; i++) {
    vector <double> w;
    for (int j = 0; j < m; j++) {
        double w_0 = dW_domegad(a, omega_d, e_1, e_2, E1[i], E2[j]);
        w.push_back(w_0);
    }
    dWdo.push_back(w);
}
double dWdOmega = mu_1 * mu_2 / (4 * pi * pi) * int2d(dWdo, E1, E2);

// dW/de1
vector <vector <double>> dWde10;
for (int i = 0; i < n; i++) {
    vector <double> w;
    for (int j = 0; j < m; j++) {
        double w_0 = dW_de1(a, omega_d, e_1, e_2, E1[i], E2[j]);
        w.push_back(w_0);
    }
    dWde10.push_back(w);
}
double dWde1 = mu_1 * mu_2 / (4 * pi * pi) * int2d(dWde10, E1, E2);

// dW/de2
vector <vector <double>> dWde20;
for (int i = 0; i < n; i++) {
    vector <double> w;
    for (int j = 0; j < m; j++) {
        double w_0 = dW_de2(a, omega_d, e_1, e_2, E1[i], E2[j]);
        w.push_back(w_0);
    }
    dWde20.push_back(w);
}
double dWde2 = mu_1 * mu_2 / (4 * pi * pi) * int2d(dWde20, E1, E2);

// W

vector <vector <double>> W;
for (int i = 0; i < n; i++) {
    vector <double> w;
    for (int j = 0; j < m; j++) {
        double w_0 = W0(a, omega_d, e_1, e_2, E1[i], E2[j]);
        w.push_back(w_0);
    }
    W.push_back(w);
}
double W_avg = mu_1 * mu_2 / (4 * pi * pi) * int2d(W, E1, E2);
cout << setprecision(12) << W_avg << endl;


// dW/dOmega
double delta = 0.00001;
double dW_avg[2] = {0, 0};
for (int k = 0; k < 2; k++) {
    vector <vector <double>> W;
    for (int i = 0; i < n; i++) {
        vector <double> w;
        for (int j = 0; j < m; j++) {
            double w_0 = W0(a, omega_d + delta * (2 * (double)k - 1), e_1, e_2, E1[i], E2[j]);
            w.push_back(w_0);
        }
        W.push_back(w);
    }
    dW_avg[k] = mu_1 * mu_2 / (4 * pi * pi) * int2d(W, E1, E2);
}
double dW_avg_r = (dW_avg[1] - dW_avg[0]) / (2 * delta);
cout << endl << "d(W)/d(omega_d)" << endl;
cout << setprecision(12) << dW_avg_r << endl << dWdOmega << endl << "relative difference = ";
cout << setprecision(3) << (dW_avg_r - dWdOmega) / dWdOmega << endl;

// dW/de1
double delta_e1 = 0.00001;
//double dW_avg[2] = { 0, 0 };
for (int k = 0; k < 2; k++) {
    vector <vector <double>> W;
    for (int i = 0; i < n; i++) {
        vector <double> w;
        for (int j = 0; j < m; j++) {
            double w_0 = W0(a, omega_d, e_1 + delta_e1 * (2 * (double)k - 1), e_2, E1[i], E2[j]);
            w.push_back(w_0);
        }
        W.push_back(w);
    }
    dW_avg[k] = mu_1 * mu_2 / (4 * pi * pi) * int2d(W, E1, E2);
}
double dW_avg_e1 = (dW_avg[1] - dW_avg[0]) / (2 * delta_e1);
cout << endl << "d(W)/d(e1)" << endl;
cout << setprecision(12) << dW_avg_e1 << endl << dWde1 << endl << "relative difference = ";
cout << setprecision(3) << (dW_avg_e1 - dWde1) / dWde1 << endl;

// dW/de2
double delta_e2 = 0.00001;
//double dW_avg[2] = { 0, 0 };
for (int k = 0; k < 2; k++) {
    vector <vector <double>> W;
    for (int i = 0; i < n; i++) {
        vector <double> w;
        for (int j = 0; j < m; j++) {
            double w_0 = W0(a, omega_d, e_1, e_2 + delta_e2 * (2 * (double)k - 1), E1[i], E2[j]);
            w.push_back(w_0);
        }
        W.push_back(w);
    }
    dW_avg[k] = mu_1 * mu_2 / (4 * pi * pi) * int2d(W, E1, E2);
}
double dW_avg_e2 = (dW_avg[1] - dW_avg[0]) / (2 * delta_e2);
cout << endl << "d(W)/d(e2)" << endl;
cout << setprecision(12) << dW_avg_e2 << endl << dWde2 << endl << "relative difference = ";
cout << setprecision(3) << (dW_avg_e2 - dWde2) / dWde2 << endl;

*/