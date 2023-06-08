#pragma once

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const double pi = 3.14159265358979323846;

struct vect_2d {
	double x = 0, y = 0;
};

struct vect_3d {
	double x = 0, y = 0, z = 0;
};

struct parameters {
	double mu = 0.01, a = 0.5, mu_1 = 0.5, sigma = 0.5;
};

struct dWd {
	vector <double> E1, E2;
	int n = 101, m = 101;
};

// from params.cpp
double Sigma(double a, double mu_1, double e1, double e2);
double e1(double a, double mu_1, double sigma, double e2);
double e2(double a, double mu_1, double sigma, double e1);
double e1_min(double a, double mu_1, double sigma);
double e2_min(double a, double mu_1, double sigma);
double e1_max(double a, double mu_1, double sigma);
double e2_max(double a, double mu_1, double sigma);
//from functions.cpp
vect_2d r1(double a, double omega_d, double E1, double e1);
vect_2d dr1_domegad(double a, double omega_d, double E1, double e1);
vect_2d dr1_de1(double a, double omega_d, double E1, double e1);
vect_2d r2(double E2, double e2);
vect_2d dr2_de2(double E2, double e2);
double W0(double a, double omega_d, double e1, double e2, double E1, double E2);
double dW_de1(double a, double omega_d, double e1, double e2, double E1, double E2);
double dW_de2(double a, double omega_d, double e1, double e2, double E1, double E2);
double dW_domegad(double a, double omega_d, double e1, double e2, double E1, double E2);
//from kepler_pr.cpp
double scal_2d(vect_2d r, vect_2d d);
double mod_2d(vect_2d r);
vect_2d r1_0(double a, double e, double omega_d);
vect_2d r2_0(double a, double e);
vect_2d p1_0(double a, double mu, double e, double omega_d);
vect_2d p2_0(double a, double mu, double e, double omega_d);
double a_from_pr(vect_2d r, vect_2d p, double mu);
double e_from_pr(vect_2d r, vect_2d p, double mu);
double true_anomaly_from_pr(vect_2d r, vect_2d p, double mu);
double omega_from_pr(vect_2d r, vect_2d p, double mu);
double H_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2, parameters prm);
double G_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2);
//from integration.cpp
double int2d(vector <vector <double>> f, vector <double> x, vector <double> y);
double int1d_simps(vector <double> f, vector <double> x);
