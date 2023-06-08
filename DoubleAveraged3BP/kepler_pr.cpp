#include "Header.h"

double scal_2d(vect_2d r, vect_2d d) {
	return r.x * d.x + r.y * d.y;
}

double mod_2d(vect_2d r) {
	return r.x * r.x + r.y * r.y;
}

// start values from kepler to (p,r)

vect_2d r1_0(double a, double e, double omega_d) {
	vect_2d r;
	r.x = a * (1 - e) * cos(omega_d);
	r.y = a * (1 - e) * sin(omega_d);
	return r;
}

vect_2d r2_0(double a, double e) {
	vect_2d r;
	r.x = a * (1 - e);
	r.y = 0;
	return r;
}

vect_2d p1_0(double a, double mu, double e, double omega_d) {
	vect_2d r;
	r.x = -mu * sin(omega_d) * sqrt((1 - e * e) / a);
	r.y = mu * cos(omega_d) * sqrt((1 - e * e) / a);
	return r;
}

vect_2d p2_0(double a, double mu, double e, double omega_d) {
	vect_2d r;
	r.x = 0;
	r.y = sqrt((1 - e * e) / a);
	return r;
}

// start values from (p,r) to kepler

double a_from_pr(vect_2d r, vect_2d p, double mu) {
	return 1 / (mod_2d(p) / (mu * mu) - 2 / sqrt(mod_2d(r)));
}

double e_from_pr(vect_2d r, vect_2d p, double mu) {
	double a = a_from_pr(r, p, mu);
	return sqrt(1 - 1 / (a * mu * mu) * (mod_2d(p) * mod_2d(r) - scal_2d(r, p) * scal_2d(r, p)));
}

double true_anomaly_from_pr(vect_2d r, vect_2d p, double mu) {
	double a = a_from_pr(r, p, mu);
	double e = e_from_pr(r, p, mu);
	double cos = 1 / e * (a * (1 - e * e) / sqrt(mod_2d(r)) - 1);
	double sin = scal_2d(r, p) / (sqrt(mod_2d(r)) * mu * e) * sqrt(a * (1 - e * e));
	if (abs(sin) < 0.000001) {
		return cos / abs(cos) * (asin(sin) - pi/2) + pi/2;
	}
	return acos(cos) * sin / abs(sin);
}

double omega_from_pr(vect_2d r, vect_2d p, double mu) {
	double tr_a = true_anomaly_from_pr(r, p, mu);
	double sin = asin(r.y / sqrt(mod_2d(r)));
	double cos = acos(r.x / sqrt(mod_2d(r)));
	if (abs(sin) < 0.000001) {
		return cos / abs(cos) * (asin(sin) - pi / 2) + pi / 2 - tr_a;
	}
	return acos(cos) * sin / abs(sin) - tr_a;
}

double H_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2, parameters prm) {
	double H01 = mod_2d(p1) / (2 * prm.mu_1) - prm.mu_1 / sqrt(mod_2d(r1));
	double H02 = mod_2d(p2) / (2 * (1 - prm.mu_1)) - (1 - prm.mu_1) / sqrt(mod_2d(r2));
	vect_2d dr = { r1.x - r2.x, r1.y - r2.y };
	vect_2d sum_p = { p1.x + p2.x, p1.y + p2.y };
	double V = prm.mu_1 * (1 - prm.mu_1) / (sqrt(mod_2d(dr)));
	V += - prm.mu_1 / sqrt(mod_2d(r1)) - (1 - prm.mu_1) / sqrt(mod_2d(r2));
	V += -1 / (2 * (1 - prm.mu)) * mod_2d(sum_p);
	return H01 + H02 - prm.mu * V;
}

double G_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2) {
	return (r1.x * p1.y - r1.y * p1.x) + (r2.x * p2.y - r2.y * p2.x);
}