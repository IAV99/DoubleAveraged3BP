#include "Header.h"
#include "Integrator.h"

void evolution_averaged(const parameters& p, const dWd& dW, const double t, const Array& x, Array* dx) {
    double mu_2 = 1 - p.mu_1;
    vector <vector <double>> dWde10;
    for (int i = 0; i < dW.n; i++) {
        vector <double> w;
        for (int j = 0; j < dW.m; j++) {
            double w_0 = dW_de1(p.a, x[2], x[0], x[1], dW.E1[i], dW.E2[j]);
            w.push_back(w_0);
        }
        dWde10.push_back(w);
    }
    double dWde1 = p.mu_1 * mu_2 / (4 * pi * pi) * int2d(dWde10, dW.E1, dW.E2);
    vector <vector <double>> dWde20;
    for (int i = 0; i < dW.n; i++) {
        vector <double> w;
        for (int j = 0; j < dW.m; j++) {
            double w_0 = dW_de2(p.a, x[2], x[0], x[1], dW.E1[i], dW.E2[j]);
            w.push_back(w_0);
        }
        dWde20.push_back(w);
    }
    double dWde2 = p.mu_1 * mu_2 / (4 * pi * pi) * int2d(dWde20, dW.E1, dW.E2);
    vector <vector <double>> dWdo;
    for (int i = 0; i < dW.n; i++) {
        vector <double> w;
        for (int j = 0; j < dW.m; j++) {
            double w_0 = dW_domegad(p.a, x[2], x[0], x[1], dW.E1[i], dW.E2[j]);
            w.push_back(w_0);
        }
        dWdo.push_back(w);
    }
    double dWdOmega = p.mu_1 * mu_2 / (4 * pi * pi) * int2d(dWdo, dW.E1, dW.E2);

    (*dx)[0] = - sqrt(1 - x[0] * x[0]) * dWdOmega/ (sqrt(p.a) * x[0] * p.mu_1);
    (*dx)[1] = sqrt(1 - x[1] * x[1]) * dWdOmega / (x[1] * mu_2);
    (*dx)[2] = (sqrt(1 - x[0] * x[0]) * dWde1 / (sqrt(p.a) * x[0] * p.mu_1) - sqrt(1 - x[1] * x[1]) * dWde2 / (x[1] * mu_2));
}

int main() {
    double mu = 1e-4;
    vector <double> E1;
    int n = 101;
    char answ;
    // arrays of eccentric anomaly
    cout << "rewrite E1? [y/n]" << endl;
    cin >> answ;
    if (answ == 'y') {
        ofstream out;
        out.open("E1_E2.txt");
        cout << "input n" << endl;
        cin >> n;
        out << n << endl;;
        for (int i = 0; i < n; i++) {
            E1.push_back((double)i * 2 * pi / ((double)n - 1));
            out << E1[i] << " ";
        }
        out << endl;
        out.close();
    }
    else {
        ifstream in("E1_E2.txt");
        if (in.is_open()) {
            in >> n;
            for (int i = 0; i < n; i++) {
                double a;
                in >> a;
                E1.push_back(a);
            }
        }
        in.close();
    }
    const dWd dW = {E1, E1, n, n};
    // parameters
    double e_1 = 0.7, omega_d = 0, a = 0.3, mu_1 = 0.5, sigma = 0.8, e_2 = e2(a, mu_1, sigma, e_1);
    cout << "use new parameters? [y/n]" << endl;

    cin >> answ;
    if (answ == 'y') {
        cout << "input mu, a, sigma, mu_1" << endl;
        cin >> mu >> a >> sigma >> mu_1;
        cout << "input omega_d" << endl;
        cin >> omega_d;
        cout << "[ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
        cout << "input e1" << endl;
        cin >> e_1;
        while (e_1 > e1_max(a, mu_1, sigma) || e_1 < e1_min(a, mu_1, sigma)) {
            cout << "e1 not in [ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
            cout << "input e1" << endl;
            cin >> e_1;
        }
        e_2 = e2(a, mu_1, sigma, e_1);

        ofstream out;
        out.open("parameters.txt");
        out << mu << " " << a << " " << sigma << " " << mu_1 << endl;
        out << omega_d << " " << e_1 << endl;
        out.close();

        while ((a * (1 - e_1 * e_1) + 1 / a * (1 - e_2 * e_2)) - 2 * (1 - e_1 * e_2 * cos(omega_d)) < 0) {
            cout << "parameters: a = " << a << "; mu_1 = " << mu_1 << "; sigma = " << sigma << endl;
            cout << "orbits intersect!" << endl;
            cout << "input omega_d" << endl;
            cin >> omega_d;
            cout << "[ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
            cout << "input e1" << endl;
            cin >> e_1;
            while (e_1 > e1_max(a, mu_1, sigma) || e_1 < e1_min(a, mu_1, sigma)) {
                cout << "e1 not in [ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
                cout << "input e1" << endl;
                cin >> e_1;
            }
            e_2 = e2(a, mu_1, sigma, e_1);
        }
    }
    else {
        ifstream in("parameters.txt");
        if (in.is_open()) {
            in >> mu >> a >> sigma >> mu_1;
            in >> omega_d >> e_1;
        }
        in.close();
    }
    double mu_2 = 1 - mu_1;
    cout << "parameters: a = " << a << "; mu_1 = " << mu_1 << "; sigma = " << sigma << endl;
    cout << "start from: e1 = " << e_1 << "; e2 = " << e_2 << "; omega_d = " << omega_d << endl;
    const parameters p = {mu, a, mu_1, sigma};

    Integrator integrator;
    Array x;
    //initial condition vector
    x[0] = e_1;
    x[1] = e_2;
    x[2] = omega_d;

    //integrator parameters
    double t = 0;
    double h = 0.001;

    double t_end = 10000;
    double t_save = 1;
    int save_counter = 1, counter = 0;

    ofstream out;
    out.open("results_for_e1_e2_omega.txt");
    out << setprecision(12) << p.mu << ", " << p.a << ", " << p.mu_1 << ", " << p.sigma << endl;

    vector <vector <double>> W;
    for (int i = 0; i < n; i++) {
        vector <double> w;
        for (int j = 0; j < n; j++) {
            double w_0 = W0(p.a, x[2], x[0], x[1], dW.E1[i], dW.E1[j]);
            w.push_back(w_0);
        }
        W.push_back(w);
    }
    double W_avg0 = mu_1 * mu_2 / (4 * pi * pi) * int2d(W, E1, E1);

    cout << "t,   e1,   e2,   omega,   sigma - sigma0,   W_t - W0" << endl;

    while (t < t_end) {
        integrator.rk78(&t, &x, &h, 1e-9, 1e-6, 1e-1, [&p, &dW](const double t, const Array& x, Array* dx) { evolution_averaged(p, dW, t, x, dx); });

        vector <vector <double>> W;
        for (int i = 0; i < n; i++) {
            vector <double> w;
            for (int j = 0; j < n; j++) {
                double w_0 = W0(p.a, x[2], x[0], x[1], dW.E1[i], dW.E1[j]);
                w.push_back(w_0);
            }
            W.push_back(w);
        }
        double W_avg = mu_1 * mu_2 / (4 * pi * pi) * int2d(W, E1, E1);

        if (counter % save_counter == 0) {
            out << setprecision(12) << t / p.mu << ", " << x[0] << ", " << x[1] << ", " << x[2] << ", ";
            out << Sigma(p.a, p.mu_1, x[0], x[1]) << ", " << W_avg << endl;

            cout << setprecision(3) << t / p.mu << ", " << x[0] << ", " << x[1] << ", " << x[2] << ", ";
            cout << Sigma(p.a, p.mu_1, x[0], x[1]) - sigma << ", " << W_avg - W_avg0 << endl;
        }
        counter++;
    }
    out.close();
    return 0;
}
