#include <iostream>
#include <cstdio>
#include <cmath>

using namespace std;

int ARG_COUNT = 20;

struct Consts {
    double D, C, Q, K, E, R;
    double ro, lambda, alpha;
} consts;

struct Conditions {
    double x_min, t_min, x_max, t_max, dt, dx;
    double T0, Tm, X0, Xn;
} conditions;

void zeroing_out(int d1, int d2, double **array) {
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            array[i][j] = 0;
        }
    }
}

void zeroing_out(int d1, double *array) {
    for (int i = 0; i < d1; i++) {
        array[i] = 0;
    }
}

double get_w(double x, double t) {
    return consts.K * pow(x, consts.alpha - 1) * exp(-consts.E / (consts.R * t));
}

double *als, *bls, *cls;
double *csucc, *freesucc;

void solve_diagonal_java(int xl, double **m, double *d, double *tx) {
    zeroing_out(xl, als);
//    als[0] = -1;
    for (int i = 1; i < xl; i++) {
        als[i] = m[i][0];
    }

    zeroing_out(xl, cls);
    cls[0] = m[0][0];
    for (int i = 1; i < xl; i++) {
        cls[i] = m[i][1];
    }

    zeroing_out(xl, bls);
    bls[0] = m[0][1];
//    bls[xl - 1] = -1;
    for (int i = 1; i < xl - 1; i++) {
        bls[i] = m[i][2];
    }

    // suck dick
    zeroing_out(xl, csucc);
    csucc[0] = cls[0];
    for (int i = 1; i < xl; i++) {
        csucc[i] = cls[i] - als[i] / csucc[i - 1] * bls[i - 1];
    }

    zeroing_out(xl, freesucc);
    freesucc[0] = d[0];
    for (int i = 1; i < xl; i++) {
        freesucc[i] = d[i] - als[i] / csucc[i - 1] * freesucc[i - 1];
    }

    tx[xl - 1] = freesucc[xl - 1] / csucc[xl - 1];
    for (int i = xl - 2; i >= 0; i--) {
        tx[i] = (freesucc[i] - bls[i] * tx[i + 1]) / csucc[i];
    }
}

void solve_diagonal(int xl, int tl, double **m, double *d, double *tx) {
    int n = xl - 1;
    zeroing_out(xl, als);
    zeroing_out(xl, bls);

    als[0] = -m[0][1] / m[0][0];
    bls[0] = -d[0] / m[0][0];

    for (int i = 1; i < n; i++) {
        double a = m[i][i - 1];
        double b = m[i][i];
        double c = m[i][i + 1];
        als[i] = -c / (a * als[i - 1] + b);
        bls[i] = (d[i] - a * bls[i - 1]) / (a * als[i - 1] + b);
    }

    tx[n] = (d[n] - m[n][n - 1] * m[n - 1][n - 1]) / (m[n][n - 1] * als[n - 1] + m[n][n]);

    for (int i = n - 1; i > -1; i--) {
        tx[i] = als[i] * tx[i + 1] + bls[i];
    }
}

void solve(int xl, int tl, double **X, double **T) {
    double **m = new double *[xl];
    for (int i = 0; i < xl; i++) {
        m[i] = new double[3];
    }
    double *b = new double[xl];
    double *tx = new double[xl];

    als = new double[xl];
    bls = new double[xl];
    cls = new double[xl];

    csucc = new double[xl];
    freesucc = new double[xl];

    for (int j = 0; j < tl - 1; j++) {
        // count T
        zeroing_out(xl, xl, m);
        zeroing_out(xl, b);

        for (int i = 0; i < xl; i++) {
            double w = get_w(X[j][i], T[j][i]) * X[j][i];
            double k1 = consts.lambda / (conditions.dx * conditions.dx);
            double k2 = consts.ro * consts.C / conditions.dt;
            b[i] = (k2 * T[j][i]) + (consts.ro * consts.Q * w);

            m[i][0] = -k1;
            m[i][1] = k2 + (2 * k1);
            m[i][2] = -k1;
        }

        m[0][0] = 1;
        m[0][1] = 0;
        m[xl - 1][0] = -1;
        m[xl - 1][1] = 1;

        b[0] = conditions.Tm;
        b[xl - 1] = 0;

        zeroing_out(xl, tx);
        solve_diagonal_java(xl, m, b, tx);
        for (int i = 0; i < xl; i++) {
            T[j + 1][i] = tx[i];
        }

        // count X
        zeroing_out(xl, 3, m);
        zeroing_out(xl, b);
        for (int i = 0; i < xl; i++) {
            b[i] = X[j][i] / conditions.dt;
            double kappa = consts.D / (conditions.dx * conditions.dx);
            double w = get_w(X[j][i], T[j][i]);

            m[i][0] = -kappa;
            m[i][1] = (2 * kappa) + (1 / conditions.dt) + w;
            m[i][2] = -kappa;
        }
        m[0][0] = 1;
        m[0][1] = 0;
        m[xl - 1][0] = -1;
        m[xl - 1][1] = 1;

        b[0] = 1;
        b[xl - 1] = 0;

        zeroing_out(xl, tx);
        solve_diagonal_java(xl, m, b, tx);
        for (int i = 0; i < xl; i++) {
            X[j + 1][i] = tx[i];
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != ARG_COUNT) {
        cout << "argc != " << ARG_COUNT << endl;
        exit(0);
    }

    conditions.x_min = stod(argv[1]);
    conditions.x_max = stod(argv[2]);
    conditions.dx = stod(argv[3]);

    conditions.t_min = stod(argv[4]);
    conditions.t_max = stod(argv[5]);
    conditions.dt = stod(argv[6]);

    conditions.X0 = stod(argv[7]);
    conditions.Xn = stod(argv[8]);

    conditions.T0 = stod(argv[9]);
    conditions.Tm = stod(argv[10]);

    consts.D = stod(argv[11]);
    consts.C = stod(argv[12]);
    consts.Q = stod(argv[13]);

    consts.ro = stod(argv[14]);
    consts.lambda = stod(argv[15]);
    consts.alpha = stod(argv[16]);

    consts.K = stod(argv[17]);
    consts.E = stod(argv[18]);
    consts.R = stod(argv[19]);

    int xl = int((conditions.x_max - conditions.x_min) / conditions.dx);
    int tl = int((conditions.t_max - conditions.t_min) / conditions.dt);

    cout << tl << " " << xl << endl;

    double **X, **T;
    X = new double *[tl];
    T = new double *[tl];
    for (int i = 0; i < tl; i++) {
        X[i] = new double[xl];
        T[i] = new double[xl];
    }

    for (int j = 0; j < xl; j++) {
        for (int i = 0; i < tl; i++) {
            X[i][j] = 0;
            T[i][j] = 0;
        }
        X[0][j] = conditions.X0;
        T[0][j] = conditions.T0;
    }
    X[0][0] = conditions.Xn;
    T[0][0] = conditions.Tm;

    solve(xl, tl, X, T);

    freopen("X.out", "w+", stdout);

    for (int i = 0; i < tl; i++) {
        for (int j = 0; j < xl; j++) {
            printf("%.12f ", X[i][j]);
        }
        cout << endl;
    }

    freopen("T.out", "w+", stdout);

    for (int i = 0; i < tl; i++) {
        for (int j = 0; j < xl; j++) {
            cout << T[i][j] << " ";
        }
        cout << endl;
    }
}