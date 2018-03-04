#include <iostream>
#include <cstdio>
#include <cmath>
#include "algorithm"

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

double get_w(double x, double t) {
    return consts.K * pow(x, consts.alpha - 1) * exp(-consts.E / (consts.R * t));
}

double *als, *bls, *cls;
double *csucc, *freesucc;

void solve_diagonal(int xl, double **m, double *d, double *tx) {
    als[0] = -1;
    for (int i = 1; i < xl; i++) {
        als[i] = m[i][0];
    }

    cls[0] = m[0][0];
    for (int i = 1; i < xl; i++) {
        cls[i] = m[i][1];
    }

    bls[0] = m[0][1];
    bls[xl - 1] = -1;
    for (int i = 1; i < xl - 1; i++) {
        bls[i] = m[i][2];
    }

    csucc[0] = cls[0];
    for (int i = 1; i < xl; i++) {
        csucc[i] = cls[i] - als[i] / csucc[i - 1] * bls[i - 1];
    }

    freesucc[0] = d[0];
    for (int i = 1; i < xl; i++) {
        freesucc[i] = d[i] - als[i] / csucc[i - 1] * freesucc[i - 1];
    }

    tx[xl - 1] = freesucc[xl - 1] / csucc[xl - 1];
    for (int i = xl - 2; i >= 0; i--) {
        tx[i] = (freesucc[i] - bls[i] * tx[i + 1]) / csucc[i];
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

        solve_diagonal(xl, m, b, tx);
        for (int i = 0; i < xl; i++) {
            T[j + 1][i] = tx[i];
        }

        // count X
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

        solve_diagonal(xl, m, b, tx);
        for (int i = 0; i < xl; i++) {
            X[j + 1][i] = tx[i];
        }
    }
}

pair<double **, double **> start(double x_min, double x_max, double dx, double t_min, double t_max, double dt,
                                 double X0, double Xn, double T0, double Tm, double D, double C, double Q, double ro,
                                 double lambda, double alpha, double K, double E, double R) {
    conditions.x_min = x_min;
    conditions.x_max = x_max;
    conditions.dx = dx;

    conditions.t_min = t_min;
    conditions.t_max = t_max;
    conditions.dt = dt;

    conditions.X0 = X0;
    conditions.Xn = Xn;

    conditions.T0 = T0;
    conditions.Tm = Tm;

    consts.D = D;
    consts.C = C;
    consts.Q = Q;

    consts.ro = ro;
    consts.lambda = lambda;
    consts.alpha = alpha;

    consts.K = K;
    consts.E = E;
    consts.R = R;

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

    return make_pair(X, T);
}

int main(int argc, char *argv[]) {
    if (argc != ARG_COUNT) {
        cout << "argc != " << ARG_COUNT << endl;
        exit(0);
    }

    auto p = start(stod(argv[1]), stod(argv[2]), stod(argv[3]), stod(argv[4]), stod(argv[5]), stod(argv[6]),
                   stod(argv[7]), stod(argv[8]), stod(argv[9]), stod(argv[10]), stod(argv[11]), stod(argv[12]),
                   stod(argv[13]), stod(argv[14]), stod(argv[15]), stod(argv[16]), stod(argv[17]), stod(argv[18]),
                   stod(argv[19]));

    int xl = int((conditions.x_max - conditions.x_min) / conditions.dx);
    int tl = int((conditions.t_max - conditions.t_min) / conditions.dt);

    double **X = p.first;
    double **T = p.second;

    freopen("X.out", "w+", stdout);

    for (int i = 0; i < tl; i++) {
        for (int j = 0; j < xl; j++) {
            printf("%.12f ", X[i][j]);
        }
        printf("\n");
    }

    freopen("T.out", "w+", stdout);

    for (int i = 0; i < tl; i++) {
        for (int j = 0; j < xl; j++) {
            printf("%.3f ", T[i][j]);
        }
        printf("\n");
    }
}