#include<iostream>
#include<fstream>
#include <vector>
using namespace std;


void RunPhysic(vector<double> &T)
{
    double Tc, Te, Tw, Tn, Ts;
    double FluxC, FluxE, FluxW, FluxN, FluxS;
    double dx = 0;
    int Tb = 240;
    double Tb0 = 0;

    int i, j;
    int imax = 40;
    int jmax = 40;
    double dh = 1;
    double k = 1;
    double delt = 0.2;
    double FaceArea = 1;
    ofstream f("result1.xls");
    for (i = 0; i < imax; i++) {
        Tc = T[i];
        dx = dh;

        if (i == imax - 1) {
            Te = Tb0;
            dx *= 0.5;
        }
        else
            Te = T[i + 1];
        FluxE = (-k * FaceArea) / dx;

        if (i == 0) {
            Tw = Tb;
            dx *= 0.5;
        }
        else
            Tw = T[i - 1];
        FluxW = (-k * FaceArea) / dx;

        FluxC = FluxE + FluxW;

        T[i] = Tc + delt * (FluxC * Tc - (FluxE * Te + FluxW * Tw));
        f << T[i] << "\t";
        f << endl;
    }
    f.close();
}


void TriDiagMtx() {
    int n = 16, size = 4;
    int n_link = size - 1;
    double Tbl = 1, Tbr = 5;
    double k_ = 1;
    double dx = 0;
    double dh = 1;
    vector<double> mtx(n, 0);// ae(4, 0), aw(4, 0), ap(4, 0);
    vector<double> rhs(size, 0), S(size, 0);
    /*
    for (int i = 0; i < n; i++) {
        dx = dh;
        if (i == 0) {
            dx *= 0.5;
        }
        if (i == n - 1) {
            dx *= 0.5;
        }
        ae[i] = k / dx;
        aw[i] = k / dx;
        ap[i] = ae[i] + aw[i];
        mtx[i*size + i] = ap[i];
        mtx[i*size + i + 1] = ae[i];
        mtx[i*size + i - 1] = aw[i];
    }
     */
    double ap(0), ae(0), aw(0);
    //diagonal elements
    double Tb = 0;
    for (int i = 0; i < size; i++) {
        dx = dh;
        Tb = 0;
        if (i == 0) {
            dx *= 0.5;
            Tb = Tbl;
        }
        if (i == size - 1) {
            dx *= 0.5;
            Tb = Tbr;
        }
        ap = k_ / dx;
        rhs[i] = Tb / dx;
        mtx[i * size + i] = ap;
    }
    //links
    for (int i = 0; i < n_link; i++) {
        auto flux = k_ / dh;
        mtx[i * (size + 1) + 1] = flux;
        mtx[(i + 1) * size + i] = -flux;
    }
    cout << "Matrix:";
    for (int i = 0; i < n; i++) {
        if (i % size == 0) {
            cout << endl;
        }
        cout << mtx[i] << "\t";
    }
    cout << endl << endl << "Rhs:" << endl;
    for (int i = 0; i < size; i++) {
        cout << rhs[i] << endl;
    }
    // решение методом Гауса линейной системы уравнений J * S = F

    // прямой ход

    for (int k = 0; k < size; k++) {
        for (int j = k + 1; j < size; j++) {
            double d = mtx[j * size + k] / mtx[k * size + k];
            for (int i = k; i < size; i++) {
                mtx[j * size + i] = mtx[j * size + i] - d * mtx[k * size + i];
            }
            rhs[j] = rhs[j] - d * rhs[k];
        }
    }
    // обратный ход
    S[size - 1] = rhs[size - 1] / mtx[(size - 1) * size + size - 1];
    for (int k = size - 2; k >= 0; k--) {
        double d(0);
        for (int j = k + 1; j < size; j++) {
            double ss = mtx[k * size + j] * S[j];
            d += ss;
        }
        S[k] = (rhs[k] - d) / mtx[k * size + k];
    }
    cout << endl << "Result:" << endl;
    for (int i = 0; i < size; i++) {
        cout << S[i] << endl;
    }
}

int main() {
    int n = 40;
    double dh = 1;
    double dx = 0;
    double a = dh / 2;
    vector<double> T(n, 0), c(n, 0), b(n, 0), p(n, 0), q(n, 0);
    RunPhysic(T);
    b[n] = 0;
    for (int i = 0; i < n - 1; i++) {
        b[i] = 1;
    }
    c[0] = 0;
    for (int i = 1; i < n; ++i) {
        c[i] = 1;
    }
    p[0] = b[0] / a;
    q[0] = 0;
    ofstream f("result.xls");
    for (int i = 1; i < n; i++) {
        dx = dh;
        if (i == n - 1) {
            dx *= 0.5;
        }
        p[i] = b[i] / (a - c[i] * p[i-1]);
        q[i] = c[i] * q[i-1] / (a - c[i] * p[i-1]);
    }
    T[n] = q[n];
    for (int i = 0; i < n - 1; i++) {
        T[i] = p[i] * T[i + 1] + q[i];
    }
    for (int i = 0; i < n; i++) {
        f << T[i] << "\t";
        f << endl;
    }
    f.close();
    TriDiagMtx();
}
