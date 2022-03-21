#include<iostream>
#include<fstream>
#include <cmath>
#include <vector>
using namespace std;


void RunPhysic(vector<double> &T)
{
    double Tc, Te, Tw, Tn, Ts;
    double FluxC, FluxE, FluxW, FluxN, FluxS;
    double dx = 0;
    int Tb = 240;
    double Tb0 = 0;
    (void)Tb0;

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


void Static() {
    const int size = 5;
    const int n = size * size;
    const int n_link = size - 1;
    const double Tbl = 1, Tbr = 5;
    const double k_ = 1;
    const double dh = 1;
    double dx = 0;
    vector<double> mtx(n, 0);// ae(4, 0), aw(4, 0), ap(4, 0);
    auto ind=[size](int i){
        return i+size*i;
    };
    auto indup=[size](int i){
        return i * (size + 1) + 1;
    };
    auto inddown=[size](int i){
        return (i + 1) * size + i;
    };
    vector<double> rhs(size, 0), S(size, 0);
    ofstream f("Static.xls");
    double ap(0), ae(0), aw(0);
    //diagonal elements
    double Tb;
    for (int i = 0; i < size; i++) {
        dx = dh;
        Tb = 0;
        ap = 2 * k_ / dx;
        if (i == 0) {
            //dx *= 0.5;
            Tb = 2 * Tbl;
            ap = 3 * k_ / dx;
        }
        if (i == size - 1) {
            //dx *= 0.5;
            Tb = 2 * Tbr;
            ap = 3 * k_ / dx;
        }
        rhs[i] = Tb / dx;
        mtx[ind(i)] = ap;
    }
    //links
    for (int i = 0; i < n_link; i++) {
        auto flux = k_ / dh;
        mtx[indup(i)] = -flux;
        mtx[inddown(i)] = -flux;
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
    // Gauss method for system J * S = F
    // forward stroke
    for (int k = 0; k < size; k++) {
        for (int j = k + 1; j < size; j++) {
            double d = mtx[j * size + k] / mtx[k * size + k];
            for (int i = k; i < size; i++) {
                mtx[j * size + i] = mtx[j * size + i] - d * mtx[k * size + i];
            }
            rhs[j] = rhs[j] - d * rhs[k];
        }
    }
    // reverse stroke
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
        //f << S[i] << "\t";
    }
    for (int i = 0; i < size; i++) {
        f << scientific<<S[i] << "\t" << endl;
    }
    f << endl;
    f.close();
}

void GaussSolver(const int size, double* mtx, double* rhs, double* S){
    // Gauss method for system J * S = F
    // forward stroke
    ofstream f("rhs.xls");
    for (int k = 0; k < size; k++) {
        for (int j = k + 1; j < size; j++) {
            double d = mtx[j * size + k] / mtx[k * size + k];
            for (int i = k; i < size; i++) {
                mtx[j * size + i] = mtx[j * size + i] - d * mtx[k * size + i];
            }
            rhs[j] = rhs[j] - d * rhs[k];
            f << scientific << rhs[j] << "\t" << endl;
        }
    }
    // reverse stroke
    S[size - 1] = rhs[size - 1] / mtx[(size - 1) * size + size - 1];
    for (int k = size - 2; k >= 0; k--) {
        double d(0);
        for (int j = k + 1; j < size; j++) {
            double ss = mtx[k * size + j] * S[j];
            d += ss;
        }
        S[k] = (rhs[k] - d) / mtx[k * size + k];
    }
    f.close();
}

double norm(vector<double> &a)
{
    double sumSqr = 0;
    for(int i = 0; i < a.size(); i++){
        sumSqr += a[i]*a[i];
    }
    return sqrt(sumSqr);
}


void NotStatic() {
    const int size = 4;
    const int n = size * size;
    const int n_link = size - 1;
    //const int time = size - 2;
    const double dh = 1;
    const double pc = 1;

    //source is taking into account
    const double Sp = 0;
    const double Sc = 0;

    const double T0 = 0;
    double dt = 0.1;
    double timing = 0;
    double T = 10;
    double eps = 0.3;
    double Tbl = 0, Tbr = 5;
    double k_ = 1;
    double dx = 0;
    double b = 0;
    double Tb = 0;
    double ap(0), ap0(0), ae(0), aw(0);
    vector<double> mtx(n, 0), rhs(size, 0), S(size, 0), temp(size, 0), diff(size, 0);// ae(4, 0), aw(4, 0), ap(4, 0);
    auto ind=[size](int i){
        return i+size*i;
    };
    auto indup=[size](int i){
        return i * (size + 1) + 1;
    };
    auto inddown=[size](int i){
        return (i + 1) * size + i;
    };
    ofstream f("NotStatic.xls");
    //diagonal elements
    while (timing < T) {
        timing += dt;
        for (int i = 0; i < size; i++) {
            dx = dh;
            Tb = 0;
            ap0 = pc * dx * 0.5 / dt;
            ap = 2 * k_ / dx + ap0 - Sp * dx * 0.5;
            b = Sc * 0.5 * dx + ap0 * temp[i];
            if (i == 0) {
                //dx *= 0.5;
                Tb = 2 * Tbl;
                ap = 3 * k_ / dx + ap0 - Sp * dx * 0.5;
            }
            if (i == size - 1) {
                //dx *= 0.5;
                Tb = 2 * Tbr;
                ap = 3 * k_ / dx + ap0 - Sp * dx * 0.5;
            }
            rhs[i] = Tb / dx + b;
            mtx[ind(i)] = ap;
        }
        //links
        for (int i = 0; i < n_link; i++) {
            auto flux = k_ / dh;
            mtx[indup(i)] = -flux;
            mtx[inddown(i)] = -flux;
        }
        cout << "Time: " << timing << endl << endl;
        if (timing == dt) {
            cout << "Matrix:";
            for (int i = 0; i < n; i++) {
                if (i % size == 0) {
                    cout << endl;
                }
                cout << mtx[i] << "\t";
            }
            cout << endl << endl;
        cout << "Rhs:" << endl;
        for (int i = 0; i < size; i++) {
            cout << rhs[i] << endl;
        }
        }
        cout << endl;
        GaussSolver(size, mtx.data(), rhs.data(), S.data());
        for (int i = 0; i < size; i++) {
            diff[i] = S[i] - temp[i];
        }
        for (int i = 0; i < size; i++) {
            {
                if (norm(diff) == eps) {
                    break;
                }
                temp[i] = S[i];
            }
        }
        cout << "Temp:" << endl;
        for (int i = 0; i < size; i++) {
            cout << temp[i] << endl;
        }
    }

    cout << endl << "Result:" << endl;
    for (int i = 0; i < size; i++) {
        cout << temp[i] << endl;
        //f << S[i] << "\t";
    }
    for (int i = 0; i < size; i++) {
        f << scientific << temp[i] << "\t" << endl;
    }
    f << endl;
    f.close();
}

void TwoDim() {
    const int nx = 5;
    const int ny = 5;
    const int size = nx * ny;
    const int n = size * size;
    const int n_link = size - 1;
    //const int time = size - 2;
    const double dh = 1;
    const double pc = 1;

    //source is taking into account
    const double Sp = 0;
    const double Sc = 0;

    double dt = 0.1;
    double timing = 0;
    double T = 10;
    double eps = 0.3;
    double Tbl = 5, Tbr = 5, Tbd = 5, Tbu = 5;
    double k_ = 1;
    double dx = 0;
    double b = 0;
    double Tb1 = 0;
    double Tb2 = 0l;
    double ap(0), ap0(0), ae(0), aw(0);
    vector<double> mtx(n, 0), rhs(size, 0), S(size, 0), temp(size, 0), diff(size, 0);// ae(4, 0), aw(4, 0), ap(4, 0);
    auto ind=[size](int i){
        return i+size*i;
    };
    auto indup=[size](int i){
        return i * (size + 1) + 1;
    };
    auto inddown=[size](int i){
        return (i + 1) * size + i;
    };
    auto inds=[size](int i){
        return (i - nx + 1) * size + i + 1;
    };
    auto indn=[size](int i){
        return (i + nx - 1) * size + size + i;
    };
    while (timing < T) {
        timing += dt;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                dx = dh;
                ap0 = pc * pow(dx, 2) * 0.25 / dt;
                ap = 2 * k_ + ap0 - Sp * pow(dx, 2) * 0.25;
                b = Sc * 0.25 * pow(dx, 2) + ap0 * temp[i];
                if (i == 0) {
                    //dx *= 0.5;
                    Tb1 = 2 * Tbl;
                    ap = 3 * k_ + ap0 - Sp * pow(dx, 2) * 0.25;
                }
                if (i == size - 1) {
                    //dx *= 0.5;
                    Tb1 = 2 * Tbr;
                    ap = 3 * k_ + ap0 - Sp * pow(dx, 2) * 0.25;
                }
                if (j == 0) {
                    Tb2 = 2 * Tbd;
                    ap = 3 * k_ + ap0 - Sp * pow(dx, 2) * 0.25;
                }
                if (j == size - 1) {
                    Tb2 = 2 * Tbu;
                    ap = 3 * k_ + ap0 - Sp * pow(dx, 2) * 0.25;
                }
                if (i < nx) {
                    rhs[i] = Tb2 / dx + b;
                    if (i == 0 || i == nx - 1) {
                        rhs[i] += Tb1 / dx;
                    }
                }
                if (i == nx) {
                    rhs[i] = Tb1 / dx + b;
                }
                if (i != 0 && i % nx == 0 && i != size - 1 && i != nx) {
                    rhs[i] = Tb1 / dx + b;
                    rhs[i - 1] = Tb1 / dx + b;
                }
                if (i == (size - nx)) {
                    rhs[i] += Tb2 / dx;
                }
                if (i > (size - nx)) {
                    rhs[i] = Tb2 / dx + b;
                    if (i == size - 1) {
                        rhs[i] += Tb1 / dx;
                    }
                }
                mtx[ind(i)] = ap;
            }
        }
        //links
        for (int i = 0; i < n_link; i++) {
            auto flux = k_ * 0.5;
            mtx[indup(i)] = -flux;
            mtx[inddown(i)] = -flux;
            mtx[inds(i)] = -flux;
            mtx[indn(i)] = -flux;
        }
        cout << "Time: " << timing << endl << endl;
        if (timing == dt) {
            cout << "Matrix:";
            for (int i = 0; i < n; i++) {
                if (i % size == 0) {
                    cout << endl;
                }
                cout << mtx[i] << "\t";
            }
            cout << endl << endl;
            cout << "Rhs:" << endl;
            for (int i = 0; i < size; i++) {
                cout << rhs[i] << endl;
            }
        }
        cout << endl;
        GaussSolver(size, mtx.data(), rhs.data(), S.data());
        for (int i = 0; i < size; i++) {
            diff[i] = S[i] - temp[i];
        }
        for (int i = 0; i < size; i++) {
            {
                if (norm(diff) == eps) {
                    break;
                }
                temp[i] = S[i];
            }
        }
        cout << "Temp:" << endl;
        //for (int i = 0; i < size; i++) {
        //    cout << temp[i] << endl;
        //}
        for (int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                cout << temp[i * nx + j] << "\t";
            }
            cout << "\n";
        }
    }
    cout << endl << "Result:" << endl;
    for (int i = 0; i < size; i++) {
        cout << temp[i] << endl;
        //f << S[i] << "\t";
    }
}
/*
void tochn (double alpha, double beta, double k, double pi, double tau, double h, int N, double T, int l) {
    double homogeneous, inhomogeneous, timing=0;
    double u[N];
    int n=10;

    while (timing<T) {

        for (int i=0; i<=N; i++)
        {
            //решение однородного уравнение (𝛏't=𝛏''xx)
            homogeneous=alpha*timing*(1-h*i)+beta*timing*h*i+(exp(-pow(pi*k,2)*timing))*sin(pi*k*h*i);

            //решение неоднородного уравнение (𝛈't=𝛈''xx - 𝛍...)
            for(int j=1; j<=n; j++)
            {
                inhomogeneous=inhomogeneous+(2/(pow(pi*j,3)))*(alpha-(pow(-1, j))*beta)*sin(pi*j*h*i)*(exp(-pow(pi*j,2)*timing)-1);
            }
            u[i]=homogeneous+inhomogeneous;
            u[0]=alpha*timing;
            u[N]=beta*timing;
            ex[l][0]=u[0];
            ex[l][N]=u[N];

            cout << u[i] << " ";
        }
        timing=timing+tau;
        l=l+1;
    }
}
*/

int main() {
    //Static();
    //NotStatic();
    TwoDim();
}