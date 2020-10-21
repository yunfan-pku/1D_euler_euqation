#include"MSWv1.hpp"

MSWv1::MSWv1(int n) :SWFS(n), rho_f(nx_f), mome_f(nx_f), e_f(nx_f), u_f(nx_f), p_f(nx_f), c_f(nx_f), Ap_f(nx_f), An_f(nx_f), wt_f(nx_f), Fmsw(nx_f) {}
void MSWv1::Au(int n) {
    double ic2 = 1 / (c_f[n] * c_f[n]);
    double u2 = u_f[n] * u_f[n];
    double u3 = u2 * u_f[n];
    double A00, A01, A02, A11, A12;
    if (u_f[n] > 0)
    {
        A00 = u_f[n] - u3 * (gamma - 1) * ic2 * 0.5;
        A01 = u2 * (gamma - 1) * ic2;
        A02 = (1 - gamma) * u_f[n] * ic2;
        Ap_f[n][0][0] += A00;
        Ap_f[n][0][1] += A01;
        Ap_f[n][0][2] += A02;
        Ap_f[n][1][0] += u_f[n] * A00;
        Ap_f[n][1][1] += u_f[n] * A01;
        Ap_f[n][1][2] += u_f[n] * A02;
        Ap_f[n][2][0] += 0.5 * u3 - u3 * u2 * (gamma - 1) * 0.25 * ic2;
        Ap_f[n][2][1] += u2 * A01 * 0.5;
        Ap_f[n][2][2] += u2 * A02 * 0.5;
    }
    else
    {
        A00 = u_f[n] - u3 * (gamma - 1) * ic2 * 0.5;
        A01 = u2 * (gamma - 1) * ic2;
        A02 = (1 - gamma) * u_f[n] * ic2;
        An_f[n][0][0] += A00;
        An_f[n][0][1] += A01;
        An_f[n][0][2] += A02;
        An_f[n][1][0] += u_f[n] * A00;
        An_f[n][1][1] += u_f[n] * A01;
        An_f[n][1][2] += u_f[n] * A02;
        An_f[n][2][0] += 0.5 * u3 - u3 * u2 * (gamma - 1) * 0.25 * ic2;
        An_f[n][2][1] += u2 * A01 * 0.5;
        An_f[n][2][2] += u2 * A02 * 0.5;
    }
}


void MSWv1::Aupc(int n) {
    double ic2 = 1 / (c_f[n] * c_f[n]);
    double u2 = u_f[n] * u_f[n];
    double u3 = u2 * u_f[n];
    double cpu = u_f[n] + c_f[n];
    double A00, A01, A02, A11, A12;
    double beta = gamma - 1;
    if (u_f[n] + c_f[n] > 0)
    {
        A00 = -cpu * u_f[n] * (2 * c_f[n] + u_f[n] * (1 - gamma)) * ic2 * 0.25;
        A01 = cpu * (c_f[n] + u_f[n] * (1 - gamma)) * ic2 * 0.5;
        A02 = (gamma - 1) * cpu * ic2 * 0.5;
        Ap_f[n][0][0] += A00;
        Ap_f[n][0][1] += A01;
        Ap_f[n][0][2] += A02;
        Ap_f[n][1][0] += cpu * A00;
        Ap_f[n][1][1] += cpu * A01;
        Ap_f[n][1][2] += cpu * A02;
        Ap_f[n][2][0] += cpu * u_f[n] * (u_f[n] * beta - 2 * c_f[n]) * (2 * c_f[n] * c_f[n] + 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.125 * ic2 / beta;
        Ap_f[n][2][1] += cpu * (c_f[n] - u_f[n] * beta) * (2 * c_f[n] * c_f[n] + 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2 / beta;
        Ap_f[n][2][2] += cpu * (2 * c_f[n] * c_f[n] + 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2;
    }
    else
    {
        A00 = -cpu * u_f[n] * (2 * c_f[n] + u_f[n] * (1 - gamma)) * ic2 * 0.25;
        A01 = cpu * (c_f[n] + u_f[n] * (1 - gamma)) * ic2 * 0.5;
        A02 = (gamma - 1) * cpu * ic2 * 0.5;
        An_f[n][0][0] += A00;
        An_f[n][0][1] += A01;
        An_f[n][0][2] += A02;
        An_f[n][1][0] += cpu * A00;
        An_f[n][1][1] += cpu * A01;
        An_f[n][1][2] += cpu * A02;
        An_f[n][2][0] += cpu * u_f[n] * (u_f[n] * beta - 2 * c_f[n]) * (2 * c_f[n] * c_f[n] + 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.125 * ic2 / beta;
        An_f[n][2][1] += cpu * (c_f[n] - u_f[n] * beta) * (2 * c_f[n] * c_f[n] + 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2 / beta;
        An_f[n][2][2] += cpu * (2 * c_f[n] * c_f[n] + 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2;
    }
}

void MSWv1::Aumc(int n) {
    double ic2 = 1 / (c_f[n] * c_f[n]);
    double u2 = u_f[n] * u_f[n];
    double u3 = u2 * u_f[n];
    double umc = u_f[n] - c_f[n];
    double A00, A01, A02, A11, A12;
    double beta = gamma - 1;
    if (u_f[n] - c_f[n] > 0)
    {
        A00 = umc * u_f[n] * (2 * c_f[n] + u_f[n] * beta) * ic2 * 0.25;
        A01 = -umc * (c_f[n] + u_f[n] * beta) * ic2 * 0.5;
        A02 = beta * umc * ic2 * 0.5;
        Ap_f[n][0][0] += A00;
        Ap_f[n][0][1] += A01;
        Ap_f[n][0][2] += A02;
        Ap_f[n][1][0] += umc * A00;
        Ap_f[n][1][1] += umc * A01;
        Ap_f[n][1][2] += umc * A02;
        Ap_f[n][2][0] += umc * u_f[n] * (u_f[n] * beta + 2 * c_f[n]) * (2 * c_f[n] * c_f[n] - 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.125 * ic2 / beta;
        Ap_f[n][2][1] += -umc * (c_f[n] + u_f[n] * beta) * (2 * c_f[n] * c_f[n] - 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2 / beta;
        Ap_f[n][2][2] += umc * (2 * c_f[n] * c_f[n] - 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2;
    }
    else
    {
        A00 = umc * u_f[n] * (2 * c_f[n] + u_f[n] * beta) * ic2 * 0.25;
        A01 = -umc * (c_f[n] + u_f[n] * beta) * ic2 * 0.5;
        A02 = beta * umc * ic2 * 0.5;
        An_f[n][0][0] += A00;
        An_f[n][0][1] += A01;
        An_f[n][0][2] += A02;
        An_f[n][1][0] += umc * A00;
        An_f[n][1][1] += umc * A01;
        An_f[n][1][2] += umc * A02;
        An_f[n][2][0] += umc * u_f[n] * (u_f[n] * beta + 2 * c_f[n]) * (2 * c_f[n] * c_f[n] - 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.125 * ic2 / beta;
        An_f[n][2][1] += -umc * (c_f[n] + u_f[n] * beta) * (2 * c_f[n] * c_f[n] - 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2 / beta;
        An_f[n][2][2] += umc * (2 * c_f[n] * c_f[n] - 2 * c_f[n] * u_f[n] * beta + u2 * beta) * 0.25 * ic2;
    }
}


void MSWv1::update()
{
    SWFS::update();
    SWFS::compfl();
    for (int i = 1;i < nx_f - 1;i++)
        for (int j = 0;j < 3;j++)
            for (int k = 0;k < 3;k++)
            {
                An_f[i][j][k] = 0;
                Ap_f[i][j][k] = 0;
            }
    for (int i = 1;i < nx_f - 1;i++)
    {
        if (i == 98)
        {
            Fmsw[i][0] = 0;
        }
        rho_f[i] = (rho[i - 1] + rho[i]) * 0.5;
        mome_f[i] = (mome[i - 1] + mome[i]) * 0.5;
        e_f[i] = (e[i - 1] + e[i]) * 0.5;
        u_f[i] = mome_f[i] / rho_f[i];
        p_f[i] = (gamma - 1) * (e_f[i] - 0.5 * u_f[i] * mome_f[i]);
        c_f[i] = sqrt(gamma * p_f[i] / rho_f[i]);
        Au(i);
        Aupc(i);
        Aumc(i);
    }
}

void MSWv1::compfl()
{
    double pg;
    for (int i = 1;i < nx_f - 1;i++)
    {
        pg = (p[i + 1] - p[i]) / (p[i + 1] > p[i] ? p[i] : p[i + 1]);
        if (i == 98)
        {
            Fmsw[i][0] = 0;
        }
        wt_f[i] = 1 / (1 + pg * pg);
        wt_f[i] = 1;
        for (int j = 0;j < 3;j++)
            Fmsw[i][j] = Ap_f[i][j][0] * rho[i - 1] + Ap_f[i][j][1] * mome[i - 1] + Ap_f[i][j][2] * e[i - 1] + An_f[i][j][0] * rho[i] + An_f[i][j][1] * mome[i] + An_f[i][j][2] * e[i];
        for (int j = 0;j < 3;j++)
            F[i][j] = wt_f[i] * Fmsw[i][j] + (1 - wt_f[i]) * F[i][j];


    }
}

void MSWv1::solve()
{
    update();
    compfl();
    execute();
}
/*
void MSWv1::execute()
{
    double Dt = dt();
    double Dx;
    for (int i = 1;i < nx_c - 1;i++)
    {
        Dx = x_g[i + 1] - x_g[i];
        rho[i] += -Dt / Dx * (F[i + 1][0] - F[i][0]);
        mome[i] += -Dt / Dx * (F[i + 1][1] - F[i][1]);
        e[i] += -Dt / Dx * (F[i + 1][2] - F[i][2]);
    }
}
*/