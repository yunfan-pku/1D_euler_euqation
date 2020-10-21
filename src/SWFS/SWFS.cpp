#include"SWFS.hpp"

void SWFS::fu(int n)
{
    double temp = u[n] * (rho[n] - p[n] / (c[n] * c[n]));
    if (u[n] > 0)
    {
        Fp_c[n][0] = temp;
        Fp_c[n][1] = temp * u[n];
        Fp_c[n][2] = temp * u[n] * u[n] * 0.5;
    }
    else
    {
        Fn_c[n][0] = temp;
        Fn_c[n][1] = temp * u[n];
        Fn_c[n][2] = temp * u[n] * u[n] * 0.5;
    }

}
void SWFS::fupc(int n)
{
    double temp = (u[n] + c[n]) * 0.5 * p[n] / (c[n] * c[n]);
    if (u[n] + c[n] > 0)
    {
        Fp_c[n][0] = temp;
        Fp_c[n][1] = temp * (u[n] + c[n]);
        Fp_c[n][2] = temp * (u[n] * u[n] * 0.5 + u[n] * c[n] + c[n] * c[n] / (gamma - 1));
    }
    else
    {
        Fn_c[n][0] = temp;
        Fn_c[n][1] = temp * (u[n] + c[n]);
        Fn_c[n][2] = temp * (u[n] * u[n] * 0.5 + u[n] * c[n] + c[n] * c[n] / (gamma - 1));
    }
}
void SWFS::fumc(int n)
{
    double temp = (u[n] - c[n]) * 0.5 * p[n] / (c[n] * c[n]);
    if (u[n] - c[n] > 0)
    {
        Fp_c[n][0] = temp;
        Fp_c[n][1] = temp * (u[n] - c[n]);
        Fp_c[n][2] = temp * (u[n] * u[n] * 0.5 - u[n] * c[n] + c[n] * c[n] / (gamma - 1));
    }
    else
    {
        Fn_c[n][0] = temp;
        Fn_c[n][1] = temp * (u[n] - c[n]);
        Fn_c[n][2] = temp * (u[n] * u[n] * 0.5 - u[n] * c[n] + c[n] * c[n] / (gamma - 1));
    }

}
void SWFS::update()
{
    Mesh::update();
    for (int i = 0;i < nx_c;i++)
    {
        if (i == 98)
            Fp_c[i][0] = 1;
        for (int j = 0;j < 3;j++)
        {
            Fp_c[i][j] = 0.0;
            Fn_c[i][j] = 0.0;
        }
        fu(i);
        fupc(i);
        fumc(i);
    }
}
double SWFS::dt()
{
    double dtmin = 100.0;
    double temp;
    for (int i = 0;i < nx_c;i++)
    {
        temp = cfl * (x_g[i + 1] - x_g[i]) / (fabs(u[i]) + c[i]);
        if (temp < dtmin)
            dtmin = temp;
    }
    return dtmin;
}
void SWFS::compfl()
{
    for (int i = 1;i < nx_f - 1;i++)
    {
        for (int j = 0;j < 3;j++)
        {
            F[i][j] = Fp_c[i - 1][j] + Fn_c[i][j];
        }
    }
}
void SWFS::execute()
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
void SWFS::solve()
{
    update();
    compfl();
    execute();
}