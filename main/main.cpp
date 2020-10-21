#include<iostream>
#include<fstream>
#include"MSWv1.hpp"
#include"MSWv2.hpp"
using namespace std;
int main()
{
    const double gamma = 1.4;
    const int N = 81;
    MSWv2 a(N);
    //SWFS a(N);
    a.cfl = 0.9;
    for (int i = 0;i <= N;i++)
    {
        a.x_g[i] = 2.0 / N * i;
    }
    a.update_grid();
    for (int i = 0;i < N;i++)
    {
        if (a.x_c[i] <= 1.0)
        {
            a.rho[i] = 1;
            a.e[i] = 1.0 / (gamma - 1);
        }
        else
        {
            a.rho[i] = 2;
            a.e[i] = 2.0 / (gamma - 1);
        }
        a.mome[i] = 0;
    }
    for (int i = 0;i < 40;i++)
    {
        a.update();
        a.compfl();
        a.execute();
    }
    ofstream fout("out2.csv");
    for (int i = 0;i < a.nx_c;i++)
    {
        fout << a.x_c[i] << "," << a.rho[i] << "," << a.p[i] << endl;
    }

    return 0;
}
