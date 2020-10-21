#ifndef MSWV1_H
#define MSWV1_H
#include"SWFS.hpp"

struct MSWv1 :public SWFS
{
    std::vector<double> rho_f, mome_f, e_f, u_f, p_f, c_f, wt_f;
    std::vector<M3> Ap_f, An_f;
    std::vector<A3> Fmsw;
    MSWv1(int n);
    void update();
    void compfl();
    void Au(int n);
    void Aupc(int n);
    void Aumc(int n);
    void solve();
    // void execute();

};

#endif