#ifndef SWFS_H
#define SWFS_H
#include"mesh.hpp"

struct SWFS:public Mesh
{
    SWFS(int n):Mesh(n),Fp_c(nx_c),Fn_c(nx_c),cfl(0.9){}
    std::vector<A3> Fp_c,Fn_c;
    double cfl;
    double dt();
    void update();
    void fu(int n);
    void fupc(int n);
    void fumc(int n);
    void execute();
    void compfl();
    void solve();
};

#endif