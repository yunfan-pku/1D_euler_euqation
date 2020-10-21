#ifndef MESH_H
#define MESH_H
#include<vector>
#include<cmath>
struct A3
{
    double a_[3];
    double& operator[](int i){return a_[i];}
};
struct M3
{
    A3 a_[3];
    A3& operator[](int i){return a_[i];}
};

struct Mesh
{
    const double gamma = 1.4;
    Mesh(int);
    int nx_c, nx_f;
    std::vector<double> x_g, x_c;
    std::vector<double> rho, mome, e, u, p,c;
    std::vector<A3> F; 
    void update_grid();
    void update();
};


#endif 