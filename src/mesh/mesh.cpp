
#include"mesh.hpp"

Mesh::Mesh(int nx_i):nx_c(nx_i),nx_f(nx_c+1),x_g(nx_f),x_c(nx_c),rho(nx_c),mome(nx_c),e(nx_c),u(nx_c),p(nx_c),c(nx_c),F(nx_f){}
void Mesh::update_grid()
{
    for(int i=0;i<nx_c;i++)
    {
        x_c[i]=(x_g[i]+x_g[i+1])*0.5;
    }
}

void Mesh::update()
{
    for(int i=0;i<nx_c;i++)
    {
        u[i]=mome[i]/rho[i];
        p[i]=(gamma-1)*(e[i]-0.5*u[i]*mome[i]);
        c[i]=sqrt(gamma*p[i]/rho[i]);
    }
}