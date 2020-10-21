%module Euler1D


%{

#include"mesh.hpp"
#include"SWFS.hpp"
#include"MSWv1.hpp"
#include"MSWv2.hpp"

%}
%include"mesh.hpp"
%include"SWFS.hpp"
%include"MSWv1.hpp"
%include"MSWv2.hpp"
%include <std_vector.i>

%template(DoubleVector) std::vector<double>;
%template(DoubleA3Vector) std::vector<A3>;
%template(DoubleM3Vector) std::vector<M3>;