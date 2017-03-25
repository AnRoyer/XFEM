#ifndef struct_h
#define struct_h

#include <complex>

typedef struct Param{
    double k_0;
    double k_1;
    double x_bnd;
    std::complex<double> wave;
}Param;

typedef struct Physical{
    std::vector<int> tagDir;
    std::vector<int> tagInf;
}Physical;

#endif /* struct_h */
