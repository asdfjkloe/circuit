#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <armadillo>

class geometry {
public:
    // parameters
    double eps_cnt;
    double eps_ox;
    double l_sc;
    double l_sox;
    double l_sg;
    double l_g;
    double l_dg;
    double l_dox;
    double l_dc;
    double r_cnt;
    double d_ox;
    double r_ext;
    double dx;
    double dr;
};

static const geometry fet_geometry {
    10.0, // eps_cnt
    25.0, // eps_ox
     5.0, // l_sc
     7.0, // l_sox
     3.0, // l_sg
    18.0, // l_g
     3.0, // l_dg
     7.0, // l_dox
     5.0, // l_dc
     1.0, // r_cnt
     2.0, // d_ox
     2.0, // r_ext
     0.2, // dx
     0.1  // dr
};

static const geometry tfet_geometry {
    10.0, // eps_cnt
    25.0, // eps_ox
     5.0, // l_sc
    22.0, // l_sox
     3.0, // l_sg
    10.0, // l_g
    25.0, // l_dg
     0.0, // l_dox
     5.0, // l_dc
     1.0, // r_cnt
     2.0, // d_ox
     2.0, // r_ext
     0.2, // dx
     0.1  // dr
};

#endif // GEOMETRY_HPP

