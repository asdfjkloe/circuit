#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <armadillo>
#include <string>

class geometry {
public:
    // parameters
    std::string name;
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

    // total lengths
    double l;
    double R;

    // x lattice
    int N_sc;        // # of points in source contact
    int N_sox;       // # of points in source oxide
    int N_sg;        // # of points between source and gate
    int N_g;         // # of points in gate
    int N_dg;        // # of points between drain and gate
    int N_dox;       // # of points in drain oxide
    int N_dc;        // # of points in drain contact
    int N_x;         // total # of points
    arma::vec x;     // x lattice points

    // r lattice
    int M_cnt;       // # of points in nanotube
    int M_ox;        // # of points in oxide
    int M_ext;       // # of points over oxide
    int M_r;         // total # of points
    arma::vec r;     // r lattice points

    // ranges
    arma::span sc;   // source contact area
    arma::span sox;  // source oxide area
    arma::span sg;   // area between source and gate
    arma::span g;    // gate area;
    arma::span dg;   // area between drain and gate
    arma::span dox;  // drain oxide area
    arma::span dc;   // drain contact area
    arma::span sc2;  // source contact area twice
    arma::span sox2; // source oxide area twice;
    arma::span sg2;  // area between source and gate
    arma::span g2;   // gate area twice
    arma::span dg2;  // area between drain and gate twice
    arma::span dox2; // drain oxide area twice
    arma::span dc2;  // drain contact twice

    inline geometry(const std::string & n,
                    double eps_cnt,
                    double eps_ox,
                    double l_sc,
                    double l_sox,
                    double l_sg,
                    double l_g,
                    double l_dg,
                    double l_dox,
                    double l_dc,
                    double r_cnt,
                    double d_ox,
                    double r_ext,
                    double dx,
                    double dr);
};

static const geometry fet_geometry(
    "FET", // name
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
);

static const geometry tfet_geometry(
    "TFET", // name
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
);

//----------------------------------------------------------------------------------------------------------------------

geometry::geometry(const std::string & n,
                   double eps_cnt_,
                   double eps_ox_,
                   double l_sc_,
                   double l_sox_,
                   double l_sg_,
                   double l_g_,
                   double l_dg_,
                   double l_dox_,
                   double l_dc_,
                   double r_cnt_,
                   double d_ox_,
                   double r_ext_,
                   double dx_,
                   double dr_)
    : name(n),
      eps_cnt(eps_cnt_),
      eps_ox(eps_ox_),
      l_sc(l_sc_),
      l_sox(l_sox_),
      l_sg(l_sg_),
      l_g(l_g_),
      l_dg(l_dg_),
      l_dox(l_dox_),
      l_dc(l_dc_),
      r_cnt(r_cnt_),
      d_ox(d_ox_),
      r_ext(r_ext_),
      dx(dx_),
      dr(dr_) {
    // total lengths
    l = l_sc + l_sox + l_sg + l_g + l_dg + l_dox + l_dc;
    R = r_cnt + d_ox + r_ext;

    // x lattice
    N_sc  = round(l_sc  / dx);
    N_sox = round(l_sox / dx);
    N_sg  = round(l_sg  / dx);
    N_g   = round(l_g   / dx);
    N_dg  = round(l_dg  / dx);
    N_dox = round(l_dox / dx);
    N_dc  = round(l_dc  / dx);
    N_x   = N_sc + N_sox + N_sg + N_g + N_dg + N_dox + N_dc;
    x     = arma::linspace(0.5 * dx, l - 0.5 * dx, N_x);

    // r lattice
    M_cnt = round(r_cnt / dr);
    M_ox  = round(d_ox  / dr);
    M_ext = round(r_ext / dr);
    M_r   = M_cnt + M_ox + M_ext;
    r     = arma::linspace(0.5 * dr, R - 0.5 * dr, M_r);

    // ranges
    sc   = arma::span(        0,   - 1 + N_sc );
    sox  = arma::span( sc.b + 1,  sc.b + N_sox);
    sg   = arma::span(sox.b + 1, sox.b + N_sg );
    g    = arma::span( sg.b + 1,  sg.b + N_g  );
    dg   = arma::span(  g.b + 1,   g.b + N_dg );
    dox  = arma::span( dg.b + 1,  dg.b + N_dox);
    dc   = arma::span(dox.b + 1, dox.b + N_dc );
    sc2  = arma::span( sc.a * 2,  sc.b * 2 + 1);
    sox2 = arma::span(sox.a * 2, sox.b * 2 + 1);
    sg2  = arma::span( sg.a * 2,  sg.b * 2 + 1);
    g2   = arma::span(  g.a * 2,   g.b * 2 + 1);
    dg2  = arma::span( dg.a * 2,  dg.b * 2 + 1);
    dox2 = arma::span(dox.a * 2, dox.b * 2 + 1);
    dc2  = arma::span( dc.a * 2,  dc.b * 2 + 1);
}

#endif // GEOMETRY_HPP

