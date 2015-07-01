#ifndef DEVICE_PARAMS_HPP
#define DEVICE_PARAMS_HPP

#include "geometry.hpp"
#include "model.hpp"

class device_params : public geometry, public model {
public:
    std::string name;

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

    // hopping parameters
    double t1;       // hopping between orbitals in same unit cell
    double t2;       // hopping between orbitals in neighbouring unit cells
    double tc1;      // hopping between orbitals in same unit cell in contact area
    double tc2;      // hopping between orbitals in neighbouring unit cells in contact area
    double tcc;      // hopping between contact and central area
    arma::vec t_vec; // vector with t-values

    inline device_params(const std::string & n, const geometry & g, const model & m);

    inline std::string to_string();
};

static const device_params nfet("nfet", fet_geometry, nfet_model);
static const device_params pfet("pfet", fet_geometry, pfet_model);
static const device_params ntfet("ntfet", tfet_geometry, ntfet_model);
static const device_params ptfet("ptfet", tfet_geometry, ptfet_model);

//----------------------------------------------------------------------------------------------------------------------

device_params::device_params(const std::string & n_, const geometry & g_, const model & m_)
    : geometry(g_), model(m_), name(n_) {
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

    // hopping parameters
    t1  = 0.25 * E_g  * (1 + sqrt(1 + 2 * c::h_bar2 / (dx*dx * 1E-18 * m_eff * E_g  * c::e)));
    t2  = 0.25 * E_g  * (1 - sqrt(1 + 2 * c::h_bar2 / (dx*dx * 1E-18 * m_eff * E_g  * c::e)));
    tc1 = 0.25 * E_gc * (1 + sqrt(1 + 2 * c::h_bar2 / (dx*dx * 1E-18 * m_efc * E_gc * c::e)));
    tc2 = 0.25 * E_gc * (1 - sqrt(1 + 2 * c::h_bar2 / (dx*dx * 1E-18 * m_efc * E_gc * c::e)));
    tcc = 2.0 / (1.0 / t2 + 1.0 / tc2);

    // create t_vec
    t_vec = arma::vec(N_x * 2 - 1);
    bool b = true;
    for (unsigned i = sc2.a; i < sc2.b; ++i) {
        t_vec(i) = b ? tc1 : tc2;
        b = !b;
    }
    t_vec(sc2.b) = tcc;
    b = true;
    for (unsigned i = sox2.a; i < dox2.b; ++i) {
        t_vec(i) = b ? t1 : t2;
        b = !b;
    }
    t_vec(dox2.b) = tcc;
    b = true;
    for (unsigned i = dc2.a; i < dc2.b; ++i) {
        t_vec(i) = b ? tc1 : tc2;
        b = !b;
    }
}

std::string device_params::to_string() {
    using namespace std;

    stringstream ss;

    ss << "name    = " << name    << endl;

    ss << endl << "; model" << endl;
    ss << model::to_string();

    ss << endl << "; geometry" << endl;
    ss << geometry::to_string();

    return ss.str();
}

#endif

