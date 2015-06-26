#ifndef MODEL_HPP
#define MODEL_HPP

#include <armadillo>
#include <array>
#include <map>
#include <string>

#include "constant.hpp"
#include "geometry.hpp"

class model {
public:
    // parameters
    std::string name;
    double E_g;
    double m_eff;
    double E_gc;
    double m_efc;
    std::array<double, 3> F;

    // hopping parameters
    double t1;       // hopping between orbitals in same unit cell
    double t2;       // hopping between orbitals in neighbouring unit cells
    double tc1;      // hopping between orbitals in same unit cell in contact area
    double tc2;      // hopping between orbitals in neighbouring unit cells in contact area
    double tcc;      // hopping between contact and central area
    arma::vec t_vec; // vector with t-values

    inline model(const geometry & g,
                 const std::string & n,
                 double E_g,
                 double m_eff,
                 double E_gc,
                 double m_efc,
                 const std::array<double, 3> & F);
};

static const model nfet_model (
    fet_geometry,     // geometry
    "NFET",           // name
    0.62,             // E_g
    0.01 * c::m_e,    // m_eff
    0.62,             // E_gc
    0.01 * c::m_e,    // m_efc
{   0.62 / 2 + 0.015, // F[S]
    0.00,             // F[G]
    0.62 / 2 + 0.015  // F[D]
});

static const model pfet_model (
    fet_geometry,     // geometry
    "PFET",           // name
    nfet_model.E_g,   // E_g
    nfet_model.m_eff, // m_eff
    nfet_model.E_gc,  // E_gc
    nfet_model.m_efc, // m_efc
{  -nfet_model.F[S],  // F[S]
   -nfet_model.F[G],  // F[G]
   -nfet_model.F[D]   // F[D]
});

static const model ntfet_model (
    tfet_geometry,    // geometry
    "NTFET",          // name
    0.62,             // E_g
    0.05 * c::m_e,    // m_eff
    0.62,             // E_gc
    0.05 * c::m_e,    // m_efc
{  -0.62 / 2 - 0.015, // F[S] (p++)
    0.00,             // F[G]
   +0.62 / 2 + 0.001  // F[D] (n+)
});

static const model ptfet_model (
    tfet_geometry,     // geometry
    "PTFET",           // name
    ntfet_model.E_g,   // E_g
    ntfet_model.m_eff, // m_eff
    ntfet_model.E_gc,  // E_gc
    ntfet_model.m_efc, // m_efc
{  -ntfet_model.F[S],  // F[S] (n++)
   -ntfet_model.F[G],  // F[G]
   -ntfet_model.F[D]   // F[D] (p+)
});

//----------------------------------------------------------------------------------------------------------------------

model::model(const geometry & g,
             const std::string & n,
             double E_g_,
             double m_eff_,
             double E_gc_,
             double m_efc_,
             const std::array<double, 3> & F_)
    : name(n),
      E_g(E_g_),
      m_eff(m_eff_),
      E_gc(E_gc_),
      m_efc(m_efc_),
      F(F_) {
    // hopping parameters
    t1  = 0.25 * E_g  * (1 + sqrt(1 + 2 * c::h_bar2 / (g.dx*g.dx * 1E-18 * m_eff * E_g  * c::e)));
    t2  = 0.25 * E_g  * (1 - sqrt(1 + 2 * c::h_bar2 / (g.dx*g.dx * 1E-18 * m_eff * E_g  * c::e)));
    tc1 = 0.25 * E_gc * (1 + sqrt(1 + 2 * c::h_bar2 / (g.dx*g.dx * 1E-18 * m_efc * E_gc * c::e)));
    tc2 = 0.25 * E_gc * (1 - sqrt(1 + 2 * c::h_bar2 / (g.dx*g.dx * 1E-18 * m_efc * E_gc * c::e)));
    tcc = 2.0 / (1.0 / t2 + 1.0 / tc2);

    // create t_vec
    t_vec = arma::vec(g.N_x * 2 - 1);
    bool b = true;
    for (unsigned i = g.sc2.a; i < g.sc2.b; ++i) {
        t_vec(i) = b ? tc1 : tc2;
        b = !b;
    }
    t_vec(g.sc2.b) = tcc;
    b = true;
    for (unsigned i = g.sox2.a; i < g.dox2.b; ++i) {
        t_vec(i) = b ? t1 : t2;
        b = !b;
    }
    t_vec(g.dox2.b) = tcc;
    b = true;
    for (unsigned i = g.dc2.a; i < g.dc2.b; ++i) {
        t_vec(i) = b ? tc1 : tc2;
        b = !b;
    }
}

#endif

