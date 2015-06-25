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

