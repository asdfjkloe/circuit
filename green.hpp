#ifndef GREEN_HPP
#define GREEN_HPP

#include <armadillo>

#include "geometry.hpp"
#include "model.hpp"
#include "potential.hpp"
#include "util/gnuplot.hpp"
#include "util/inverse.hpp"

static inline void self_energy(const model & m, const potential & phi, double E, arma::cx_double & Sigma_s, arma::cx_double & Sigma_d) {
    using namespace arma;
    using namespace std;

    // kinetic energy in source and drain
    auto E_s = E - phi.s();
    auto E_d = E - phi.d();

    // shortcuts
    double t12 = m.tc1 * m.tc1;
    double t22 = m.tc2 * m.tc2;

    // self energy
    Sigma_s = E_s * E_s - t12 - t22;
    Sigma_s = 0.5 * (E_s * E_s - t12 + t22 + sqrt(Sigma_s * Sigma_s + - 4 * t12 * t22)) / E_s;
    Sigma_d = E_d * E_d - t12 - t22;
    Sigma_d = 0.5 * (E_d * E_d - t12 + t22 + sqrt(Sigma_d * Sigma_d + - 4 * t12 * t22)) / E_d;

    // imaginary part must be negative
    Sigma_s.imag(-std::abs(Sigma_s.imag()));
    Sigma_d.imag(-std::abs(Sigma_d.imag()));
}

template<bool source>
static inline arma::cx_vec green_col(const model & m, const potential & phi, double E, arma::cx_double & Sigma_s, arma::cx_double & Sigma_d) {
    using namespace arma;

    self_energy(m, phi, E, Sigma_s, Sigma_d);

    // build diagonal part of hamiltonian
    auto D = conv_to<cx_vec>::from(E - phi.twice);
    D(0)            -= Sigma_s;
    D(D.size() - 1) -= Sigma_d;

    return inverse_col<source>(vec(-m.t_vec), D);
}

static inline arma::mat get_lDOS(const geometry & g, const model & m, const potential & phi, int N_grid, arma::vec & E) {
    using namespace arma;
    using namespace std::complex_literals;

    mat ret(N_grid, g.N_x);

    double phi_min = min(phi.data);
    double phi_max = max(phi.data);

    E = linspace(phi_min - 0.5 * m.E_g - 0.2, phi_max + 0.5 * m.E_g + 0.2, N_grid);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N_grid; ++i) {
        cx_double Sigma_s;
        cx_double Sigma_d;
        self_energy(m, phi, E(i), Sigma_s, Sigma_d);

        auto D = conv_to<cx_vec>::from(E(i) - phi.twice);
        D(0)            -= Sigma_s;
        D(D.size() - 1) -= Sigma_d;
        D += 0.001i;

        vec mixed = -arma::imag(inverse_diag(vec(-m.t_vec), D)) / M_PI;

        for (int j = 0; j < g.N_x; ++j) {
            ret(i, j) = mixed(2*j) + mixed(2*j+1);
        }
    }
    return ret;
}

static inline void plot_ldos(const geometry & g, const model & m, const potential & phi, const unsigned N_grid = 1000) {
    gnuplot gp;

    gp << "set title \"Logarithmic lDOS\"\n";
    gp << "set xlabel \"x / nm\"\n";
    gp << "set ylabel \"E / eV\"\n";
    gp << "set zlabel \"log(lDOS)\"\n";
    gp << "unset key\n";
    gp << "unset colorbox\n";
//    gp << "set terminal pdf rounded color enhanced font 'arial,12'\n";
//    gp << "set output 'lDOS.pdf'\n";

    arma::vec E;
    arma::mat lDOS = get_lDOS(g, m, phi, N_grid, E);
    gp.set_background(g.x, E, arma::log(lDOS));

    arma::vec vband = phi.data;
    vband(g.sc)  -= 0.5 * m.E_gc;
    vband({g.sox.a, g.dox.b}) -= 0.5 * m.E_g;
    vband(g.dc)  -= 0.5 * m.E_gc;

    arma::vec cband = phi.data;
    cband(g.sc)  += 0.5 * m.E_gc;
    cband({g.sox.a, g.dox.b}) += 0.5 * m.E_g;
    cband(g.dc)  += 0.5 * m.E_gc;

    gp.add(g.x, vband);
    gp.add(g.x, cband);
    gp.add(g.x, phi.data);

    unsigned N_s = std::round(g.N_sc + 0.5 * g.N_sox);
    arma::vec fermi_l(N_s);
    fermi_l.fill(m.F[S] + phi.s());
    arma::vec x_l = g.x(arma::span(0, N_s-1));
    gp.add(x_l, fermi_l);

    unsigned N_d = std::round(g.N_dc + 0.5 * g.N_dox);
    arma::vec fermi_r(N_d);
    fermi_r.fill(m.F[D] + phi.d());
    arma::vec x_r = g.x(arma::span(g.N_x-N_d, g.N_x-1));
    gp.add(x_r, fermi_r);

    gp << "set style line 1 lt 1 lc rgb RWTH_Orange lw 2\n";
    gp << "set style line 2 lt 1 lc rgb RWTH_Orange lw 2\n";
    gp << "set style line 3 lt 1 lc rgb RWTH_Rot lw 2\n";
    gp << "set style line 4 lc rgb RWTH_Schwarz lw 1 lt 3\n";
    gp << "set style line 5 lc rgb RWTH_Schwarz lw 1 lt 3\n";

    gp.plot();
}

#endif
