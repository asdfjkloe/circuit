#ifndef CHARGE_DENSITY_HPP_HEADER
#define CHARGE_DENSITY_HPP_HEADER

#include <armadillo>

#include "geometry.hpp"
#include "model.hpp"
#include "util/fermi.hpp"
#include "util/integral.hpp"

// forward declarations
class potential;
class wave_packet;

class charge_density {
public:
    static constexpr double E_min = -1.2;
    static constexpr double E_max = 1.2;
    static constexpr double rel_tol = 8e-3;
    static constexpr int initial_waypoints = 30;

    arma::vec lv;    // density due to left valence band
    arma::vec rv;    // density due to right valence band
    arma::vec lc;    // density due to left conduction band
    arma::vec rc;    // density due to right conduction band
    arma::vec total; // total density

    inline charge_density();
    inline charge_density(const geometry & g, const model & m, const potential & phi, arma::vec E0[4], arma::vec W[4]);
    inline charge_density(const geometry & g, const model & m, const potential & phi, const wave_packet psi[4]);

private:
    static inline arma::vec get_bound_states(const geometry & g, const model & m, const potential & phi);
    static inline arma::vec get_bound_states(const geometry & g, const model & m, const potential & phi, double E0, double E1);
    static inline arma::vec get_intervals(const arma::vec & E_bound, double E0, double E1);

    template<bool source>
    static inline arma::vec get_A(const geometry & g, const model & m, const potential & phi, double E);

    static inline arma::vec get_n0(const geometry & g, const model & m);
};

#endif

//----------------------------------------------------------------------------------------------------------------------

#ifndef CHARGE_DENSITY_HPP_BODY
#define CHARGE_DENSITY_HPP_BODY

charge_density::charge_density() {
}

charge_density::charge_density(const geometry & g, const model & m, const potential & phi, arma::vec E0[4], arma::vec W[4]) {
    using namespace arma;

    // get a prediction of bound-state levels in the channel region
    vec E_bound = get_bound_states(g, m, phi);

    // get integration intervals
    vec i_lv = get_intervals(E_bound, phi.s() + E_min, phi.s() - 0.5 * m.E_gc);
    vec i_rv = get_intervals(E_bound, phi.d() + E_min, phi.d() - 0.5 * m.E_gc);
    vec i_lc = get_intervals(E_bound, phi.s() + 0.5 * m.E_gc, phi.s() + E_max);
    vec i_rc = get_intervals(E_bound, phi.d() + 0.5 * m.E_gc, phi.d() + E_max);

    // integrand function definitions
    auto I_l = [&] (double E) -> vec {
        // spectral function from source
        vec A = get_A<true>(g, m, phi, E);

        // fermi distribution in source
        double f = fermi(E - phi.s(), m.F[S]);

        // multiply with f or (f - 1) depending on branching point
        for (int i = 0; i < g.N_x; ++i) {
            A(i) *= fermi<true>(f, E - phi(i));
        }

        return A;
    };
    auto I_r = [&] (double E) -> vec {
        // spectral function from drain
        vec A = get_A<true>(g, m, phi, E);

        // fermi distribution in drain
        double f = fermi(E - phi.d(), m.F[D]);

        // multiply with f or (f - 1) depending on branching point
        for (int i = 0; i < g.N_x; ++i) {
            A(i) *= fermi<true>(f, E - phi(i));
        }

        return A;
    };

    // scaling
    double scale = - 0.5 * c::e / M_PI / M_PI / g.r_cnt / g.dr / g.dx;

    // calculate integrals (adaptively) and output the used energy lattice and weights
    lv = scale * integral(I_l, g.N_x, i_lv, rel_tol, c::epsilon(), E0[LV], W[LV]);
    rv = scale * integral(I_r, g.N_x, i_rv, rel_tol, c::epsilon(), E0[RV], W[RV]);
    lc = scale * integral(I_l, g.N_x, i_lc, rel_tol, c::epsilon(), E0[LC], W[LC]);
    rc = scale * integral(I_r, g.N_x, i_rc, rel_tol, c::epsilon(), E0[RC], W[RC]);

    // calculate total charge density and add the background charge introduced by dopands
    total = lv + rv + lc + rc + get_n0(g, m);
}

arma::vec charge_density::get_bound_states(const geometry & g, const model & m, const potential & phi) {
    double phi0, phi1, phi2, limit;

    // check for bound states in valence band
    phi0 = arma::min(phi.data(g.sg)) - 0.5 * m.E_g;
    phi1 = arma::max(phi.data(g.g )) - 0.5 * m.E_g;
    phi2 = arma::min(phi.data(g.dg)) - 0.5 * m.E_g;
    limit = phi0 > phi2 ? phi0 : phi2;
    if (limit < phi1) {
        return get_bound_states(g, m, phi, limit, phi1);
    }

    // check for bound states in conduction band
    phi0 = arma::max(phi.data(g.sg)) + 0.5 * m.E_g;
    phi1 = arma::min(phi.data(g.g )) + 0.5 * m.E_g;
    phi2 = arma::max(phi.data(g.dg)) + 0.5 * m.E_g;
    limit = phi0 < phi2 ? phi0 : phi2;
    if (limit > phi1) {
        return get_bound_states(g, m, phi, phi1, limit);
    }

    // no bound states found
    return arma::vec(arma::uword(0));
}

arma::vec charge_density::get_bound_states(const geometry & g, const model & m, const potential & phi, double E0, double E1) {
    using namespace arma;

    static constexpr double tol = 1e-10;

    // only look at the gated region
    span range { g.sg2.a, g.dg2.b };

    // for building the hamiltonian only in the gated region
    vec a = m.t_vec(range);   // off-diagonal
    vec a2 = a % a;           // off-diagonal squared
    vec b = phi.twice(range); // main diagonal

    // evaluate Sturm sequence to find number of eigenvalues smaller than E
    auto eval = [&] (double E) {
        int n = b.size();

        static const double eps = c::epsilon();

        // first iteration (i = 0)
        double q;
        double q0 = b[0] - E;
        int s = q0 < 0 ? 1 : 0;

        // start with i = 1
        for (int i = 1; i < n; ++i) {
            if (q0 == 0) {
                q = b[i] - E - a[i - 1] / eps;
            } else {
                q = b[i] - E - a2[i - 1] / q0;
            }

            q0 = q;
            if (q < 0) {
                ++s;
            }
        }

        return s;
    };

    double E2;
    int i0, i1;
    int s0, s1, s2;

    // s1 - s0 is the number of eigenvalues in the interval
    s0 = eval(a, a2, b, E0);
    s1 = eval(a, a2, b, E1);

    // check if no bound states in this interval
    if (s1 - s0 == 0) {
        return vec(uword(0));
    }

    // start iterative bisection method to find all eigenvalues in the interval
    unsigned n = 2;
    vec E(1025);
    ivec s(1025);
    E(0) = E0;
    E(1) = E1;
    s(0) = s0;
    s(1) = s1;

    unsigned n_bound = 0;
    vec E_bound(100);

    // stack for recursion
    std::stack<std::pair<int, int>> stack;

    // push first interval to stack
    stack.push(std::make_pair(0, 1));

    // repeat until all intervals inspected
    while (!stack.empty()) {
        const auto & i = stack.top();
        i0 = i.first;
        i1 = i.second;

        stack.pop();

        // load data
        E0 = E(i0);
        E1 = E(i1);
        s0 = s(i0);
        s1 = s(i1);

        // mid energy
        E2 = 0.5 * (E0 + E1);

        // if interval size sufficiently small enough, add new bound state (can be degenerate eigenvalue)
        if (E1 - E0 <= tol) {
            if (E_bound.size() <= n_bound) {
                E_bound.resize(n_bound * 2);
            }
            E_bound(n_bound++) = E2;
        } else {
            // evaluate s at mid energy
            s2 = eval(a, a2, b, E2);

            // add intervals to stack if they contain bound states
            if (s1 - s2 > 0) {
                stack.push(std::make_pair(n, i1));
            }
            if (s2 - s0 > 0) {
                stack.push(std::make_pair(i0, n));
            }

            // save E2 and s2
            if (E.size() <= n) {
                E.resize(2 * n - 1);
                s.resize(2 * n - 1);
            }
            E(n) = E2;
            s(n) = s2;
            ++n;
        }
    }

    // shrink to fit and return bound states
    E_bound.resize(n_bound);
    return E_bound;
}

#endif
