#ifndef DEVICE_HPP
#define DEVICE_HPP

#include <memory>

#include "constant.hpp"
#include "contact.hpp"
#include "current.hpp"
#include "geometry.hpp"
#include "potential.hpp"
#include "model.hpp"
#include "voltage.hpp"

class device {
public:
    static constexpr double dphi_threshold = 1e-9;
    static constexpr int max_iterations = 50;

    // geometry and model
    geometry g;
    model m;

    // observables
    potential phi;
    charge_density n;
    current I;

    // lattices and weights for adaptive energy integration
    arma::vec E0[4];
    arma::vec W[4];

    // pointers to contacts
    using contact_ptrs = std::array<contact_ptr, 3>;
    contact_ptrs contacts;

    // constructors
    inline device(const geometry & g, const model & m, const contact_ptrs & ct);
    inline device(const geometry & g, const model & m, const voltage & V);
    inline device(const geometry & g, const model & m);

    // solve steady state
    template<bool smooth = true>
    inline bool steady_state();

};

//----------------------------------------------------------------------------------------------------------------------

device::device(const geometry & gg, const model & mm, const contact_ptrs & ct)
    : g(gg), m(mm), contacts(ct) {
}

device::device(const geometry & gg, const model & mm, const voltage & V)
    : device(gg, mm,
             contact_ptrs {
                 std::make_shared<contact>(V[S], c::inf),
                 std::make_shared<contact>(V[G], c::inf),
                 std::make_shared<contact>(V[D], c::inf)
             }) {
}

device::device(const geometry & gg, const model & mm)
    : device(gg, mm, voltage { 0.0, 0.0, 0.0 }) {
}

template<bool smooth>
bool device::steady_state() {
    // get the right-hand-side vector in poisson's equation
    arma::vec R0 = potential::get_R0(g, m, { contacts[S]->V, contacts[G]->V, contacts[D]->V });

    // solve poisson's equation with zero charge_density
    phi = potential(g, R0);

    // maximum deviation of phi
    double dphi;

    // iteration counter
    int it;

    // was self-consistency achieved in max_iterations steps?
    bool converged = false;

    // init anderson
    anderson mr_neo(phi.data);

    // repeat until potential converged or maximum number of iterations has been reached
    for (it = 1; it <= max_iterations; ++it) {
        // update charge density
        n = charge_density(g, m, phi, E0, W);

        // update potential
        dphi = phi.update(g, R0, n, mr_neo);

        // check for convergence (i.e. deviation is smaller than threshold value)
        converged = dphi < dphi_threshold;
        if (converged) {
            break;
        }

        // straighten out the potential in the contact regions to improve convergence (can be turned off)
        if (smooth && (it < 3)) {
            phi.smooth(g, m);
        }
    }

    // get current
    I = current(g, m, phi);

    std::cout << it << " iterations, reldev=" << dphi/dphi_threshold << ", " << converged ? "converged!" : "DIVERGED!!!";
    std::cout << ", n_E = " << E0[0].size() + E0[1].size() + E0[2].size() + E0[3].size() << std::endl;

    return converged;
}

#endif

