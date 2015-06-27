#ifndef DEVICE_HPP
#define DEVICE_HPP

#include <memory>

#include "constant.hpp"
#include "contact.hpp"
#include "current.hpp"
#include "device_params.hpp"
#include "potential.hpp"
#include "voltage.hpp"

class device {
public:
    static constexpr double dphi_threshold = 1e-9; // convergence threshold
    static constexpr int max_iterations = 50;      // maximum number of iterations in convergence
    static constexpr unsigned mem = 2000;          // length of memory integral

    // parameters
    device_params p;

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

    // waves for time evolution
    wave_packet psi[4];

    // constructors
    inline device(const device_params & p, const contact_ptrs & ct);
    inline device(const device_params & p, const voltage & V);
    inline device(const device_params & p);

    // solve steady state
    template<bool smooth = true>
    inline bool steady_state();

    // initialize time evolution (call after steady_state)
    inline void init_time_evolution();

    // simulate next time step
    inline void time_step();

private:
    arma::cx_mat u;
    arma::cx_mat L;
    arma::cx_mat q;
    arma::cx_mat qsum;
    arma::cx_mat H_eff;
    arma::cx_mat old_L;

    inline void calc_q();
};

//----------------------------------------------------------------------------------------------------------------------

device::device(const device_params & pp, const contact_ptrs & ct)
    : p(pp), contacts(ct) {
}

device::device(const device_params & pp, const voltage & V)
    : device(pp,
             contact_ptrs {
                 std::make_shared<contact>(V[S], c::inf),
                 std::make_shared<contact>(V[G], c::inf),
                 std::make_shared<contact>(V[D], c::inf)
             }) {
}

device::device(const device_params & pp)
    : device(pp, voltage { 0.0, 0.0, 0.0 }) {
}

template<bool smooth>
bool device::steady_state() {
    // get the right-hand-side vector in poisson's equation
    arma::vec R0 = potential::get_R0(p, { contacts[S]->V, contacts[G]->V, contacts[D]->V });

    // solve poisson's equation with zero charge_density
    phi = potential(p, R0);

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
        n = charge_density(p, phi, E0, W);

        // update potential
        dphi = phi.update(p, R0, n, mr_neo);

        // check for convergence (i.e. deviation is smaller than threshold value)
        converged = dphi < dphi_threshold;
        if (converged) {
            break;
        }

        // straighten out the potential in the contact regions to improve convergence (can be turned off)
        if (smooth && (it < 3)) {
            phi.smooth(p);
        }
    }

    // get current
    I = current(p, phi);

    std::cout << "(" << p.name << ") steady_state: " <<  it << " iterations, reldev=" << dphi/dphi_threshold;
    std::cout << ", " << (converged ? "converged!" : "DIVERGED!!!");
    std::cout << ", n_E = " << E0[0].size() + E0[1].size() + E0[2].size() + E0[3].size() << std::endl;

    return converged;
}

void device::init_time_evolution() {
    using namespace arma;

    // initialize waves
    psi[LV] = wave_packet(p, S, mem, E0[LV], W[LV], phi);
    psi[RV] = wave_packet(p, D, mem, E0[RV], W[RV], phi);
    psi[LC] = wave_packet(p, S, mem, E0[LC], W[LC], phi);
    psi[RC] = wave_packet(p, D, mem, E0[RC], W[RC], phi);

    // precalculate q values
    calc_q();

    // build constant part of Hamiltonian
    H_eff = cx_mat(2 * p.N_x, 2 * p.N_x);
    H_eff.fill(0);
    H_eff.diag(+1) = conv_to<cx_vec>::from(p.t_vec);
    H_eff.diag(-1) = conv_to<cx_vec>::from(p.t_vec);

    // setup L
    L.fill(1.0);
}

void device::time_step() {

}

void device::calc_q() {
    using namespace arma;
    using namespace std::complex_literals;

    std::cout << "(" << p.name << ") calculating q-values... ";
    std::flush(std::cout);

    // shortcuts
    static constexpr double g = c::g;
    static constexpr double g2 = g * g;
    const double t1 = p.tc1;
    const double t12 = t1 * t1;
    const double t2 = p.tc2;
    const double t22 = t2 * t2;
    static const cx_mat22 eye2 = { 1, 0, 0, 1 };

    // get q values dependent on potential in lead
    auto get_q = [&] (double phi0) {
        // storage
        cx_vec qq(mem + 4);
        std::vector<cx_mat22> P(mem + 4);

        // hamiltonian in lead
        mat22 h = { phi0, t1, t1, phi0 };

        // coupling hamiltonian
        mat22 Vau = { 0, t2, 0, 0 };

        // set first 3 values of q and p to 0
        std::fill(qq.begin(), qq.begin() + 3, 0);
        std::fill(P.begin(), P.begin() + 3, cx_mat22{ 0, 0, 0, 0 });

        // first actual q value (with pq-formula)
        auto a = (1.0 + 2i * g * phi0 + g2 * (t12 - t22 - phi0*phi0)) / g2 / (1.0 + 1i * g * phi0);
        qq(3) = - 0.5 * a + sqrt(0.25 * a * a + t22 / g2);

        // first actual p value
        P[3] = inv(eye2 + 1i * g * h + cx_mat22{ g2 * qq(3), 0, 0, 0 });

        // calculate A & C parameters
        cx_mat22 A = eye2 + 1i * g * h + g2 * Vau.t() * P[3] * Vau;
        auto C = A(0,0) * A(1,1) - A(0,1) * A(1,0);

        // loop over all time steps
        for (uword i = 4; i < mem + 4; ++i) {
            // perform sum
            cx_mat22 R = { 0, 0, 0, 0 };
            for (uword k = 4; k < i; ++k) {
                R += (P[k] + 2 * P[k - 1] + P[k - 2]) * Vau * P[i - k + 3];
            }

            // calculate B parameter
            cx_mat22 B = (eye2 - 1i * g * h) * P[i - 1] - g2 * Vau.t() * ((2 * P[i - 1] + P[i - 2]) * Vau * P[3] + R);

            // calculate next p values
            P[i](1,1) = (A(1,0) * B(0,1) - A(0,0) * B(1,1)) / (g2 * t22 * P[3](0,1) * A(1,0) - C);
            P[i](0,1) = (B(1,1) - A(1,1) * P[i](1,1)) / A(1,0);
            P[i](0,0) = (A(1,1) * B(0,0) - A(0,1) * B(1,0) - g2 * t22 * P[3](0,0) * P[i](1,1) * A(1,1)) / C;
            P[i](1,0) = (B(1,0) - A(1,0) * P[i](0,0)) / A(1,1);

            // calculate next q value
            qq(i) = t22 * P[i](1,1);
        }

        return qq;
    };

    // calculate and save q values
    q.resize(mem + 1, 2);
    q.col(S) = get_q(phi.s()).rows(3, mem + 3);
    q.col(D) = get_q(phi.d()).rows(3, mem + 3);

    // sum of two following q-values reversed
    qsum.resize(mem, 2);
    for (unsigned i = 0; i < mem; ++i) {
        qsum(i, S) = q(mem - i, S) + q(mem - 1 - i, S);
        qsum(i, D) = q(mem - i, D) + q(mem - 1 - i, D);
    }
    std::cout << "done!" << std::endl;
}

#endif

