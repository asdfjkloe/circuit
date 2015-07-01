#ifndef DEVICE_HPP
#define DEVICE_HPP

#include <iomanip>
#include <memory>
#include <string>
#include <iostream>

#include "constant.hpp"
#include "contact.hpp"
#include "current.hpp"
#include "device_params.hpp"
#include "potential.hpp"
#include "voltage.hpp"
#include "util/inverse.hpp"
#include "util/system.hpp"
#include "util/gnuplot.hpp"

class device {
public:
    static constexpr double dphi_threshold = 1e-5; // convergence threshold
    static constexpr int max_iterations = 50;      // maximum number of iterations before abortion
    static constexpr unsigned mem = 2000;          // maximum length of the memory integral

    // name
    std::string name;

    // parameters
    device_params p;

    // current time step
    unsigned m;

    // voltage history (not read by time_step)
    std::vector<voltage<3>> V;

    // observables
    std::vector<potential> phi;
    std::vector<charge_density> n;
    std::vector<current> I;

    // lattices and weights for adaptive integration of the charge-density
    arma::vec E0[4];
    arma::vec W[4];

    // pointers to contact-objects
    using contact_ptrs = std::array<contact_ptr, 3>;
    contact_ptrs contacts;

    // wavefunctions used in time evolutions
    wave_packet psi[4];

    // constructors
    inline device(const std::string & n, const device_params & p, const contact_ptrs & ct);
    inline device(const std::string & n, const device_params & p, const voltage<3> & V);
    inline device(const std::string & n, const device_params & p);

    // solve steady state
    template<bool smooth = true>
    inline bool steady_state();

    // initialize time evolution (only call after steady state was solved!)
    inline void init_time_evolution(int N_t);

    // propagate wavefunctions to the next time step
    inline bool time_step();

    // update contacts (current)
    inline void update_contacts();

    template<bool plots = false>
    inline void save();
private:
    // variables used by the time-evolution algorithm
    arma::cx_mat u;
    arma::cx_mat L;
    arma::cx_mat q;
    arma::cx_mat qsum;
    arma::cx_mat H_eff;
    arma::cx_mat old_L;
    arma::cx_mat cx_eye;

    // precalculate the (solely geometry-dependent) q-values
    inline void calc_q();
};

template<bool csv = false>
static inline arma::mat transfer(const device_params & p, const arma::vec & V_d, double V_g0, double V_g1, int N);

//----------------------------------------------------------------------------------------------------------------------

device::device(const std::string & n, const device_params & pp, const contact_ptrs & ct)
    : name(n), p(pp), m(0), contacts(ct) {
}

device::device(const std::string & n, const device_params & pp, const voltage<3> & V)
    : device(n, pp,
             contact_ptrs {
                 std::make_shared<contact>(V[S], c::inf),
                 std::make_shared<contact>(V[D], c::inf),
                 std::make_shared<contact>(V[G], c::inf)
             }) {
}

device::device(const std::string & n, const device_params & pp)
    : device(n, pp, voltage<3> { 0.0, 0.0, 0.0 }) {
}

template<bool smooth>
bool device::steady_state() {
    // initialize observables
    phi.resize(1);
    n.resize(1);
    I.resize(1);
    V.resize(1);

    // save voltage
    V[0] = { contacts[S]-> V, contacts[D]->V, contacts[G]->V };

    // get the right-hand-side vector in poisson's equation
    arma::vec R0 = potential::get_R0(p, V[0]);

    // solve poisson's equation with zero charge_density
    phi[0] = potential(p, R0);

    // maximum deviation of phi
    double dphi;

    // iteration counter
    int it;

    // was self-consistency achieved in max_iterations steps?
    bool converged = false;

    // init anderson
    anderson mr_neo(phi[0].data);

    // repeat until potential converged or maximum number of iterations has been reached
    for (it = 1; it <= max_iterations; ++it) {
        // update charge density
        n[0] = charge_density(p, phi[0], E0, W);

        // update potential
        dphi = phi[0].update(p, R0, n[0], mr_neo);

        // check for convergence (i.e. deviation is smaller than threshold value)
        converged = dphi < dphi_threshold;
        if (converged) {
            break;
        }

        // straighten out the potential in the contact regions to improve convergence (can be turned off)
        if (smooth && (it < 3)) {
            phi[0].smooth(p);
        }
    }

    // get current
    I[0] = current(p, phi[0]);

    std::cout << "(" << name << ") steady_state: " <<  it << " iterations, reldev=" << dphi/dphi_threshold;
    std::cout << ", " << (converged ? "" : "DIVERGED!!!");
    std::cout << ", n_E = " << E0[0].size() + E0[1].size() + E0[2].size() + E0[3].size() << std::endl;

    return converged;
}

void device::init_time_evolution(int N_t) {
    using namespace arma;

    m = 0;

    // resize observables
    phi.resize(N_t);
    n.resize(N_t);
    I.resize(N_t);

    // initialize waves
    psi[LV] = wave_packet(p, S, mem, E0[LV], W[LV], phi[0]);
    psi[RV] = wave_packet(p, D, mem, E0[RV], W[RV], phi[0]);
    psi[LC] = wave_packet(p, S, mem, E0[LC], W[LC], phi[0]);
    psi[RC] = wave_packet(p, D, mem, E0[RC], W[RC], phi[0]);

    // precalculate q values
    calc_q();

    // build constant part of Hamiltonian
    H_eff = cx_mat(2 * p.N_x, 2 * p.N_x);
    H_eff.fill(0);
    H_eff.diag(+1) = conv_to<cx_vec>::from(p.t_vec);
    H_eff.diag(-1) = conv_to<cx_vec>::from(p.t_vec);

    // setup u
    u.resize(N_t, 2);

    // setup L
    L.resize(N_t, 2);
    L.fill(1.0);

    // complex unity matrix
    cx_eye = eye<cx_mat>(2 * p.N_x, 2 * p.N_x);
}

bool device::time_step() {
    using namespace arma;
    using namespace std::complex_literals;

    // shortcut
    static constexpr double g = c::g;

    // next time step
    ++m;

    // get voltage and save it in history
    V[m] = { contacts[S]->V, contacts[D]->V, contacts[G]->V };

    // estimate charge distribution in following step from previous values
    n[m].total = (m == 1) ? n[0].total : (2 * n[m - 1].total - n[m - 2].total);

    // prepare right side of poisson equation
    vec R0 = potential::get_R0(p, V[m]);

    // first guess for potential
    phi[m] = potential(p, R0, n[m]);

    // maximum deviation of phi
    double dphi;

    // iteration counter
    int it;

    // was self-consistency achieved in max_iterations steps?
    bool converged = false;

    // init anderson
    anderson mr_neo(phi[m].data);

    // current data becomes old data
    for (int i = 0; i < 4; ++i) {
        psi[i].remember();
    }
    old_L = L;

    // self-consistency loop
    for (it = 0; it < max_iterations; ++it) {
        // diagonal of H with self-energy
        H_eff.diag() = conv_to<cx_vec>::from(0.5 * (phi[m].twice + phi[m - 1].twice));
        H_eff(            0,             0) -= 1i * g * q(0, S);
        H_eff(2 * p.N_x - 1, 2 * p.N_x - 1) -= 1i * g * q(0, D);

        // crank-nicolson propagator
        cx_mat U_eff = solve(cx_eye + 1i * g * H_eff, cx_eye - 1i * g * H_eff);

        // inv
        cx_mat inv(2 * p.N_x, 2);
        inv.col(S) = inverse_col< true>(cx_vec(1i * g * p.t_vec), cx_vec(1.0 + 1i * g * H_eff.diag()));
        inv.col(D) = inverse_col<false>(cx_vec(1i * g * p.t_vec), cx_vec(1.0 + 1i * g * H_eff.diag()));

        // u
        u(m - 1, S) = 0.5 * (phi[m].s() + phi[m - 1].s()) - phi[0].s();
        u(m - 1, D) = 0.5 * (phi[m].d() + phi[m - 1].d()) - phi[0].d();
        u.row(m - 1) = (1.0 - 0.5i * g * u.row(m - 1)) / (1.0 + 0.5i * g * u.row(m - 1));

        // L
        L.rows(0, m - 1) = old_L.rows(0, m - 1) * diagmat(u.row(m - 1) % u.row(m - 1));

        if (m == 1) {
            for (int i = 0; i < 4; ++i) {
                psi[i].memory_init();
                psi[i].source_init(p, u, q);
                psi[i].propagate(U_eff, inv);
                psi[i].update_E(p, phi[m], phi[0]);
            }
        } else {
            cx_mat affe;
            if (m <= mem + 1) {
                affe = - g * g * L.rows(0, m - 2) % qsum.rows(mem + 1 - m, mem - 1) / u.rows(0, m - 2) * diagmat(1.0 / u.row(m - 1));
            } else {
                affe = - g * g * L.rows(m - mem - 1, m - 2) % qsum.rows(0, mem - 1) / u.rows(m - mem - 1, m - 2) * diagmat(1.0 / u.row(m - 1));
            }

            // propagate wave functions
            for (int i = 0; i < 4; ++i) {
                psi[i].memory_update(affe, m);
                psi[i].source_update(u, L, qsum, m);
                psi[i].propagate(U_eff, inv);
                psi[i].update_E(p, phi[m], phi[0]);
            }
        }

        // update n
        n[m] = charge_density(p, phi[m], psi);

        // update potential
        dphi = phi[m].update(p, R0, n[m], mr_neo);

        // check for convergence (i.e. deviation is smaller than threshold value)
        converged = dphi < dphi_threshold;
        if (converged) {
            break;
        }
    }

    // update psi_sum
    for (int i = 0; i < 4; ++i) {
        psi[i].update_sum(m);
    }

    // get current
    I[m] = current(p, phi[m], psi);

    std::cout << "(" << name << ") timestep " << m << ": t=" << std::setprecision(5) << std::fixed << m * c::dt * 1e12
              << "ps, " << it + 1 << " iterations, reldev=" << dphi / dphi_threshold
              << ", " << (converged ? "" : "DIVERGED!!!") << std::endl;

    return converged;
}

void device::update_contacts() {
    contacts[S]->update(I[m].s(), c::dt);
    contacts[D]->update(I[m].d(), c::dt);
}

template<bool plots>
void device::save() {

    std::cout << "(" << name << ") saving time-dependent observables... ";
    std::flush(std::cout);

    int N_t = m + 1;
    arma::vec t = arma::linspace(0, N_t * c::dt, N_t);

    arma::mat phi_mat(p.N_x, N_t);
    arma::mat n_mat(p.N_x, N_t);
    arma::mat I_mat(p.N_x, N_t);
    arma::mat V_mat(N_t, 3);

    for (unsigned i = 0; i < N_t; ++i) {
        phi_mat.col(i) = phi[i].data;
        n_mat.col(i) = n[i].total;
        I_mat.col(i) = I[i].total;
        V_mat(i, S) = V[i][S];
        V_mat(i, D) = V[i][D];
        V_mat(i, G) = V[i][G];
    }

    std::string subfolder(save_folder() + "/" + name);
    system("mkdir -p " + subfolder);

    phi_mat.save(subfolder + "/phi.arma");
    n_mat.save(subfolder + "/n.arma");
    I_mat.save(subfolder + "/I.arma");
    p.x.save(subfolder + "/xtics.arma");
    t.save(subfolder + "/ttics.arma");
    V_mat.save(subfolder + "/V.arma");

    std::ofstream params_file(subfolder + "/params.ini");
    params_file << p.to_string();
    params_file.close();

    std::ofstream capacitance_file(subfolder + "/capacitance.ini");
    capacitance_file << "C_s = " << contacts[S]->c << std::endl;
    capacitance_file << "C_d = " << contacts[D]->c << std::endl;
    capacitance_file << "C_g = " << contacts[G]->c << std::endl;

    // produce plots of purely time-dependent observables and save them as PNGs
    if (plots) {
        gnuplot gp;
        gp << "set terminal png\n";
        gp << "set xlabel 't / ps'\n";

        // source and drain current
        gp << "set format x '%1.2f'\n";
        gp << "set format y '%1.0g'\n";
        arma::vec I_s(N_t);
        arma::vec I_d(N_t);
        for (unsigned i = 0; i < N_t; ++i) {
            I_s(i) = I[i].s();
            I_d(i) = I[i].d();
        }
        gp << "set title '" << name << " time-dependent source current'\n";
        gp << "set ylabel 'I_{s} / A'\n";
        gp << "set output '" << subfolder << "/I_s.png'\n";
        gp.add(std::make_pair(t * 1e12, I_s));
        gp.plot();
        gp.reset();
        gp << "set title '" << name << " time-dependent drain voltage'\n";
        gp << "set ylabel 'I_{d} / A'\n";
        gp << "set output '" << subfolder << "/I_d.png'\n";
        gp.add(std::make_pair(t * 1e12, I_d));
        gp.plot();
    }

    std::cout << " done!\n";
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
    q.col(S) = get_q(phi[0].s()).rows(3, mem + 3);
    q.col(D) = get_q(phi[0].d()).rows(3, mem + 3);

    // sum of two following q-values reversed
    qsum.resize(mem, 2);
    for (unsigned i = 0; i < mem; ++i) {
        qsum(i, S) = q(mem - i, S) + q(mem - 1 - i, S);
        qsum(i, D) = q(mem - i, D) + q(mem - 1 - i, D);
    }
    std::cout << "done!" << std::endl;
}

static inline std::vector<current> curve(const device_params & p, const std::vector<voltage<3>> & V) {
    // solves the steady state problem for a given set of voltages and returns the corresponding currents

    auto I = std::vector<current>(V.size());

//    #pragma omp parallel for
    for (unsigned i = 0; i < V.size(); ++i) {
        device d(p.name + std::to_string(i), p, V[i]);
        d.steady_state();
        I[i] = d.I[0];
        std::cout << "thread " << omp_get_thread_num() << ": voltage point " << i+1 << "/" << V.size() << " done" << std::endl;
    }
    return I;
}

//static inline arma::mat transfer(const device_params & p, double V_d, double V_g0, double V_g1, int N) {
//    // first column is V_g, second column I(V_g, V_d)

//    V_g = arma::linspace(V_g0, V_g1, N);
//    ret = arma::mat(N, 2);
//    std::copy(V_g.begin(), V_g.end(), I.colptr(0));

//    // prepare a vector of voltage points
//    std::vector<voltage<3>> voltage_points(N);
//    for (int j = 0; j < V_g.size(); ++j) {
//        voltage_points[i] = ({ 0, V_g(j), V_d });
//    }

//    // solve
//    std::vector<current> I = curve(p, voltage_points);

//    // fill return vector
//    for (int j = 0; j < V_g.size(); ++j) {
//        ret(j, 1) = I[j].total;
//    }

//    return ret;
//}


// ToDo: Overload??
template<bool csv>
static arma::mat transfer(const device_params & p, const arma::vec & V_d, double V_g0, double V_g1, int N) {
    // first column is V_g, second column I(V_g, V_d(1)), third column I(V_g, V_d(2)) etc

    auto V_g = arma::linspace(V_g0, V_g1, N);
    arma::mat ret(N, V_d.size() + 1);
    std::copy(V_g.begin(), V_g.end(), ret.colptr(0));

    // prepare a vector of voltage points
    std::vector<voltage<3>> voltage_points(V_d.size() * N);
    for (unsigned i = 0; i < V_d.size(); ++i) {
        for (int j = 0; j < N; ++j) {
            voltage_points[i*N + j] = { 0, V_g(j), V_d(i) };
        }
    }

    // solve
    std::vector<current> I = curve(p, voltage_points);

    // fill return vector
    for (unsigned i = 0; i < V_d.size(); ++i) {
        for (int j = 0; j < N; ++j) {
            ret(j, i+1) = I[i*N + j].total(0);
        }
    }

    if (csv) {
        std::ofstream of(save_folder() + "/" + p.name + "_transfer.csv");
        for (int j = 0; j < N; ++j) {
            for (unsigned i = 0; i < ret.n_cols; ++i) {
                of << ret(j, i) << "\t";
            }
            of << std::endl;
        }
        of.close();
    }

    return ret;
}





#endif

