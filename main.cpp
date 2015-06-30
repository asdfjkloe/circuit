//#define ARMA_NO_DEBUG // no bound checks
//#define GNUPLOT_NOPLOTS

#include <armadillo>
#include <iostream>

#define CHARGE_DENSITY_HPP_BODY

#include "charge_density.hpp"
#include "circuit.hpp"
#include "constant.hpp"
#include "contact.hpp"
#include "current.hpp"
#include "device.hpp"
#include "device_params.hpp"
#include "geometry.hpp"
#include "green.hpp"
#include "inverter.hpp"
#include "model.hpp"
#include "potential.hpp"
#include "voltage.hpp"
#include "wave_packet.hpp"

#undef CHARGE_DENSITY_HPP_BODY

#include "charge_density.hpp"

using namespace arma;
using namespace std;

int main() {
    /*device p(pfet, {0.5, 0.5, 0.1});
    p.steady_state();
    plot_ldos(p.p, p.phi[0], 1000);
    plot(p.n[0].lv);
    plot(p.n[0].rv);
    plot(p.n[0].lc);
    plot(p.n[0].rc);
    std::cout << p.I[0].lv(0) << std::endl;
    std::cout << p.I[0].rv(0) << std::endl;
    std::cout << p.I[0].lc(0) << std::endl;
    std::cout << p.I[0].rc(0) << std::endl;
    std::cout << p.I[0].total(0) << std::endl;
    std::cout << p.phi[0].s() << std::endl;
    std::cout << p.phi[0].d() << std::endl;
    return 0;*/

    inverter inv(nfet, pfet, 1e-16);

    int N = 100;
    vec V_in = linspace(0.1, 0.4, N);
    vec V_out(N);

    for (int i = 0; i < N; ++i) {
        inv.steady_state({ 0.0, 0.5, V_in(i) });
        V_out(i) = inv.n().contacts[D]->V;
    }
    plot(V_out);

    return 0;
}

