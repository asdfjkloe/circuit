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
#include "ring_oscillator.hpp"
#include "voltage.hpp"
#include "wave_packet.hpp"

#undef CHARGE_DENSITY_HPP_BODY

#include "charge_density.hpp"

using namespace arma;
using namespace std;

int main() {
    ring_oscillator<9> ro(nfet, pfet, 1e-17);
    ro.steady_state({ 0.0, 0.5 });
    for (int i = 0; i < 9; ++i) {
        cout << ro.n(i).contacts[D]->V << endl;
    }

    return 0;

    inverter inv(nfet, pfet, 1e-16);

    int N = 20;
    vec V_in = linspace(0.0, 0.5, N);
    vec V_out(N);

    for (int i = 0; i < N; ++i) {
        inv.steady_state({ 0.0, 0.5, V_in(i) });
        V_out(i) = inv.n().contacts[D]->V;
    }
    plot(V_out);

    return 0;
}

