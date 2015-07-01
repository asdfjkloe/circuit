//#define ARMA_NO_DEBUG // no bound checks
//#define GNUPLOT_NOPLOTS

#include <armadillo>
#include <iostream>
#include <omp.h>
#include <xmmintrin.h>

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

static inline void setup_env() {
    // disable nested parallelism globally
    omp_set_nested(0);

    //flush denormal floats to zero for massive speedup
    //(i.e. set bits 15 and 6 in SSE control register MXCSR)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
}

int main() {
    setup_env();

    arma::vec V_d = {0.2, 0.3, 0.4, 0.5};
    transfer<true>(ntfet, V_d, 0, 0.4, 300);

//    cout << res << endl;

//    ring_oscillator<9> ro(nfet, pfet, 1e-17);
//    ro.steady_state({ 0.0, 0.5 });
//    for (int i = 0; i < 9; ++i) {
//        cout << ro.n(i).contacts[D]->V << endl;
//    }

//    return 0;

//    inverter inv(nfet, pfet, 1e-16);

//    int N = 20;
//    vec V_in = linspace(0.0, 0.5, N);
//    vec V_out(N);

//    for (int i = 0; i < N; ++i) {
//        inv.steady_state({ 0.0, 0.5, V_in(i) });
//        V_out(i) = inv.n().contacts[D]->V;
//    }
//    plot(V_out);

    return 0;
}

