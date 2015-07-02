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

static inline void setup() {
    // disable nested parallelism globally
    omp_set_nested(0);

    //flush denormal floats to zero for massive speedup
    //(i.e. set bits 15 and 6 in SSE control register MXCSR)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
}

int main() {
    setup();

//    cout << "saving results in " << save_folder(true, "nfet") << endl;

//    transfer<true>(nfet, {{0.0, +0.2, 0.0}, {0.0, +0.4, 0.0}}, +0.2, 100);
//    transfer<true>(pfet, {{0.0, -0.2, 0.0}, {0.0, -0.4, 0.0}}, -0.2, 100);
    /*
    transfer<true>(ntfet, vector<voltage<3>>{voltage<3>{0, -.2, .1}, voltage<3>{0, -.2, .2}}, .2, 200);

    output<true>(ntfet, vector<voltage<3>>{voltage<3>{0, .05, 0}, voltage<3>{0, .1, 0}}, .3, 200);
*/

//    device d1("test1", ntfet, voltage<3>{0, 0.2, 0.});
//    d1.steady_state();
//    plot_ldos(d1.p, d1.phi[0]);

//    device d2("test2", ntfet, voltage<3>{0, 0.2, 0.});
//    d2.steady_state();
//    plot_ldos(d2.p, d2.phi[0]);

    ring_oscillator<3> ro(nfet, pfet, 1e-19);
    ro.time_evolution(signal<2>(5e-11, voltage<2>{0.0, 0.5}));
    ro.save<true>();
    return 0;
}

