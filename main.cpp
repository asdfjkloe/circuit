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
    omp_set_nested(1);

    //flush denormal floats to zero for massive speedup
    //(i.e. set bits 15 and 6 in SSE control register MXCSR)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
}

int main() {
    setup();

//    cout << "saving results in " << save_folder(true, "ptfet_transfer") << endl;

//    transfer<true>(ptfet, {{0, -.1, -.2}, {0, -.2, -.2}}, .2, 200);
//    output<true>(ptfet, {{0, 0, .05}, {0, 0, .1}}, .3, 3);

//    device d("nfet", nfet, voltage<3>{0.0, 0.5, 0.3});
//    d.steady_state();
//    d.init_time_evolution(20);
//    wall_clock timer;
//    timer.tic();
//    for (int i = 1; i < 20; ++i) {
//        d.time_step();
//    }
//    cout << timer.toc() << endl;
//    return 0;

    ring_oscillator<3> ro(nfet, pfet, 5e-17);
    wall_clock timer;
    timer.tic();
    ro.time_evolution(signal<2>(1e-10, voltage<2>{0.0, 0.5}));
    cout << timer.toc() << endl;
//    ro.save<true>();

    return 0;
}

