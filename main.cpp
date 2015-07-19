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

int main(int argc, char ** argv) {
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

    omp_set_num_threads(18);
    cout << "saving results in " << save_folder(false, "TFET_RO") << endl;
    ring_oscillator<3> ro(ntfet, ptfet, 5e-18);
    ro.time_evolution(signal<2>(1e-10, voltage<2>{0.0, 0.1}));
    ro.save<true>();
    return 0;

//    omp_set_num_threads(stoi(argv[1]));

//    device d("prototype", ntfet);
//    double l_g = stod(argv[2]);
//    d.p.l_g = l_g;
//    d.p.update("ntfet_lg=" + string(argv[2]));

//    cout << "saving results in " << save_folder(false, "lg=20") << endl;
//    ofstream s(save_folder() + "geometry.ini");
//    s << d.p.to_string();
//    s.close();

//    transfer<true>(d.p, {{0, 0.1, -0.1}}, .1, 240);

//    return 0;

//    inverter inv(ntfet, ptfet, 5e-18);
//    int N = 300;
//    vec V_in = linspace(0, .1, N);
//    vec V_out(N);

//    for (int i = 0; i < N; ++i) {
//        inv.steady_state({0, 0.1, V_in(i)});
//        V_out(i) = inv.get_output(0)->V;
//    }

//    mat ret = join_horiz(V_in, V_out);
//    ret.save("tfet_inverter.csv", csv_ascii);

//    return 0;


}

