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

// set the type of device
static device_params ntype = ntfet;
static device_params ptype = ptfet;

static inline void setup() {
    // disable nested parallelism globally
    omp_set_nested(0);

    //flush denormal floats to zero for massive speedup
    //(i.e. set bits 15 and 6 in SSE control register MXCSR)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
}


static inline void trans (int argc, char ** argv) {
    // starts transfer-curve simulations with a certain gate-length

    double l_g = stod(argv[3]);
    double Vg0 = stod(argv[4]);
    double Vg1 = stod(argv[5]);
    double Vd  = stod(argv[6]);
    int N      = stod(argv[7]);

    stringstream ss;
    ss << "transfer/lg=" << l_g << "_Vd=" << Vd;
    cout << "saving results in " << save_folder(true, ss.str()) << endl;

    device d("prototype", ntype);
    d.p.l_g = l_g;
    d.p.update(ss.str());

    transfer<true>(d.p, {{0, Vd, Vg0}}, Vg1, N);

    ofstream s(save_folder() + "/parameters.ini");
    s << d.p.to_string();
    s.close();
}

static inline void outp (int argc, char ** argv) {
    // starts output-curve simulations with a certain gate-length

    double l_g = stod(argv[3]);
    double Vd0 = stod(argv[4]);
    double Vd1 = stod(argv[5]);
    double Vg  = stod(argv[6]);
    int N      = stod(argv[7]);

    stringstream ss;
    ss << "output/lg=" << l_g;
    cout << "saving results in " << save_folder(true, ss.str()) << endl;

    device d("prototype", ntype);
    d.p.l_g = l_g;
    d.p.update(ss.str());

    output<true>(d.p, {{0, Vd0, Vg}}, Vd1, N);

    ofstream s(save_folder() + "/parameters.ini");
    s << d.p.to_string();
    s.close();
}

static void inv (int argc, char ** argv) {
    // starts a static inverter simulation

    double Vin0 = stod(argv[3]);
    double Vin1 = stod(argv[4]);
    double V_dd = stod(argv[5]);
    int    N    = stoi(argv[6]);

    inverter inv(ntype, ptype);
    vec V_in = linspace(Vin0, Vin1, N);
    vec V_out(N);

    for (int i = 0; i < N; ++i) {
        inv.steady_state({0, V_dd, V_in(i)});
        V_out(i) = inv.get_output(0)->V;
    }

    stringstream ss;
    ss << "ntd_inverter/V=" << V_dd;
    cout << "saving results in " << save_folder(true, ss.str()) << endl;

    mat ret = join_horiz(V_in, V_out);
    ret.save(save_folder() + "/inverter_curve.csv", csv_ascii);

    ofstream sn(save_folder() + "/parameters_ntype.ini");
    sn << ntype.to_string();
    sn.close();

    ofstream sp(save_folder() + "/parameters_ptype.ini");
    sp << ptype.to_string();
    sp.close();
}

static inline void ro (int argc, char ** argv) {
    // starts a transient ring-oscillator simulation
    double T = stod(argv[3]);
    double C = stod(argv[4]);
    double V_dd = stod(argv[5]);

    stringstream ss;
    ss << "ring_oscillator/T=" << T << "_C=" << C << "_V=" << V_dd;
    cout << "saving results in " << save_folder(true, ss.str()) << endl;

    ring_oscillator<3> ro(ntype, ptype, C);
    ro.time_evolution(signal<2>(T, voltage<2>{0.0, V_dd}));
    ro.save<true>();

    ofstream sn(save_folder() + "/parameters_ntype.ini");
    sn << ntype.to_string();
    sn.close();

    ofstream sp(save_folder() + "/parameters_ptype.ini");
    sp << ptype.to_string();
    sp.close();
}

int main(int argc, char ** argv) {
    setup();

    // first argument is always the number of threads
    // (can not be more than the number specified when compiling openBLAS)
    omp_set_num_threads(stoi(argv[1]));

    // second argument chooses the type of simulation
    string stype(argv[2]);
    if (stype == "trans") {
        trans(argc, argv);
    } else if (stype == "outp") {
        outp(argc, argv);
    } else if (stype == "inv") {
        inv(argc, argv);
    } else if (stype == "ro") {
        ro(argc, argv);
    }

    return 0;
}


