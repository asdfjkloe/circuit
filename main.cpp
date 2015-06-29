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
#include "model.hpp"
#include "potential.hpp"
#include "voltage.hpp"
#include "wave_packet.hpp"

#undef CHARGE_DENSITY_HPP_BODY

#include "charge_density.hpp"

using namespace arma;
using namespace std;

int main() {
    circuit<3> c;
    int nfet_i = c.add_device(nfet);
    int pfet_i = c.add_device(pfet);
    c.link(nfet_i, S, S);
    c.link(nfet_i, D, pfet_i, D);
    c.link(nfet_i, G, G);
    c.link(pfet_i, S, D);
    c.link(pfet_i, G, G);

    c.contact(S)->V = 0.0;
    c.contact(G)->V = 0.2;
    c.contact(D)->V = 0.5;

    cout << "(nfet) V_s = " << c[nfet_i].contacts[S]->V << endl;
    cout << "(nfet) V_d = " << c[nfet_i].contacts[D]->V << endl;
    cout << "(nfet) V_g = " << c[nfet_i].contacts[G]->V << endl;
    cout << "(pfet) V_s = " << c[pfet_i].contacts[S]->V << endl;
    cout << "(pfet) V_d = " << c[pfet_i].contacts[D]->V << endl;
    cout << "(pfet) V_g = " << c[pfet_i].contacts[G]->V << endl;


    return 0;
}

