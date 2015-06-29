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
    circuit c;
    int nfet_i = c.add_device(nfet);
    int pfet_i = c.add_device(pfet);
    c.link(nfet_i, G, pfet_i, G);
    c.link(nfet_i, D, pfet_i, D);

    c[nfet_i].contacts[S]->V = 0.0;
    c[nfet_i].contacts[D]->c = 1e-16;
    c[nfet_i].contacts[G]->V = 0.2;
    c[pfet_i].contacts[S]->V = 0.5;

    cout << c[nfet_i].contacts[S]-> V << endl;
    cout << c[nfet_i].contacts[D]-> V << endl;
    cout << c[nfet_i].contacts[G]-> V << endl;
    cout << c[pfet_i].contacts[S]-> V << endl;
    cout << c[pfet_i].contacts[D]-> V << endl;
    cout << c[pfet_i].contacts[G]-> V << endl;

    return 0;
}

