//#define ARMA_NO_DEBUG // no bound checks
//#define GNUPLOT_NOPLOTS

#include <armadillo>
#include <iostream>

#define CHARGE_DENSITY_HPP_BODY

#include "charge_density.hpp"
#include "constant.hpp"
#include "contact.hpp"
#include "current.hpp"
#include "device.hpp"
#include "geometry.hpp"
#include "green.hpp"
#include "model.hpp"
#include "potential.hpp"
#include "voltage.hpp"

#undef CHARGE_DENSITY_HPP_BODY

#include "charge_density.hpp"

using namespace arma;
using namespace std;

int main() {
    device d(fet_geometry, nfet_model);
    d.steady_state();
    return 0;
}

