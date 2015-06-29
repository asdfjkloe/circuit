#ifndef INVERTER_HPP
#define INVERTER_HPP

#include "circuit.hpp"
#include "voltage.hpp"
#include "util/brent.hpp"

class inverter : private circuit<3> {
public:
    inline inverter(const device_params & n, const device_params & p, double capacitance);

    inline const device & n() const;
    inline device & n();
    inline const device & p() const;
    inline device & p();

    inline bool steady_state(const voltage & V) override;

private:
    int n_i;
    int p_i;
};

//----------------------------------------------------------------------------------------------------------------------

inverter::inverter(const device_params & n, const device_params & p, double capacitance) {
    n_i = add_device(n);
    p_i = add_device(p);

    // link devices
    link(n_i, S, S);
    link(n_i, G, G);
    link(n_i, D, p_i, D);
    link(p_i, S, D);
    link(p_i, G, G);

    devices[n_i].contacts[D]->c = capacitance;
}

const device & inverter::n() const {
    return devices[n_i];
}
device & inverter::n() {
    return devices[n_i];
}
const device & inverter::p() const {
    return devices[p_i];
}
device & inverter::p() {
    return devices[p_i];
}

bool inverter::steady_state(const voltage & V) {
    auto delta_I = [&] (double V_o) {
        n().contacts[D]->V = V_o;

        n().steady_state();
        p().steady_state();

        return n().I[0].d() + p().I[0].d();
    };

    n().contacts[S]->V = V[S];
    p().contacts[S]->V = V[D];
    n().contacts[G]->V = V[G];

    double V_out;
    bool converged = brent(delta_I, V[S], V[D], 0.0005, V_out);
    std::cout << "V_out = " << V_out << std::endl;

    return converged;
}

#endif

