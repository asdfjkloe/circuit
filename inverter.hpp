#ifndef INVERTER_HPP
#define INVERTER_HPP

#include "circuit.hpp"
#include "voltage.hpp"
#include "util/brent.hpp"

class inverter : public circuit<3, 1> {
public:
    inline inverter(const device_params & n, const device_params & p, double capacitance);

    inline const device & n() const;
    inline device & n();
    inline const device & p() const;
    inline device & p();

    inline bool steady_state(const voltage<3> & V) override;
    using circuit<3, 1>::time_step;
    using circuit<3, 1>::time_evolution;
    using circuit<3, 1>::save;

private:
    int n_i;
    int p_i;
};

//----------------------------------------------------------------------------------------------------------------------

inverter::inverter(const device_params & n, const device_params & p, double capacitance) {
    n_i = add_device(n.name, n);
    p_i = add_device(p.name, p);

    // link devices
    link_input(n_i, S, S);  // to ground
    link_input(n_i, G, G);  // to V_in
    link_output(n_i, D, 0); // to V_out
    link_input(p_i, S, D);  // to V_dd
    link_input(p_i, G, G);  // to V_in
    link_output(p_i, D, 0); // to V_out

    // set capacitance
    outputs[0]->c = capacitance;
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

bool inverter::steady_state(const voltage<3> & V) {
    auto delta_I = [&] (double V_o) {
        outputs[0]->V = V_o;

    n().steady_state();
    p().steady_state();

        return n().I[0].d() + p().I[0].d();
    };

    // set input voltages
    inputs[S]->V = V[S];
    inputs[D]->V = V[D];
    inputs[G]->V = V[G];

    V_out.resize(1);
    bool converged = brent(delta_I, V[S], V[D], 0.0001, V_out[0][0]);
    std::cout << "\nV_in = " << V[G] << " -> V_out = " << V_out[0][0];
    std::cout << (converged ? "" : " ERROR!!!") << std::endl << std::endl;

    return converged;
}

#endif

