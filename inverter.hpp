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
    using circuit<3>::time_step;

    template<bool plots = false>
    inline void save();

private:
    int n_i;
    int p_i;
};

//----------------------------------------------------------------------------------------------------------------------

inverter::inverter(const device_params & n, const device_params & p, double capacitance) {
    n_i = add_device(n);
    p_i = add_device(p);

    // link devices
    link(n_i, S, S); // to ground
    link(n_i, G, G); // to input
    link(n_i, D, p_i, D); // common output terminal
    link(p_i, S, D); // to V_dd
    link(p_i, G, G); // to input

    // set capacitance
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

    contacts[S]->V = V[S];
    contacts[D]->V = V[D];
    contacts[G]->V = V[G];

    double V_out;
    bool converged = brent(delta_I, V[S], V[D], 0.0005, V_out);
    std::cout << "V_out = " << V_out;
    std::cout << ", " << (converged ? "" : "ERROR!!!") << std::endl;

    return converged;
}

template<bool plots>
void inverter::save() {
    n().save<plots>();
    p().save<plots>();

    V_out.save(save_folder() + "/V_out.arma");
    std::ofstream just_C(save_folder() + "/C.txt");
    just_C << capacitance;
    just_C.close();

    if (plots) {
        // make a plot of V_out and save it as a PNG
        gnuplot gp;
        gp << "set terminal png\n";
        gp << "set title 'Inverter output voltage'\n";
        gp << "set xlabel 't / ps'\n";
        gp << "set ylabel 'V_{out} / V'\n";
        gp << "set format x '%1.2f'\n";
        gp << "set format y '%1.2f'\n";
        gp << "set output '" << save_folder() << "/V_out.png'\n";
        gp.add(std::make_pair(sg.t * 1e12, V_out));
        gp.plot();
    }
}

#endif

