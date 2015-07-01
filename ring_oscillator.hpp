#ifndef RING_OSCILLATOR_HPP
#define RING_OSCILLATOR_HPP

#include "circuit.hpp"

template<int N>
class ring_oscillator : private circuit<2, 1> {
public:
    inline ring_oscillator(const device_params & n, const device_params & p, double capacitance);

    inline const device & n(int i) const;
    inline device & n(int i);
    inline const device & p(int i) const;
    inline device & p(int i);

    inline bool steady_state(const voltage<2> & V) override;
    using circuit<2, 1>::time_step;

    template<bool plots = false>
    inline void save();
private:
    std::array<int, N> n_i;
    std::array<int, N> p_i;
};

//----------------------------------------------------------------------------------------------------------------------

template<int N>
ring_oscillator<N>::ring_oscillator(const device_params & n_, const device_params & p_, double capacitance) {
    // create devices
    for (int i = 0; i < N; ++i) {
        // give each device a distinct name
        std::stringstream ss;
        ss << "_" << i;
        std::string suffix = ss.str();

        // add devices and save indices
        n_i[i] = add_device(n_.name + suffix, n_);
        p_i[i] = add_device(p_.name + suffix, p_);
    }

    // link devices
    for (int i = 0; i < N - 1; ++i) {
        link_input(n_i[i], S, S);   // to ground
        link_input(p_i[i], S, D);   // to V_dd
        link(n_i[i], D, p_i[i], D); // common output port
    }
    link_output(n_i[N - 1], D, 0); // to circuit output
    link_output(p_i[N - 1], D, 0); // to circuit output
    for (int i = 0; i < N; ++i) {
        // link previous output to current input
        // last output is fed back to first input
        link(n_i[i], G, n_i[(i + N - 1) % N], D);
        link(p_i[i], G, p_i[(i + N - 1) % N], D);
    }

    // set capacitance
    for (int i = 0; i < N; ++i) {
        n(i).contacts[G]->c = capacitance;
    }
}

template<int N>
const device & ring_oscillator<N>::n(int i) const {
    return devices[n_i[i]];
}
template<int N>
device & ring_oscillator<N>::n(int i) {
    return devices[n_i[i]];
}
template<int N>
const device & ring_oscillator<N>::p(int i) const {
    return devices[p_i[i]];
}
template<int N>
device & ring_oscillator<N>::p(int i) {
    return devices[p_i[i]];
}

template<int N>
bool ring_oscillator<N>::steady_state(const voltage<2> & V) {
    // NOTE: there is no steady state for a ring oscillator,
    //       but a self-consistent solution is needed as a starting point for time-evolutions!
    int i;
    auto delta_I = [&] (double V_o) {
        n(i).contacts[D]->V = V_o;

        n(i).steady_state();
        p(i).steady_state();

        return n(i).I[0].d() + p(i).I[0].d();
    };

    // set input voltages
    inputs[S]->V = V[S];
    inputs[D]->V = V[D];

    // starting point
    outputs[0]->V = V[S];

    // solve each inverter seperately, don't go back to the start
    for (i = 0; i < N; ++i) {
        double V_o;
        bool converged = brent(delta_I, V[S], V[D], 0.0005, V_o);
        std::cout << "i = " << i << "; V_out = " << V_o;
        std::cout << (converged ? "" : ", ERROR!!!") << std::endl;
        if (!converged) {
            return false;
        }
    }

    // save output voltage
    V_out.resize(1);
    V_out[0][0] = outputs[0]->V;

    return true;
}

#endif

