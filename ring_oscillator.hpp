#ifndef RING_OSCILLATOR_HPP
#define RING_OSCILLATOR_HPP

#include "circuit.hpp"

template<int N>
class ring_oscillator : private circuit<2> {
public:
    inline ring_oscillator(const device_params & n, const device_params & p, double capacitance);

    inline const device & n(int i) const;
    inline device & n(int i);
    inline const device & p(int i) const;
    inline device & p(int i);

    inline bool steady_state(const std::array<double, 2> & V) override;
    using circuit<2>::time_step;

private:
    std::array<int, N> n_i;
    std::array<int, N> p_i;
};

//----------------------------------------------------------------------------------------------------------------------

template<int N>
ring_oscillator<N>::ring_oscillator(const device_params & n_, const device_params & p_, double capacitance) {
    // create devices
    for (int i = 0; i < N; ++i) {
        n_i[i] = add_device(n_);
        p_i[i] = add_device(p_);
    }

    // link devices
    for (int i = 0; i < N; ++i) {
        link(n_i[i], S, S);
        link(p_i[i], S, D);
        link(n_i[i], D, p_i[i], D);
    }
    for (int i = 0; i < N; ++i) {
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
bool ring_oscillator<N>::steady_state(const std::array<double, 2> & V) {
    int i;
    auto delta_I = [&] (double V_o) {
        n(i).contacts[D]->V = V_o;

        n(i).steady_state();
        p(i).steady_state();

        return n(i).I[0].d() + p(i).I[0].d();
    };

    contacts[S]->V = V[S];
    contacts[D]->V = V[D];

    // starting point
    n(0).contacts[G]->V = V[S];

    // solve each inverter, don't go back to the start
    for (i = 0; i < N; ++i) {
        double V_out;
        bool converged = brent(delta_I, V[S], V[D], 0.0005, V_out);
        std::cout << "i = " << i << "; V_out = " << V_out;
        std::cout << ", " << (converged ? "" : "ERROR!!!") << std::endl;
        if (!converged) {
            return false;
        }
    }

    return true;
}

#endif

