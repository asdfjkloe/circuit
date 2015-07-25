#ifndef SIGNAL_HPP
#define SIGNAL_HPP

#include <vector>

#include "constant.hpp"
#include "voltage.hpp"

template<ulint N>
class signal {
public:
    int N_t;
    double T;
    std::vector<voltage<N>> V;

    inline signal();
    inline signal(double T);
    inline signal(double T, const voltage<N> & V);

    inline voltage<N> & operator[](int index);
    inline const voltage<N> & operator[](int index) const;
};

// append two signals
template<ulint N>
static inline signal<N> operator+(const signal<N> & s1, const signal<N> & s2);

// create a linear signal: V[i] = V0 + (V1 - V0) * t_i / T
template<ulint N>
static inline signal<N> linear_signal(double T, const voltage<N> & V0, const voltage<N> & V1);

// create a sine signal: V[i] = V0 + VA * sin(f * t_i  + ph)
template<ulint N>
static inline signal<N> sine_signal(double T, const voltage<N> & V0, const voltage<N> & VA, const double f, const double ph = 0);

// create a cosine signal: V[i] = V0 + VA * cos(f * t_i + ph)
template<ulint N>
static inline signal<N> cosine_signal(double T, const voltage<N> & V0, const voltage<N> & VA, const double f, const double ph = 0);

//----------------------------------------------------------------------------------------------------------------------

template<ulint N>
signal<N>::signal() {
}
template<ulint N>
signal<N>::signal(double T_)
    :  N_t(std::round(T_ / c::dt)), T(N_t * c::dt), V(N_t) {
}
template<ulint N>
signal<N>::signal(double T_, const voltage<N> & V_)
    : signal(T_) {
    std::fill(V.begin(), V.end(), V_);
}

template<ulint N>
voltage<N> & signal<N>::operator[](int index) {
    return V[index];
}
template<ulint N>
const voltage<N> & signal<N>::operator[](int index) const {
    return V[index];
}

template<ulint N>
signal<N> operator+(const signal<N> & s1, const signal<N> & s2) {
    signal<N> s3(s1.T + s2.T);
    std::copy(s1.V.begin(), s1.V.end(), s3.V.begin());
    std::copy(s2.V.begin(), s2.V.end(), s3.V.begin() + s1.V.size());
}

template<ulint N>
signal<N> linear_signal(double T, const voltage<N> & V0, const voltage<N> & V1) {
    signal<N> s(T);

    for (int i = 0; i < s.N_t; ++i) {
        double r = ((double)i) / ((double)(s.N_t - 1));
        s[i] = V0 * (1.0 - r) + V1 * r;
    }

    return s;
}
template<ulint N>
signal<N> sine_signal(double T, const voltage<N> & V0, const voltage<N> & VA, const double f, const double ph) {
    signal<N> s(T);

    for (int i = 0; i < s.N_t; ++i) {
        double t = i * s.dt;
        s[i] = V0 + VA * std::sin(t * f + ph);
    }

    return s;
}
template<ulint N>
signal<N> cosine_signal(double T, const voltage<N> & V0, const voltage<N> & VA, const double f, const double ph) {
    return sine_signal(T, V0, VA, f, ph + 0.5 * M_PI);
}

#endif

