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

template<ulint N>
static inline signal<N> operator+(const signal<N> & s1, const signal<N> s2);

template<ulint N>
static inline signal<N> step_signal(double T, const std::vector<double> & t, const std::vector<voltage<N>> & V);
template<ulint N>
static inline signal<N> linear_signal(double T, const voltage<N> & V0, const voltage<N> & V1);
template<ulint N>
static inline signal<N> linear_signal(double T, const std::vector<double> & t, const std::vector<voltage<N>> & V);
template<ulint N>
static inline signal<N> sine_signal(double T, const voltage<N> & V0, const voltage<N> & V_a, const double f, const double ph = 0);
template<ulint N>
static inline signal<N> cosine_signal(double T, const voltage<N> & V0, const voltage<N> & V_a, const double f, const double ph = 0);

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
/*
template<ulint N>
signal<N> operator+(const signal<N> & s1, const signal<N> & s2) {
    signal<N> s3(s1.T + s2.T);
    std::copy(s1.V.begin(), s1.V.end(), s3.V.begin());
    std::copy(s2.V.begin(), s2.V.end(), s3.V.begin() + s1.V.size());
}

template<ulint N>
signal<N> step_signal(double T, const std::vector<double> & t, const std::vector<voltage<N>> & V) {
    signal<N> s(T);

    std::vector<int> idx(t.size());
    for (ulint i = 0; i < t.size(); ++i) {
        idx[i] = std::round(t[i] / c::dt);
    }

    auto it0 = std::begin(s.V);
    auto it1 = std::begin(s.V);
    for (int i = 0; i < idx.size(); ++i) {
        it0 = it1;
        it1 = std::begin(s.V) + idx[i];
        std::fill(it0, it1, V[i]);
    }
    it0 = it1;
    it1 = std::end(s.V);
    std::fill(it0, it1, V[V.size() - 1]);
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
signal<N> linear_signal(double T, const std::vector<double> & t, const std::vector<voltage<N>> & V) {
    double t0;
    double t1 = 0.0;
    voltage<N> V0;
    voltage<N> V1 { 0 };
    for (ulint i = 0; i < t.size(); ++i) {
        t0 = t1;
        t1 = t[i];
        V0 = V1;
        V1 = V[i];
        linear_signal(t1 - t0, V0, V1);
    }
}*/

#endif

