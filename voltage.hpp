#ifndef VOLTAGE_HPP
#define VOLTAGE_HPP

#include <array>

using voltage = std::array<double, 3>;

static inline voltage operator+(const voltage & a, double b) {
    return { a[0] + b, a[1] + b, a[2] + b };
}
static inline voltage operator+(double a, const voltage & b) {
    return b + a;
}
static inline voltage operator+(const voltage & a, const voltage & b) {
    return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}
static inline voltage operator-(const voltage & a) {
    return { -a[0], -a[1], -a[2] };
}
static inline voltage operator-(const voltage & a, double b) {
    return { a[0] - b, a[1] - b, a[2] - b };
}
static inline voltage operator-(double a, const voltage & b) {
    return a + (-b);
}
static inline voltage operator-(const voltage & a, const voltage & b) {
    return a + (-b);
}
static inline voltage operator*(const voltage & a, double b) {
    return { a[0] * b, a[1] * b, a[2] * b };
}
static inline voltage operator*(double a, const voltage & b) {
    return b * a;
}
static inline voltage operator*(const voltage & a, const voltage & b) {
    return { a[0] * b[0], a[1] * b[1], a[2] * b[2] };
}
static inline voltage operator/(const voltage & a, double b) {
    return a * (1.0 / b);
}
static inline voltage operator/(double a, const voltage & b) {
    return { a / b[0], a / b[1], a / b[2] };
}
static inline voltage operator/(const voltage & a, const voltage & b) {
    return { a[0] / b[0], a[1] / b[1], a[2] / b[2] };
}
template<class F>
static inline voltage func(F && f, const voltage & a) {
    voltage ret;

    ret[0] = f(a[0]);
    ret[1] = f(a[1]);
    ret[2] = f(a[2]);

    return ret;
}

#endif

