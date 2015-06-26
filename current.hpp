#ifndef CURRENT_HPP
#define CURRENT_HPP

#include <armadillo>

#include "geometry.hpp"
#include "model.hpp"
#include "potential.hpp"

class current {
public:
    arma::vec lv;    // current from left valence band
    arma::vec rv;    // current from right valence band
    arma::vec lc;    // current from left conduction band
    arma::vec rc;    // current from right conduction band
    arma::vec total; // total current

    inline current();
    inline current(const geometry & g, const model & m, const potential & phi);
    //inline current(const device & d, const wave_packet psi[4], const potential & phi);

    inline double s() const;
    inline double d() const;
};

//----------------------------------------------------------------------------------------------------------------------

current::current() {
}

current::current(const geometry & g, const model & m, const potential & phi) {
}

double current::s() const {
    return total[0];
}

double current::d() const {
    return total[total.size() - 1];
}

#endif

