#ifndef MODEL_HPP
#define MODEL_HPP

#include <armadillo>
#include <array>
#include <string>

#include "constant.hpp"

class model {
public:
    // parameters
    double E_g;
    double m_eff;
    double E_gc;
    double m_efc;
    std::array<double, 3> F;
    double shift;

    inline std::string to_string();
};

//static const model nfet_model {
//    0.62,             // E_g
//    0.02 * c::m_e,    // m_eff
//    0.62,             // E_gc
//    0.02 * c::m_e,    // m_efc
//    {
//        0.62 / 2 + 0.015, // F[S]
//        0.62 / 2 + 0.001, // F[D]
//        0.2               // F[G]
//    }
//};

//static const model nfetc_model {
//    0.62,             // E_g
//    0.01 * c::m_e,    // m_eff
//    0.4,              // E_gc
//    0.05 * c::m_e,    // m_efc
//    {
//        0.62 / 2 + 0.015, // F[S]
//        0.62 / 2 + 0.015, // F[D]
//        0.2               // F[G]
//    }
//};

//static const model pfet_model {
//    nfet_model.E_g,   // E_g
//    nfet_model.m_eff, // m_eff
//    nfet_model.E_gc,  // E_gc
//    nfet_model.m_efc, // m_efc
//    {
//        -nfet_model.F[S],  // F[S]
//        -nfet_model.F[D],  // F[D]
//        -nfet_model.F[G]   // F[G]
//    }
//};

static const model ntfetc_model { // most general model
    0.552,            // E_g (16, 0)
    0.04 * c::m_e,    // m_eff (
    0.30,             // E_gc
    0.10 * c::m_e,    // m_efc
    {
        -0.5 / 2 - 0.050, // F[S] (p++)
        +0.5 / 2 + 0.001, // F[D] (n+)
         0                // F[G]
    },
    0.00              // gate volatge shift
};

static const model ntfet_model {
    ntfetc_model.E_g,   // E_g
    ntfetc_model.m_eff, // m_eff
    ntfetc_model.E_g,   // same as E_g
    ntfetc_model.m_eff, // same as m_eff
    {
        ntfetc_model.F[S],
        ntfetc_model.F[D],
        ntfetc_model.F[G],
    },
    ntfetc_model.shift
};

static const model ptfet_model {
    ntfet_model.E_g,   // E_g
    ntfet_model.m_eff, // m_eff
    ntfet_model.E_gc,  // E_gc
    ntfet_model.m_efc, // m_efc
    { // reversed doping:
        -ntfet_model.F[S],  // F[S] (n++)
        -ntfet_model.F[D],  // F[D] (p+)
        -ntfet_model.F[G],  // F[G]
    },
    -ntfet_model.shift // characteristics are perfectly symmetrical
};

std::string model::to_string() {
    using namespace std;

    stringstream ss;

    ss << "E_g     = " << E_g   << endl;
    ss << "m_eff   = " << m_eff << endl;
    ss << "E_gc    = " << E_gc  << endl;
    ss << "m_efc   = " << m_efc << endl;
    ss << "F_s     = " << F[S]  << endl;
    ss << "F_d     = " << F[D]  << endl;
    ss << "F_g     = " << F[G]  << endl;
    ss << "shift   = " << shift << endl;

    return ss.str();
}

#endif

