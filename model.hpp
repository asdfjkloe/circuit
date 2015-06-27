#ifndef MODEL_HPP
#define MODEL_HPP

#include <armadillo>
#include <array>

#include "constant.hpp"

class model {
public:
    // parameters
    double E_g;
    double m_eff;
    double E_gc;
    double m_efc;
    std::array<double, 3> F;
};

static const model nfet_model {
    0.62,             // E_g
    0.01 * c::m_e,    // m_eff
    0.62,             // E_gc
    0.01 * c::m_e,    // m_efc
    {
        0.62 / 2 + 0.015, // F[S]
        0.62 / 2 + 0.015, // F[D]
        0.00              // F[G]
    }
};

static const model pfet_model {
    nfet_model.E_g,   // E_g
    nfet_model.m_eff, // m_eff
    nfet_model.E_gc,  // E_gc
    nfet_model.m_efc, // m_efc
    {
        -nfet_model.F[S],  // F[S]
        -nfet_model.F[D],  // F[D]
        -nfet_model.F[G]   // F[G]
    }
};

static const model ntfet_model {
    0.62,             // E_g
    0.05 * c::m_e,    // m_eff
    0.62,             // E_gc
    0.05 * c::m_e,    // m_efc
    {
        -0.62 / 2 - 0.015, // F[S] (p++)
        +0.62 / 2 + 0.001, // F[D] (n+)
         0.00              // F[G]
    }
};

static const model ptfet_model {
    ntfet_model.E_g,   // E_g
    ntfet_model.m_eff, // m_eff
    ntfet_model.E_gc,  // E_gc
    ntfet_model.m_efc, // m_efc
    {
        -ntfet_model.F[S],  // F[S] (n++)
        -ntfet_model.F[D],  // F[D] (p+)
        -ntfet_model.F[G]   // F[G]
    }
};

#endif

