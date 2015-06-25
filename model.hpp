#ifndef MODEL_HPP
#define MODEL_HPP

#include <string>
#include <array>

class model {
public:
    // parameters
    std::string name;
    double E_g;
    double m_eff;
    double E_gc;
    double m_efc;
    std::array<double, 3> F;
};

#endif

