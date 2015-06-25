#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#include <armadillo>
#include <map>

#include "charge_density.hpp"
#include "constant.hpp"
#include "geometry.hpp"
#include "model.hpp"
#include "voltage.hpp"
#include "util/anderson.hpp"

class potential {
public:
    arma::vec data;
    arma::vec twice;

    inline potential();
    inline potential(const geometry & g, const arma::vec & R);
    inline potential(const geometry & g, const arma::vec & R0, const charge_density & n);
    inline potential(const geometry & g, const model & m, const voltage & V);
    inline double update(const geometry & g, const arma::vec & R0, const charge_density & n, anderson & mr_neo);

    inline double & operator()(int index);
    inline const double & operator()(int index) const;
    inline double s() const;
    inline double d() const;

    inline void smooth(const geometry & g, const model & m);

    static inline arma::vec get_R0(const geometry & g, const model & m, const voltage & V);

private:
    enum {
        L = 0, // left
        R = 1, // right
        I = 2, // inside
        O = 3  // outside
    };

    inline void update_twice();

    template<bool minmax>
    inline void smooth(unsigned x0, unsigned x1);

    static inline arma::vec get_R(const geometry & g, const arma::vec & R0, const charge_density & n);
    static inline const arma::sp_mat & get_S(const geometry & g);
    static inline const std::array<arma::mat, 4> & get_eps(const geometry & g);
};

//----------------------------------------------------------------------------------------------------------------------

potential::potential() {
}
potential::potential(const geometry & g, const arma::vec & R) {
    using namespace arma;
    vec phi2D = spsolve(get_S(g), R);
    data = phi2D({uword((g.M_cnt - 1) * g.N_x), uword(g.M_cnt * g.N_x - 1)});

    update_twice();
}
potential::potential(const geometry & g, const arma::vec & R0, const charge_density & n)
    : potential(g, get_R(g, R0, n)) {
}
potential::potential(const geometry & g, const model & m, const voltage & V)
    : potential(g, get_R0(g, m, V)) {
}

double potential::update(const geometry & g, const arma::vec & R0, const charge_density & n, anderson & mr_neo) {
    // calculate undamped update
    potential phi_update(g, R0, n);

    // anderson mixing
    arma::vec f = phi_update.data - data;
    mr_neo.update(data, f);

    update_twice();

    // return dphi
    return max(abs(f));
}

double & potential::operator()(int index) {
    return data[index];
}
const double & potential::operator()(int index) const {
    return data[index];
}
double potential::s() const {
    return data[0];
}
double potential::d() const {
    return data[data.size() - 1];
}

void potential::smooth(const geometry & g, const model & m) {
    // smooth source region
    if (m.F[S] > 0) {
        smooth<true>(0, g.N_sc + g.N_sox + g.N_sg + g.N_g * 0.05);
    } else {
        smooth<false>(0, g.N_sc + g.N_sox + g.N_sg + g.N_g * 0.05);
    }

    // smooth drain region
    if (m.F[D] > 0) {
        smooth<true>(g.N_sc + g.N_sox + g.N_sg + g.N_g * 0.95, g.N_x);
    } else {
        smooth<false>(g.N_sc + g.N_sox + g.N_sg + g.N_g * 0.95, g.N_x);
    }

    update_twice();
}

arma::vec potential::get_R0(const geometry & g, const model & m, const voltage & V0) {
    using namespace arma;

    // shortcuts
    double dr2 = 1.0 / g.dr / g.dr;
    double dx2 = 1.0 / g.dx / g.dx;

    // get reference to 4 eps matrices (one for each direction)
    const std::array<arma::mat, 4> & eps = get_eps(g);

    // add built-in voltages
    std::array<double, 3> V = - (V0 + m.F);

    // init return vector
    vec R0(g.N_x * g.M_r);
    R0.fill(0.0);

    // create R0 from geometry and voltages
    int i, j, k;
    double r, rp;

    // add dirichlet boundary conditions where necessary
    j = g.M_cnt - 1;
    r = j * g.dr + 0.5 * g.dr;
    rp = r + 0.5 * g.dr;
    k = j * g.N_x;
    for (i = 0; i < g.N_sc; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[S];
    }
    k = j * g.N_x + g.N_x - g.N_dc;
    for (i = g.N_x - g.N_dc; i < g.N_x; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[D];
    }
    for (j = g.M_cnt; j < g.M_cnt + g.M_ox - 1; ++j) {
        r = j * g.dr + 0.5 * g.dr;
        R0(k) -= dx2 * r * eps[L](g.N_sc, j) * V[S];
        k += g.N_sox + g.N_sg + g.N_g + g.N_dg + g.N_dox - 1;
        R0(k++) -= dx2 * r * eps[R](g.N_x - g.N_dc - 1, j) * V[D];
    }
    r = j * g.dr + 0.5 * g.dr;
    rp = r + 0.5 * g.dr;
    R0(k) -= dx2 * r * eps[L](g.N_sc, j) * V[S];
    for (i = g.N_sc; i < g.N_sc + g.N_sox; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[S];
    }
    k += g.N_sg;
    for (i = g.N_sc + g.N_sox + g.N_sg; i < g.N_sc + g.N_sox + g.N_sg + g.N_g; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[G];
    }
    k += g.N_dg;
    for (i = g.N_x - g.N_dc - g.N_dox; i < g.N_x - g.N_dc; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[D];
    }
    R0(k - 1) -= dx2 * r * eps[R](g.N_x - g.N_dc - 1, j) * V[D];

    for (j = g.M_cnt + g.M_ox; j < g.M_r; ++j) {
        r = j * g.dr + 0.5 * g.dr;
        R0(k) -= dx2 * r * eps[L](g.N_sc + g.N_sox, j) * V[S];
        k += g.N_sg - 1;
        R0(k++) -= dx2 * r * eps[R](g.N_sc + g.N_sox + g.N_sg - 1, j) * V[G];
        R0(k) -= dx2 * r * eps[L](g.N_sc + g.N_sox + g.N_sg + g.N_g, j) * V[G];
        k += g.N_dg - 1;
        R0(k++) -= dx2 * r * eps[R](g.N_x - g.N_dc - g.N_dox - 1, j) * V[D];
    }

    // shrink to fit and return
    R0.resize(k);
    return R0;
}

void potential::update_twice() {
    // duplicate each entry
    for (unsigned i = 0; i < data.size(); ++i) {
        twice(2 * i    ) = data(i);
        twice(2 * i + 1) = data(i);
    }
}

template<bool minmax>
void potential::smooth(unsigned x0, unsigned x1) {
    using namespace arma;

    if (minmax) {
        for (unsigned i = x0; i < x1 - 1; ++i) {
            if (data(i+1) >= data(i)) {
                continue;
            }
            for (unsigned j = i + 1; j < x1; ++j) {
                if (data(j) >= data(i)) {
                    data({i+1, j-1}).fill(data(i));
                    break;
                }
            }
        }
        for (unsigned i = x1 - 1; i >= x0 + 1; --i) {
            if (data(i-1) >= data(i)) {
                continue;
            }
            for (unsigned j = i - 1; j >= 1; --j) {
                if (data(j) >= data(i)) {
                    data({j+1, i-1}).fill(data(i));
                    break;
                }
            }
        }
    } else {
        for (unsigned i = x0; i < x1 - 1; ++i) {
            if (data(i+1) <= data(i)) {
                continue;
            }
            for (unsigned j = i + 1; j < x1; ++j) {
                if (data(j) <= data(i)) {
                    data({i+1, j-1}).fill(data(i));
                    break;
                }
            }
        }
        for (unsigned i = x1 - 1; i >= x0 + 1; --i) {
            if (data(i-1) <= data(i)) {
                continue;
            }
            for (unsigned j = i - 1; j >= 1; --j) {
                if (data(j) <= data(i)) {
                    data({j+1, i - 1}).fill(data(i));
                    break;
                }
            }
        }
    }
}

arma::vec potential::get_R(const geometry & g, const arma::vec & R0, const charge_density & n) {
    using namespace arma;
    vec R = R0;

    // add non-constant part to R0 due to charge_density at the boundary of the cnt
    R({uword((g.M_cnt - 1) * g.N_x), uword(g.M_cnt * g.N_x - 1)}) += n.total * g.r_cnt * 1e9; // 10^9 because of m->nm in eps_0
    return R;
}
const arma::sp_mat & potential::get_S(const geometry & g) {
    // check if S was already calculated for this geometry
    static std::map<std::string, arma::sp_mat> S;
    auto it = S.find(g.name);
    if (it != std::end(S)) {
        return it->second;
    }

    // calculate new S
    using namespace arma;

    // shortcuts
    double dr2 = 1.0 / g.dr / g.dr;
    double dx2 = 1.0 / g.dx / g.dx;

    // get reference to 4 eps matrices (one for each direction)
    const std::array<arma::mat, 4> & eps = get_eps(g);

    // helper objects to create sparse matrix faster
    umat indices(2, g.M_r * g.N_x * 5);
    vec values(     g.M_r * g.N_x * 5);
    uword N_v = 0;
    auto set_value = [&indices, &values, &N_v] (uword i, uword j, double val) {
        indices(0, N_v) = i;
        indices(1, N_v) = j;
        values(N_v++) = val;
    };

    // start values
    int k     = 0;     // current main diagonal element
    int i0    = 0;     // left i limit
    int i1    = g.N_x; // right i limit
    int delta = g.N_x; // distance to next vertical coupling off diagonal

    // main loop over rows
    for (int j = 0; j < g.M_r; ++j) {
        // radius
        double r = (j + 0.5) * g.dr;
        double rm = r - 0.5 * g.dr;
        double rp = r + 0.5 * g.dr;

        // loop over columns
        for (int i = i0; i < i1; ++i) {
            // matrix entries
            double diag    = - dx2 * (r  * eps[L](i, j) + r  * eps[R](i, j))
                             - dr2 * (rp * eps[O](i, j) + rm * eps[I](i, j));
            double left    = dx2 * r  * eps[L](i, j);
            double right   = (i > i0) ? dx2 * r  * eps[R](i - 1, j) : 0;
            double inside  = dr2 * rm * eps[I](i, j);
            double outside = (j >  0) ? dr2 * rm * eps[O](i, j - 1) : 0;

            // horizontal von Neumann boundary conditions
            if (i == 1) {
                right *= 2;
            }
            if (i == g.N_x - 1) {
                left *= 2;
            }

            // vertical von Neumann boundary conditions
            if (j == 1) {
                outside -= dr2 * 0.5 * g.dr * eps[O](i, 0);
                outside *= 2;
            }
            if (j == g.M_r - 1) {
                inside += dr2 * 0.5 * g.dr * eps[I](i, g.M_r - 1);
                inside *= 2;
            }

            // store values
            set_value(k, k, diag);
            if (i > i0) {
                set_value(k, k - 1, left);
                set_value(k - 1, k, right);
            }
            if (j > 0) {
                set_value(k, k - delta, inside);
                set_value(k - delta, k, outside);
            }

            // next diag element
            ++k;
        }

        // cut off source, drain contacts
        if (j == g.M_cnt - 1) {
            i0 = g.N_sc;
            i1 = g.N_x - g.N_dc;
            delta = g.N_x - g.N_sc;
        }
        if (j == g.M_cnt) {
            delta = g.N_x - g.N_sc - g.N_dc;
        }

        // cut off gate contact
        if (j == g.M_cnt + g.M_ox - 1) {
            i0 = g.N_sc + g.N_sox;
            i1 = i0 + g.N_sg;
            delta = g.N_x - g.N_sc - g.N_dc - g.N_sox;
        }
        if ((j == g.M_cnt + g.M_ox) && (i0 == g.N_sc + g.N_sox)) {
            delta = g.N_x - g.N_sc - g.N_dc - g.N_sox - g.N_g;
        }
        if ((j == g.M_cnt + g.M_ox) && (i0 == g.N_sc + g.N_sox + g.N_sg + g.N_g)) {
            delta = g.N_x - g.N_sc - g.N_dc - g.N_sox - g.N_g - g.N_dox;
        }
        if (j >= g.M_cnt + g.M_ox) {
            if (i0 == g.N_sc + g.N_sox) {
                i0 = i1 + g.N_g;
                i1 = i0 + g.N_dg;
                --j; // repeat i loop for second part (right side of gate)
            } else {
                i1 = i0 - g.N_g;
                i0 = i1 - g.N_sg;
            }
        }
    }

    // fit to size
    indices.resize(2, N_v);
    values.resize(N_v);

    // create, save and return sparse matrix S
    S[g.name] = sp_mat(indices, values);
    return S[g.name];
}
const std::array<arma::mat, 4> & potential::get_eps(const geometry & g) {
    // check if eps was already calculated for this geometry
    static std::map<std::string, std::array<arma::mat, 4>> eps;
    auto it = eps.find(g.name);
    if (it != std::end(eps)) {
        return it->second;
    }

    // create new eps
    using namespace arma;
    std::array<arma::mat, 4> e;

    for (int i = 0; i < 4; ++i) {
        e[i] = arma::mat(g.N_x, g.M_r);
        e[i].fill(c::eps_0);
    }

    int i, j;

    // cnt
    for (j = 0; j < g.M_cnt - 1; ++j) {
        for (i = 0; i < g.N_x; ++i) {
            e[L](i, j) = g.eps_cnt * c::eps_0;
            e[R](i, j) = g.eps_cnt * c::eps_0;
            e[I](i, j) = g.eps_cnt * c::eps_0;
            e[O](i, j) = g.eps_cnt * c::eps_0;
        }
    }

    // cnt border
    for (i = 0; i < g.N_sc; ++i) {
        e[L](i, j) = 0.5 * (1.0 + g.eps_cnt) * c::eps_0;
        e[R](i, j) = 0.5 * (1.0 + g.eps_cnt) * c::eps_0;
        e[I](i, j) = g.eps_cnt * c::eps_0;
        e[O](i, j) = c::eps_0;
    }
    e[L](i, j) = 0.5 * (1.0 + g.eps_cnt) * c::eps_0;
    e[R](i, j) = 0.5 * (g.eps_ox + g.eps_cnt) * c::eps_0;
    e[I](i, j) = g.eps_cnt * c::eps_0;
    e[O](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
    for (++i; i < g.N_x - g.N_dc - 1; ++i) {
        e[L](i, j) = 0.5 * (g.eps_ox + g.eps_cnt) * c::eps_0;
        e[R](i, j) = 0.5 * (g.eps_ox + g.eps_cnt) * c::eps_0;
        e[I](i, j) = g.eps_cnt * c::eps_0;
        e[O](i, j) = g.eps_ox * c::eps_0;
    }
    e[L](i, j) = 0.5 * (g.eps_ox + g.eps_cnt) * c::eps_0;
    e[R](i, j) = 0.5 * (1.0 + g.eps_cnt) * c::eps_0;
    e[I](i, j) = g.eps_cnt * c::eps_0;
    e[O](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
    for (++i; i < g.N_x; ++i) {
        e[L](i, j) = 0.5 * (1.0 + g.eps_cnt) * c::eps_0;
        e[R](i, j) = 0.5 * (1.0 + g.eps_cnt) * c::eps_0;
        e[I](i, j) = g.eps_cnt * c::eps_0;
        e[O](i, j) = c::eps_0;
    }

    // oxide
    for (++j; j < g.M_cnt + g.M_ox - 1; ++j) {
        i = g.N_sc;
        e[L](i, j) = c::eps_0;
        e[R](i, j) = g.eps_ox * c::eps_0;
        e[I](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
        e[O](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
        for (++i; i < g.N_x - g.N_dc - 1; ++i) {
            e[L](i, j) = g.eps_ox * c::eps_0;
            e[R](i, j) = g.eps_ox * c::eps_0;
            e[I](i, j) = g.eps_ox * c::eps_0;
            e[O](i, j) = g.eps_ox * c::eps_0;
        }
        e[L](i, j) = g.eps_ox * c::eps_0;
        e[R](i, j) = c::eps_0;
        e[I](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
        e[O](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
    }

    // oxide border
    i = g.N_sc;
    e[L](i, j) = c::eps_0;
    e[R](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
    e[I](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
    e[O](i, j) = c::eps_0;
    for (++i; i < g.N_x - g.N_dc - 1; ++i) {
        e[L](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
        e[R](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
        e[I](i, j) = g.eps_ox * c::eps_0;
        e[O](i, j) = c::eps_0;
    }
    e[L](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
    e[R](i, j) = c::eps_0;
    e[I](i, j) = 0.5 * (1.0 + g.eps_ox) * c::eps_0;
    e[O](i, j) = c::eps_0;

    // save and return eps
    eps[g.name] = e;
    return eps[g.name];
}

#endif

