#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#include <armadillo>
#include <map>

#include "charge_density.hpp"
#include "constant.hpp"
#include "device_params.hpp"
#include "voltage.hpp"
#include "util/anderson.hpp"

class potential {
public:
    arma::vec data;
    arma::vec twice;

    inline potential();
    inline potential(const device_params & p, const arma::vec & R);
    inline potential(const device_params & p, const arma::vec & R0, const charge_density & n);
    inline potential(const device_params & p, const voltage<3> & V);
    inline double update(const device_params & p, const arma::vec & R0, const charge_density & n, anderson & mr_neo);

    inline double & operator()(int index);
    inline const double & operator()(int index) const;
    inline double s() const;
    inline double d() const;

    inline void smooth(const device_params & p);

    static inline arma::vec get_R0(const device_params & p, const voltage<3> & V);

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

    static inline arma::vec get_R(const device_params & p, const arma::vec & R0, const charge_density & n);
    static inline const arma::sp_mat & get_S(const device_params & p);
    static inline const std::array<arma::mat, 4> & get_eps(const device_params & p);
};

//----------------------------------------------------------------------------------------------------------------------

potential::potential() {
}
potential::potential(const device_params & p, const arma::vec & R) {
    using namespace arma;

    // superlu is apparantly not thread-safe
    vec phi2D;
    //#pragma omp critical
    //{
        phi2D = spsolve(get_S(p), R);
    //}

    data = phi2D.rows((p.M_cnt - 1) * p.N_x, p.M_cnt * p.N_x - 1);

    update_twice();
}
potential::potential(const device_params & p, const arma::vec & R0, const charge_density & n)
    : potential(p, get_R(p, R0, n)) {
}
potential::potential(const device_params & p, const voltage<3> & V)
    : potential(p, get_R0(p, V)) {
}

double potential::update(const device_params & p, const arma::vec & R0, const charge_density & n, anderson & mr_neo) {
    // calculate undamped update
    potential phi_update(p, R0, n);

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

void potential::smooth(const device_params & p) {
    // smooth source region
    if (p.F[S] > 0) {
        smooth<true>(0, p.N_sc + p.N_sox + p.N_sg + p.N_g * 0.05);
    } else {
        smooth<false>(0, p.N_sc + p.N_sox + p.N_sg + p.N_g * 0.05);
    }

    // smooth drain region
    if (p.F[D] > 0) {
        smooth<true>(p.N_sc + p.N_sox + p.N_sg + p.N_g * 0.95, p.N_x);
    } else {
        smooth<false>(p.N_sc + p.N_sox + p.N_sg + p.N_g * 0.95, p.N_x);
    }

    update_twice();
}

arma::vec potential::get_R0(const device_params & p, const voltage<3> & V0) {
    using namespace arma;

    // shortcuts
    double dr2 = 1.0 / p.dr / p.dr;
    double dx2 = 1.0 / p.dx / p.dx;

    // get reference to 4 eps matrices (one for each direction)
    const std::array<arma::mat, 4> & eps = get_eps(p);

    // add built-in voltages
    voltage<3> V = - (V0 + p.F);

    // init return vector
    vec R0(p.N_x * p.M_r);
    R0.fill(0.0);

    // create R0 from geometry and voltages
    int i, j, k;
    double r, rp;

    // add dirichlet boundary conditions where necessary
    j = p.M_cnt - 1;
    r = j * p.dr + 0.5 * p.dr;
    rp = r + 0.5 * p.dr;
    k = j * p.N_x;
    for (i = 0; i < p.N_sc; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[S];
    }
    k = j * p.N_x + p.N_x - p.N_dc;
    for (i = p.N_x - p.N_dc; i < p.N_x; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[D];
    }
    for (j = p.M_cnt; j < p.M_cnt + p.M_ox - 1; ++j) {
        r = j * p.dr + 0.5 * p.dr;
        R0(k) -= dx2 * r * eps[L](p.N_sc, j) * V[S];
        k += p.N_sox + p.N_sg + p.N_g + p.N_dg + p.N_dox - 1;
        R0(k++) -= dx2 * r * eps[R](p.N_x - p.N_dc - 1, j) * V[D];
    }
    r = j * p.dr + 0.5 * p.dr;
    rp = r + 0.5 * p.dr;
    R0(k) -= dx2 * r * eps[L](p.N_sc, j) * V[S];
    for (i = p.N_sc; i < p.N_sc + p.N_sox; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[S];
    }
    k += p.N_sg;
    for (i = p.N_sc + p.N_sox + p.N_sg; i < p.N_sc + p.N_sox + p.N_sg + p.N_g; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[G];
    }
    k += p.N_dg;
    for (i = p.N_x - p.N_dc - p.N_dox; i < p.N_x - p.N_dc; ++i) {
        R0(k++) -= dr2 * rp * eps[O](i, j) * V[D];
    }
    R0(k - 1) -= dx2 * r * eps[R](p.N_x - p.N_dc - 1, j) * V[D];

    for (j = p.M_cnt + p.M_ox; j < p.M_r; ++j) {
        r = j * p.dr + 0.5 * p.dr;
        R0(k) -= dx2 * r * eps[L](p.N_sc + p.N_sox, j) * V[S];
        k += p.N_sg - 1;
        R0(k++) -= dx2 * r * eps[R](p.N_sc + p.N_sox + p.N_sg - 1, j) * V[G];
        R0(k) -= dx2 * r * eps[L](p.N_sc + p.N_sox + p.N_sg + p.N_g, j) * V[G];
        k += p.N_dg - 1;
        R0(k++) -= dx2 * r * eps[R](p.N_x - p.N_dc - p.N_dox - 1, j) * V[D];
    }

    // shrink to fit and return
    R0.resize(k);
    return R0;
}

void potential::update_twice() {
    twice.resize(data.size() * 2);

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
                    data.rows(i+1, j-1).fill(data(i));
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
                    data.rows(j+1, i-1).fill(data(i));
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
                    data.rows(i+1, j-1).fill(data(i));
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
                    data.rows(j+1, i - 1).fill(data(i));
                    break;
                }
            }
        }
    }
}

arma::vec potential::get_R(const device_params & p, const arma::vec & R0, const charge_density & n) {
    using namespace arma;
    vec R = R0;

    // add non-constant part to R0 due to charge_density at the boundary of the cnt
    R.rows((p.M_cnt - 1) * p.N_x, p.M_cnt * p.N_x - 1) += n.total * p.r_cnt * 1e9; // 10^9 because of m->nm in eps_0
    return R;
}
const arma::sp_mat & potential::get_S(const device_params & p) {
    // check if S was already calculated for this device
    static thread_local std::map<std::string, arma::sp_mat> S;
    auto it = S.find(p.name);
    if (it != std::end(S)) {
        return it->second;
    }

    // calculate new S
    using namespace arma;

    // shortcuts
    double dr2 = 1.0 / p.dr / p.dr;
    double dx2 = 1.0 / p.dx / p.dx;

    // get reference to 4 eps matrices (one for each direction)
    const std::array<arma::mat, 4> & eps = get_eps(p);

    // helper objects to create sparse matrix faster
    umat indices(2, p.M_r * p.N_x * 5);
    vec values(     p.M_r * p.N_x * 5);
    uword N_v = 0;
    auto set_value = [&indices, &values, &N_v] (uword i, uword j, double val) {
        indices(0, N_v) = i;
        indices(1, N_v) = j;
        values(N_v++) = val;
    };

    // start values
    int k     = 0;     // current main diagonal element
    int i0    = 0;     // left i limit
    int i1    = p.N_x; // right i limit
    int delta = p.N_x; // distance to next vertical coupling off diagonal

    // main loop over rows
    for (int j = 0; j < p.M_r; ++j) {
        // radius
        double r = (j + 0.5) * p.dr;
        double rm = r - 0.5 * p.dr;
        double rp = r + 0.5 * p.dr;

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
            if (i == p.N_x - 1) {
                left *= 2;
            }

            // vertical von Neumann boundary conditions
            if (j == 1) {
                outside -= dr2 * 0.5 * p.dr * eps[O](i, 0);
                outside *= 2;
            }
            if (j == p.M_r - 1) {
                inside += dr2 * 0.5 * p.dr * eps[I](i, p.M_r - 1);
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
        if (j == p.M_cnt - 1) {
            i0 = p.N_sc;
            i1 = p.N_x - p.N_dc;
            delta = p.N_x - p.N_sc;
        }
        if (j == p.M_cnt) {
            delta = p.N_x - p.N_sc - p.N_dc;
        }

        // cut off gate contact
        if (j == p.M_cnt + p.M_ox - 1) {
            i0 = p.N_sc + p.N_sox;
            i1 = i0 + p.N_sg;
            delta = p.N_x - p.N_sc - p.N_dc - p.N_sox;
        }
        if ((j == p.M_cnt + p.M_ox) && (i0 == p.N_sc + p.N_sox)) {
            delta = p.N_x - p.N_sc - p.N_dc - p.N_sox - p.N_g;
        }
        if ((j == p.M_cnt + p.M_ox) && (i0 == p.N_sc + p.N_sox + p.N_sg + p.N_g)) {
            delta = p.N_x - p.N_sc - p.N_dc - p.N_sox - p.N_g - p.N_dox;
        }
        if (j >= p.M_cnt + p.M_ox) {
            if (i0 == p.N_sc + p.N_sox) {
                i0 = i1 + p.N_g;
                i1 = i0 + p.N_dg;
                --j; // repeat i loop for second part (right side of gate)
            } else {
                i1 = i0 - p.N_g;
                i0 = i1 - p.N_sg;
            }
        }
    }

    // fit to size
    indices.resize(2, N_v);
    values.resize(N_v);

    // create, save and return sparse matrix S
    S[p.name] = sp_mat(indices, values);
    return S[p.name];
}
const std::array<arma::mat, 4> & potential::get_eps(const device_params & p) {
    // check if eps was already calculated for this geometry
    static thread_local std::map<std::string, std::array<arma::mat, 4>> eps;
    auto it = eps.find(p.name);
    if (it != std::end(eps)) {
        return it->second;
    }

    // create new eps
    using namespace arma;
    std::array<arma::mat, 4> e;

    for (int i = 0; i < 4; ++i) {
        e[i] = arma::mat(p.N_x, p.M_r);
        e[i].fill(c::eps_0);
    }

    int i, j;

    // cnt
    for (j = 0; j < p.M_cnt - 1; ++j) {
        for (i = 0; i < p.N_x; ++i) {
            e[L](i, j) = p.eps_cnt * c::eps_0;
            e[R](i, j) = p.eps_cnt * c::eps_0;
            e[I](i, j) = p.eps_cnt * c::eps_0;
            e[O](i, j) = p.eps_cnt * c::eps_0;
        }
    }

    // cnt border
    for (i = 0; i < p.N_sc; ++i) {
        e[L](i, j) = 0.5 * (1.0 + p.eps_cnt) * c::eps_0;
        e[R](i, j) = 0.5 * (1.0 + p.eps_cnt) * c::eps_0;
        e[I](i, j) = p.eps_cnt * c::eps_0;
        e[O](i, j) = c::eps_0;
    }
    e[L](i, j) = 0.5 * (1.0 + p.eps_cnt) * c::eps_0;
    e[R](i, j) = 0.5 * (p.eps_ox + p.eps_cnt) * c::eps_0;
    e[I](i, j) = p.eps_cnt * c::eps_0;
    e[O](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
    for (++i; i < p.N_x - p.N_dc - 1; ++i) {
        e[L](i, j) = 0.5 * (p.eps_ox + p.eps_cnt) * c::eps_0;
        e[R](i, j) = 0.5 * (p.eps_ox + p.eps_cnt) * c::eps_0;
        e[I](i, j) = p.eps_cnt * c::eps_0;
        e[O](i, j) = p.eps_ox * c::eps_0;
    }
    e[L](i, j) = 0.5 * (p.eps_ox + p.eps_cnt) * c::eps_0;
    e[R](i, j) = 0.5 * (1.0 + p.eps_cnt) * c::eps_0;
    e[I](i, j) = p.eps_cnt * c::eps_0;
    e[O](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
    for (++i; i < p.N_x; ++i) {
        e[L](i, j) = 0.5 * (1.0 + p.eps_cnt) * c::eps_0;
        e[R](i, j) = 0.5 * (1.0 + p.eps_cnt) * c::eps_0;
        e[I](i, j) = p.eps_cnt * c::eps_0;
        e[O](i, j) = c::eps_0;
    }

    // oxide
    for (++j; j < p.M_cnt + p.M_ox - 1; ++j) {
        i = p.N_sc;
        e[L](i, j) = c::eps_0;
        e[R](i, j) = p.eps_ox * c::eps_0;
        e[I](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
        e[O](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
        for (++i; i < p.N_x - p.N_dc - 1; ++i) {
            e[L](i, j) = p.eps_ox * c::eps_0;
            e[R](i, j) = p.eps_ox * c::eps_0;
            e[I](i, j) = p.eps_ox * c::eps_0;
            e[O](i, j) = p.eps_ox * c::eps_0;
        }
        e[L](i, j) = p.eps_ox * c::eps_0;
        e[R](i, j) = c::eps_0;
        e[I](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
        e[O](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
    }

    // oxide border
    i = p.N_sc;
    e[L](i, j) = c::eps_0;
    e[R](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
    e[I](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
    e[O](i, j) = c::eps_0;
    for (++i; i < p.N_x - p.N_dc - 1; ++i) {
        e[L](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
        e[R](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
        e[I](i, j) = p.eps_ox * c::eps_0;
        e[O](i, j) = c::eps_0;
    }
    e[L](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
    e[R](i, j) = c::eps_0;
    e[I](i, j) = 0.5 * (1.0 + p.eps_ox) * c::eps_0;
    e[O](i, j) = c::eps_0;

    // save and return eps
    eps[p.name] = e;
    return eps[p.name];
}

#endif

