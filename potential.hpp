#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#include <armadillo>
#include <map>

#include "charge_density.hpp"
#include "constant.hpp"
#include "device_params.hpp"
#include "voltage.hpp"
#include "util/anderson.hpp"
#include "util/gnuplot.hpp"

class potential {
public:
    arma::vec data;
    arma::vec twice;

    inline potential();
    inline potential(const device_params & p, const arma::vec & R);
    inline potential(const device_params & p, const arma::vec & R0, const charge_density & n);
    inline potential(const device_params & p, const voltage<3> & V);
    inline potential(const device_params & p, const voltage<3> & V, const charge_density & n);
    inline double update(const device_params & p, const arma::vec & R0, const charge_density & n, anderson & mr_neo);

    inline double & operator()(int index);
    inline const double & operator()(int index) const;
    inline double s() const;
    inline double d() const;

    inline void smooth(const device_params & p);

    static inline arma::vec get_R0(const device_params & p, const voltage<3> & V);

    static inline void plot2D(const device_params & p, const voltage<3> & V, const charge_density & n);
    static inline void plot2D(const device_params & p, const voltage<3> & V);

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
    #pragma omp critical
    {
        phi2D = spsolve(get_S(p), R);
    }

    data = phi2D.rows((p.M_cnt - 1) * p.N_x, p.M_cnt * p.N_x - 1);

    update_twice();
}
potential::potential(const device_params & p, const arma::vec & R0, const charge_density & n)
    : potential(p, get_R(p, R0, n)) {
}
potential::potential(const device_params & p, const voltage<3> & V)
    : potential(p, get_R0(p, V)) {
}
potential::potential(const device_params & p, const voltage<3> & V, const charge_density & n)
    : potential(p, get_R(p, get_R0(p, V), n)) {
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

void potential::plot2D(const device_params & p, const voltage<3> & V, const charge_density & n) {
    // solve poisson's
    arma::vec phi2Dv = spsolve(get_S(p), get_R(p, get_R0(p, V), n));

    arma::mat phi2D(p.N_x, p.M_r);
    int k = 0;
    for (int j = 0; j < p.M_r; ++j) {
        for (int i = 0; i < p.N_x; ++i) {
            if (j < p.M_cnt) {
                phi2D(i, j) = phi2Dv(k++);
            } else if (j < p.M_cnt + p.M_ox) {
                if (i < p.N_sc) {
                    phi2D(i, j) = -(V[S] + p.F[S]);
                } else if (i >= p.N_x - p.N_dc) {
                    phi2D(i, j) = -(V[D] + p.F[D]);
                } else {
                    phi2D(i, j) = phi2Dv(k++);
                }
            } else {
                if (i < p.N_sc + p.N_sox) {
                    phi2D(i, j) = -(V[S] + p.F[S]);
                } else if (i >= p.N_x - p.N_dc - p.N_dox) {
                    phi2D(i, j) = -(V[D] + p.F[D]);
                } else if ((i >= p.N_sc + p.N_sox + p.N_sg) && (i < p.N_sc + p.N_sox + p.N_sg + p.N_g)) {
                    phi2D(i, j) = -(V[G]+ p.F[G]);
                } else {
                    phi2D(i, j) = phi2Dv(k++);
                }
            }
        }
    }
    phi2D = arma::join_horiz(arma::fliplr(phi2D), phi2D).t();

    gnuplot gp;

    // gnuplot setup
    gp << "set palette defined ( 0 '#D73027', 1 '#F46D43', 2 '#FDAE61', 3 '#FEE090', 4 '#E0F3F8', 5 '#ABD9E9', 6 '#74ADD1', 7 '#4575B4' )\n";
    gp << "set title \"Potential cross-section for V_{s} = " << V[S] << " V, V_{g} = " << V[G] << " V and V_{d} = " << V[D] << " V\"\n";
    gp << "set xlabel \"x / nm\"\n";
    gp << "set ylabel \"r / nm\"\n";
    gp << "set zlabel \"Phi / V\"\n";
    gp << "unset key\n";

//    gp << "set terminal pdf\nset output 'potential2D.pdf'\n";

    // indicate cnt area
    gp << "set obj rect from " << 0 << "," << p.r_cnt << " to " << p.l << "," << -p.r_cnt << "front fillstyle empty\n";
    gp << "set label \"CNT\" at " << 0.5 * p.l << "," << 0 << " center front\n";

    // indicate oxide area
    double x0 = p.l_sc;
    double x1 = p.l - p.l_dc;
    double y0 = p.r_cnt + p.d_ox;
    double y1 = p.r_cnt;
    gp << "set obj rect from " << x0 << "," << y0 << " to " << x1 << "," << y1 << "front fillstyle empty\n";
    gp << "set label \"gate oxide\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * (y1 - y0) - y1<< " center front\n";
    gp << "set obj rect from " << x0 << "," << -y1 << " to " << x1 << "," << -y0 << "front fillstyle empty\n";
    gp << "set label \"gate oxide\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * -(y1 - y0) + y1 << " center front\n";

    // indicate gate contact area
    x0 = p.l_sc + p.l_sox + p.l_sg;
    x1 = p.l - p.l_dc - p.l_dox - p.l_dg;
    y0 = p.R;
    y1 = p.r_cnt + p.d_ox;
    gp << "set obj rect from " << x0 << "," << y0 << " to " << x1 << "," << y1 << "front fillstyle empty\n";
    gp << "set label \"gate\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * (y1 - y0) - y1<< " center front\n";
    gp << "set obj rect from " << x0 << "," << -y1 << " to " << x1 << "," << -y0 << "front fillstyle empty\n";
    gp << "set label \"gate\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * -(y1 - y0) + y1 << " center front\n";

    // indicate left contact areas
    // --------------------------------------------------------------------------------------
    // top/bottom
    x0 = 0 + p.dx/2;
    x1 = p.l_sc + p.l_sox + p.dx/2;
    y0 = p.R - p.dr/2;
    y1 = y0;
    gp << "set arrow from " << x0 << "," << y0 << " to "<< x1 << "," << y1 << " nohead front\n";
    gp << "set arrow from " << x0 << "," << -y0 << " to "<< x1 << "," << -y1 << " nohead front\n"; //mirror y
    x0 = p.l - p.dx/2;
    x1 = p.l - x1 - p.l_dc - p.l_dox - p.dx/2;
    gp << "set arrow from " << x0 << "," << y0 << " to "<< x1 << "," << y1 << " nohead front\n"; //mirror x
    gp << "set arrow from " << x0 << "," << -y0 << " to "<< x1 << "," << -y1 << " nohead front\n"; //mirror both

    // outside
    x0 = 0;
    y0 = p.R - p.dr/2;
    x1 = x0;
    y1 = p.r_cnt;
    gp << "set arrow from " << x0 << "," << y0 << " to "<< x1 << "," << y1 << " nohead front\n";
    gp << "set arrow from " << x0 << "," << -y0 << " to "<< x1 << "," << -y1 << " nohead front\n"; //mirror
    x0 = p.l - p.dx/2;
    x1 = x0;
    gp << "set arrow from " << x0 << "," << y0 << " to "<< x1 << "," << y1 << " nohead front\n"; //mirror x
    gp << "set arrow from " << x0 << "," << -y0 << " to "<< x1 << "," << -y1 << " nohead front\n"; //mirror both

    // inside
    x0 = p.l_sc + p.l_sox;
    y0 = p.R - p.dr/2;
    x1 = x0;
    y1 = p.r_cnt + p.d_ox;
    gp << "set arrow from " << x0 << "," << y0 << " to "<< x1 << "," << y1 << " nohead front\n";
    gp << "set arrow from " << x0 << "," << -y0 << " to "<< x1 << "," << -y1 << " nohead front\n"; //mirror
    x0 = p.l - p.l_dc - p.l_dox;
    x1 = x0;
    gp << "set arrow from " << x0 << "," << y0 << " to "<< x1 << "," << y1 << " nohead front\n"; //mirror x
    gp << "set arrow from " << x0 << "," << -y0 << " to "<< x1 << "," << -y1 << " nohead front\n"; //mirror both

    // label
    x0 = 0;
    x1 = p.l_sc + p.l_sox;
    y0 = p.r_cnt + p.d_ox;
    y1 = p.R - p.dr/2;
    gp << "set label \"source\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * -(y1 - y0) + y1 << " center front\n";
    gp << "set label \"source\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * +(y1 - y0) - y1 << " center front\n"; //mirror y
    x0 = p.l;
    x1 = p.l - p.l_dc - p.l_dox;
    gp << "set label \"drain\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * -(y1 - y0) + y1 << " center front\n"; //mirror x
    gp << "set label \"drain\" at " << 0.5 * (x1 - x0) + x0 << "," << 0.5 * +(y1 - y0) - y1 << " center front\n"; //mirror both

    gp.set_background(p.x, arma::join_vert(arma::flipud(-p.r), p.r), phi2D);
    gp.plot();
}
void potential::plot2D(const device_params & p, const voltage<3> & V) {
    charge_density n;
    n.total.resize(p.N_x);
    n.total.fill(0.0);
    plot2D(p, V, n);
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
