#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include "device.hpp"

template<int N_in, int N_out>
class circuit {
public:
    inline circuit();

    inline const device & operator[](int index) const;
    inline device & operator[](int index);

    inline const contact_ptr & get_input(int index) const;
    inline contact_ptr & get_input(int index);
    inline const contact_ptr & get_output(int index) const;
    inline contact_ptr & get_output(int index);

    inline int add_device(const std::string & n, const device_params & p);

    inline void link(int d1, int ct1, int d2, int ct2);
    inline void link_input(int d, int ct, int in);
    inline void link_output(int d, int ct, int out);
    inline void link_inout(int in, int out);
    inline void link_outin(int out, int in);

    inline virtual bool steady_state(const std::array<double, N_in> & V) = 0;
    inline bool time_step(const std::array<double, N_in> & V);

    template<bool plots>
    inline void save();

protected:
    std::vector<device> devices;
    std::array<contact_ptr, N_in> inputs;
    std::array<contact_ptr, N_out> outputs;
    std::array<std::vector<double>, N_out> V_out;
};

//----------------------------------------------------------------------------------------------------------------------

template<int N_in, int N_out>
circuit<N_in, N_out>::circuit() {
    for (int i = 0; i < N_in; ++i) {
        inputs[i] = std::make_shared<contact>(0.0, c::inf);
    }
    for (int i = 0; i < N_out; ++i) {
        outputs[i] = std::make_shared<contact>(0.0, c::inf);
    }
}

template<int N_in, int N_out>
const device & circuit<N_in, N_out>::operator[](int index) const {
    return devices[index];
}
template<int N_in, int N_out>
device & circuit<N_in, N_out>::operator[](int index) {
    return devices[index];
}

template<int N_in, int N_out>
const contact_ptr & circuit<N_in, N_out>::get_input(int index) const {
    return inputs[index];
}
template<int N_in, int N_out>
contact_ptr & circuit<N_in, N_out>::get_input(int index) {
    return inputs[index];
}
template<int N_in, int N_out>
const contact_ptr & circuit<N_in, N_out>::get_output(int index) const {
    return outputs[index];
}
template<int N_in, int N_out>
contact_ptr & circuit<N_in, N_out>::get_output(int index) {
    return outputs[index];
}

template<int N_in, int N_out>
int circuit<N_in, N_out>::add_device(const std::string & n, const device_params & p) {
    int s = devices.size();
    devices.emplace_back(n, p);
    return s;
}

template<int N_in, int N_out>
void circuit<N_in, N_out>::link(int d1, int ct1, int d2, int ct2) {
    devices[d1].contacts[ct1] = devices[d2].contacts[ct2];
}
template<int N_in, int N_out>
void circuit<N_in, N_out>::link_input(int d, int ct, int in) {
    devices[d].contacts[ct] = inputs[in];
}
template<int N_in, int N_out>
void circuit<N_in, N_out>::link_output(int d, int ct, int out) {
    devices[d].contacts[ct] = outputs[out];
}
template<int N_in, int N_out>
void circuit<N_in, N_out>::link_inout(int in, int out) {
    inputs[in] = outputs[out];
}
template<int N_in, int N_out>
void circuit<N_in, N_out>::link_outin(int out, int in) {
    outputs[out] = inputs[in];
}

template<int N_in, int N_out>
bool circuit<N_in, N_out>::time_step(const std::array<double, N_in> & V) {
    // set input voltages
    for (int i = 0; i < N_in; ++i) {
        inputs[i]->V = V[i];
    }

    // calculate time_step for each device
    bool converged = true;
    for (int i = 0; i < devices.size(); ++i) {
        converged &= devices[i].time_step();
    }

    // update device contacts
    for (int i = 0; i < devices.size(); ++i) {
        devices[i].update_contacts();
    }

    // save output voltages
    for (int i = 0; i < N_out; ++i) {
        V_out[i].push_back(outputs[i]->V);
    }

    return converged;
}

template<int N_in, int N_out>
template<bool plots>
void circuit<N_in, N_out>::save() {
    for (int i = 0; i < devices.size(); ++i) {
        devices[i].save<plots>();
    }

    for (int i = 0; i < N_out; ++i) {
        arma::vec V(V_out[i]);
        std::stringstream ss;
        ss << save_folder() << "/V_out" << i;
        std::string file_name = ss.str();
        V.save(file_name + ".arma");

        if (plots) {
            // time vector
            arma::vec t = arma::linspace(0, V.size() * c::dt, V.size());

            // make a plot of V and save it as a PNG
            gnuplot gp;
            gp << "set terminal png\n";
            gp << "set title 'Output voltage'\n";
            gp << "set xlabel 't / ps'\n";
            gp << "set ylabel 'V_{out" << i << "} / V'\n";
            gp << "set format x '%1.2f'\n";
            gp << "set format y '%1.2f'\n";
            gp << "set output '" << file_name << ".png'\n";
            gp.add(std::make_pair(t * 1e12, V));
            gp.plot();
        }
    }
}

/*template<bool plots>
void inverter::save() {
    n().save<plots>();
    p().save<plots>();

    arma::vec t = arma::linspace(0, V_out.size() * c::dt, V_out.size());

    V_out.save(save_folder() + "/V_out.arma");

    std::ofstream capacitance_file(save_folder() + "/C.txt");
    capacitance_file << contacts[D]->c << std::endl;
    capacitance_file.close();

    if (plots) {

    }
}*/

#endif
