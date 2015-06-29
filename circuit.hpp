#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include "device.hpp"

template<int N_in>
class circuit {
public:
    inline circuit();

    inline const device & operator[](int index) const;
    inline device & operator[](int index);

    inline const contact_ptr & get_contact(int index) const;
    inline contact_ptr & get_contact(int index);

    inline int add_device(const device_params & p);

    inline void link(int d1, int ct1, int ct2);
    inline void link(int d1, int ct1, int d2, int ct2);

    inline virtual bool steady_state(const std::array<double, N_in> & V) = 0;
    inline bool time_step(const std::array<double, N_in> & V);

protected:
    std::array<contact_ptr, N_in> contacts;
    std::vector<device> devices;
};

//----------------------------------------------------------------------------------------------------------------------

template<int N_in>
circuit<N_in>::circuit() {
    for (int i = 0; i < N_in; ++i) {
        contacts[i] = std::make_shared<contact>(0.0, c::inf);
    }
}

template<int N_in>
const device & circuit<N_in>::operator[](int index) const {
    return devices[index];
}
template<int N_in>
device & circuit<N_in>::operator[](int index) {
    return devices[index];
}

template<int N_in>
const contact_ptr & circuit<N_in>::get_contact(int index) const {
    return contacts[index];
}

template<int N_in>
contact_ptr & circuit<N_in>::get_contact(int index) {
    return contacts[index];
}

template<int N_in>
int circuit<N_in>::add_device(const device_params & p) {
    int n = devices.size();
    devices.push_back(p);
    return n;
}

template<int N_in>
void circuit<N_in>::link(int d1, int ct1, int ct2) {
    devices[d1].contacts[ct1] = contacts[ct2];
}

template<int N_in>
void circuit<N_in>::link(int d1, int ct1, int d2, int ct2) {
    devices[d1].contacts[ct1] = devices[d2].contacts[ct2];
}

template<int N_in>
bool circuit<N_in>::time_step(const std::array<double, N_in> & V) {
    // save voltage
    for (int i = 0; i < N_in; ++i) {
        contacts[i]->V = V[i];
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

    return converged;
}

#endif
