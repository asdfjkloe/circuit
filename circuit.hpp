#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include "constant.hpp"
#include "contact.hpp"
#include "device.hpp"

template<int N_in>
class circuit {
public:
    inline circuit();

    inline const device & operator[](int index) const;
    inline device & operator[](int index);

    inline const contact_ptr & contact(int index) const;
    inline contact_ptr & contact(int index);

    inline int add_device(const device & d);

    inline void link(int d1, int ct1, int ct2);
    inline void link(int d1, int ct1, int d2, int ct2);

    inline bool steady_state();

private:
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
const contact_ptr & circuit<N_in>::contact(int index) const {
    return contacts[index];
}

template<int N_in>
contact_ptr & circuit<N_in>::contact(int index) {
    return contacts[index];
}

template<int N_in>
int circuit<N_in>::add_device(const device & d) {
    int n = devices.size();
    devices.push_back(d);
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
bool circuit<N_in>::steady_state() {

}

#endif
