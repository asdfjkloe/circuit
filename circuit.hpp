#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include "device.hpp"

class circuit {
public:
    inline circuit();

    inline const device & operator[](int index) const;
    inline device & operator[](int index);

    inline int add_device(const device & d);

    inline void link(int d1, int ct1, int d2, int ct2);

private:
    std::vector<device> devices;
};

//----------------------------------------------------------------------------------------------------------------------

circuit::circuit() {
}

const device & circuit::operator[](int index) const {
    return devices[index];
}

device & circuit::operator[](int index) {
    return devices[index];
}

int circuit::add_device(const device & d) {
    int n = devices.size();
    devices.push_back(d);
    return n;
}

void circuit::link(int d1, int ct1, int d2, int ct2) {
    devices[d2].contacts[ct2] = devices[d1].contacts[ct1];
}

#endif
