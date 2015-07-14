TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += .
DEPENDPATH += .

SOURCES += main.cpp

HEADERS += \
    constant.hpp \
    device.hpp \
    util/anderson.hpp \
    util/brent.hpp \
    util/integral.hpp \
    util/inverse.hpp \
    util/fermi.hpp \
    contact.hpp \
    voltage.hpp \
    potential.hpp \
    charge_density.hpp \
    geometry.hpp \
    model.hpp \
    current.hpp \
    green.hpp \
    util/gnuplot.hpp \
    util/rwth.hpp \
    circuit.hpp \
    device_params.hpp \
    wave_packet.hpp \
    inverter.hpp \
    ring_oscillator.hpp \
    util/system.hpp \
    signal.hpp \
#    util/simd/avx.hpp \
#    util/simd/avx2.hpp \
    util/simd/fma.hpp \
    util/simd/mmx.hpp \
    util/simd/simd.hpp \
    util/simd/sse.hpp \
    util/simd/sse2.hpp \
    util/simd/sse3.hpp \
    util/simd/sse4_1.hpp \
    util/simd/sse4_2.hpp \
    util/simd/ssse3.hpp

LIBS += -lblas -lgomp -lsuperlu

QMAKE_CXXFLAGS = -std=c++14 -march=native -fopenmp
QMAKE_CXXFLAGS_RELEASE = -Ofast -fno-finite-math-only
