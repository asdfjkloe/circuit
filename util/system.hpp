#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <sstream>
#include <string>
#include <unistd.h>

static inline void shell(const std::stringstream & command) {
    system(command.str().c_str());
}

static inline void system(const std::string & s) {
    system(s.c_str());
}

static inline std::string now() {
    using namespace std;
    time_t rawtime;
    struct tm * timeinfo;
    char buf[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buf, 80, "%Y-%m-%d-%H-%M-%S", timeinfo);
    return buf;
}

static inline std::string get_login_name() {
    char buf[64];
    getlogin_r(buf, sizeof(buf));
    return buf;
}

static inline const std::string & save_folder() {
    static std::string folder;

    if (folder.empty()) {
        folder = "/home/" + get_login_name() + "/" + now();
    }

    return folder;
}

#endif
