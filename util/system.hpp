#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <sstream>
#include <string>
#include <cstdlib>

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

static inline const std::string & save_folder(bool timestamp = true, const std::string & prefix = "td_results") {
    static std::string folder;

    if (folder.empty()) {
        folder = std::string(std::getenv("HOME")) + "/" + prefix + (timestamp ? "_" + now() : "");
    }

    return folder;
}

#endif
