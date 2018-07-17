#include <map>
#include <iostream>
#include "util.h"
#include <stdlib.h>

int main() {
    for (auto it = samples.begin(); it != samples.end(); it++) {
        if (it->first.find("VBF") != std::string::npos)
            std::system(("./Skim " + it->first).c_str());
    }
    return 0;
}
