#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <string>

struct Nucleus
{
    Nucleus();
    Nucleus(std::string configFileName);

    std::string Name;
    unsigned int Z;
    unsigned int A;
};

#endif /* NUCLEUS_H */
