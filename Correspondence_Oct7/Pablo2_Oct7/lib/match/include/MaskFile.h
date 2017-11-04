#ifndef MASK_FILE_H
#define MASK_FILE_H

#include "Match.h"
#include "Registry.h"

class MaskFile
{
public:
    MaskFile() {}

    bool read(const char * filename, Match * match);
    void write(const char * filename, Match * match);

private:
    Registry registry;
};

#endif

