#pragma once

#include "CatalystAdaptor.h"
#include "Femocs.h"
#include <string>

namespace CatalystAdaptor {
    void Initialize(const char *path);

    void Finalize();

    void CoProcess(femocs::Femocs &project, double time, unsigned int timeStep, bool lastTimeStep);
}
