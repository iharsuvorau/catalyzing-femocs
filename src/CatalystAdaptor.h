#pragma once

#include "CatalystAdaptor.h"
#include "Femocs.h"
#include <string>

namespace CatalystAdaptor {
    void Initialize(const char *path);

    void Finalize();

    void CoProcess(femocs::Femocs &project, const double time, const unsigned int timeStep, const bool lastTimeStep);
}
