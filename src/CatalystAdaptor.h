#pragma once
#include "CatalystAdaptor.h"
#include "Femocs.h"

namespace CatalystAdaptor {
    void Initialize(int numScripts, char *scripts[]);
    void Finalize();
    void CoProcess(femocs::Femocs &project, double time, unsigned int timeStep, bool lastTimeStep);
}
