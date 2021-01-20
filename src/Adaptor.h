#pragma once
#include "Adaptor.h"
#include "Femocs.h"

namespace Adaptor {
    void Initialize(int numScripts, char *scripts[]);
    void Finalize();
    void CoProcess(femocs::Femocs &project, double time, unsigned int timeStep, bool lastTimeStep);
}
