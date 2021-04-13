#pragma once

#include "CatalystAdaptor.h"
#include "Femocs.h"
#include <string>

namespace CatalystAdaptor {
    void Initialize(const char *path, const char *meshCellType);

    void Finalize();

    void ImportAtoms(double **atoms, int numAtoms);

    void CoProcess(femocs::Femocs &project,
                   const double time,
                   const unsigned int timeStep,
                   const bool lastTimeStep);
}
