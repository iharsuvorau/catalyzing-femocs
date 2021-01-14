#include "Adaptor.h"

#include <cstdlib>
#include <iterator>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

namespace {
vtkCPProcessor *Processor = NULL;
vtkUnstructuredGrid *vtkGrid;

void BuildVTKGrid(Grid &grid) {
  // ...
  return;
}

void UpdateVTKAttributes(Grid &grid, Attributes &attributes,
                         vtkCPInputDataDescription *idd) {
  // ...
  return;
}

void BuildVTKDataStructures(Grid &grid, Attributes &attributes,
                            vtkCPInputDataDescription *idd) {
  // ...
  return;
}
} // namespace

namespace Adaptor {
void Initialize(int numScripts, char *scripts[]) {
  if (Processor == NULL) {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
  } else {
    Processor->RemoveAllPipelines();
  }

  for (int i = 0; i < numScripts; i++) {
    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize(scripts[i]);
    Processor->AddPipeline(pipeline.GetPointer());
  }
}

void Finalize() {
  if (Processor) {
    Processor->Delete();
    Processor = NULL;
  }

  if (vtkGrid) {
    vtkGrid->Delete();
    vtkGrid = NULL;
  }
}

void CoProcess(Grid &grid, Attributes &attributes, double time,
               unsigned int timeStep, bool lastTimeStep) {
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");
  dataDescription->SetTimeData(time, timeStep);

  if (lastTimeStep == true) {
    dataDescription->ForceOutputOn();
  }

  if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
    vtkCPInputDataDescription *idd =
        dataDescription->GetInputDescriptionByName("input");
    BuildVTKDataStructures(grid, attributes, idd);
    idd->SetGrid(vtkGrid);
    Processor->CoProcess(dataDescription.GetPointer());
  }
}
} // namespace Adaptor

Grid::Grid() {}

void Grid::Initialize(const unsigned int numPoints[3],
                      const double spacing[3]) {
  // ...
}

size_t Grid::GetNumberOfPoints() { return this->Points.size() / 3; }

size_t Grid::GetNumberOfCells() { return this->Cells.size() / 8; }

double *Grid::GetPointsArray() {
  if (this->Points.empty()) {
    return NULL;
  }
  return &(this->Points[0]);
}

double *Grid::GetPoint(size_t pointID) {
  if (pointID >= this->Points.size()) {
    return NULL;
  }
  return &(this->Points[pointID * 3]);
}

unsigned int *Grid::GetCellPoints(size_t cellID) {
  if (cellID >= this->Cells.size()) {
    return NULL;
  }
  return &(this->Cells[cellID * 8]);
}

Attributes::Attributes() { this->GridPtr = NULL; }

void Attributes::Initialize(Grid *grid) { this->GridPtr = grid; }

void Attributes::UpdateFields(double time) {
  // ...
}

double *Attributes::GetVelocityArray() {
  if (this->Velocity.empty()) {
    return NULL;
  }
  return &this->Velocity[0];
}

float *Attributes::GetPressureArray() {
  if (this->Pressure.empty()) {
    return NULL;
  }
  return &this->Pressure[0];
}
