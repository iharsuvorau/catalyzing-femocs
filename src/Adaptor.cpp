#include "Adaptor.h"
#include "Femocs.h"

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

//namespace {
//    vtkCPProcessor *Processor = NULL;
//    vtkUnstructuredGrid *vtkGrid;
//
//    void BuildVTKGrid(Grid &grid) {
//        // ...
//        return;
//    }
//
//    void UpdateVTKAttributes(Grid &grid, Attributes &attributes,
//                             vtkCPInputDataDescription *idd) {
//        // ...
//        return;
//    }
//
//    void BuildVTKDataStructures(Grid &grid, Attributes &attributes,
//                                vtkCPInputDataDescription *idd) {
//        // ...
//        return;
//    }
//} // namespace

namespace Adaptor {
    vtkCPProcessor *Processor = NULL;
    vtkSmartPointer <vtkUnstructuredGrid> grid;

    void Initialize(int numScripts, char *scripts[]) {
        if (grid == NULL) {
            grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        } else {
            grid->Initialize();
        }

        if (Processor == NULL) {
            Processor = vtkCPProcessor::New();
            Processor->Initialize();
        } else {
            Processor->RemoveAllPipelines();
        }

        for (int i = 0; i < numScripts; i++) {
            vtkNew <vtkCPPythonScriptPipeline> pipeline;
            pipeline->Initialize(scripts[i]);
            Processor->AddPipeline(pipeline.GetPointer());
        }
    }

    void Finalize() {
        if (Processor) {
            Processor->Delete();
            Processor = NULL;
        }
    }

//    void CoProcess(Grid &grid, Attributes &attributes, double time,
//                   unsigned int timeStep, bool lastTimeStep) {
//        vtkNew <vtkCPDataDescription> dataDescription;
//        dataDescription->AddInput("input");
//        dataDescription->SetTimeData(time, timeStep);
//
//        if (lastTimeStep == true) {
//            dataDescription->ForceOutputOn();
//        }
//
//        if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
//            vtkCPInputDataDescription *idd =
//                    dataDescription->GetInputDescriptionByName("input");
//            BuildVTKDataStructures(grid, attributes, idd);
//            idd->SetGrid(vtkGrid);
//            Processor->CoProcess(dataDescription.GetPointer());
//        }
//    }

    void CoProcess(femocs::Femocs &project, double time, unsigned int timeStep, bool lastTimeStep) {
        vtkNew <vtkCPDataDescription> dataDescription;
        dataDescription->AddInput("input");
        dataDescription->SetTimeData(time, timeStep);

        if (lastTimeStep == true) {
            dataDescription->ForceOutputOn();
        }

        // Mesh data

        const double *nodes = NULL;
        const int n_coordinates = 3;
        const int n_nodes = project.export_data(&nodes, "nodes");
        // for (int i = 0; i < n_nodes; i++) {
        //   int I = n_coordinates * i;
        //   printf("%.3f %.3f %.3f\n", nodes[I], nodes[I + 1], nodes[I + 2]);
        // }

        const int *cells = NULL;
        const int n_nodes_per_cell = 4;
        const int n_cells = project.export_data(&cells, "quadrangles"); // export hexahedron
        // for (int i = 0; i < n_cells; i++) {
        //   int I = n_nodes_per_cell * i;
        //   printf("%d %d %d %d\n", cells[I], cells[I + 1], cells[I + 2],
        //          cells[I + 3]);
        // }

        // Field data

        double elfield_data[n_nodes] = {0};
        double temperature_data[n_nodes] = {0};
        project.export_data(elfield_data, n_nodes, "elfield_norm");
        project.export_data(temperature_data, n_nodes, "temperature");

        vtkSmartPointer <vtkDoubleArray> elfield =
                vtkSmartPointer<vtkDoubleArray>::New();
        elfield->SetName("elfield_norm");
        for (int i = 0; i < n_nodes; i++)
            elfield->InsertNextValue(elfield_data[i]);

        vtkSmartPointer <vtkDoubleArray> temperature =
                vtkSmartPointer<vtkDoubleArray>::New();
        temperature->SetName("temperature");
        for (int i = 0; i < n_nodes; i++)
            temperature->InsertNextValue(temperature_data[i]);

        vtkSmartPointer <vtkPointData> pointFieldData =
                vtkSmartPointer<vtkPointData>::New();
        pointFieldData->AddArray(elfield);
        pointFieldData->AddArray(temperature);

        // making VTK points and cells

        double coordinates[n_nodes][n_coordinates];
        for (int i = 0; i < n_nodes; i++) {
            for (int j = 0; j < n_coordinates; j++)
                coordinates[i][j] = nodes[n_coordinates * i + j];
        }
        vtkSmartPointer <vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for (int i = 0; i < n_nodes; i++)
            points->InsertPoint(i, coordinates[i]);

        vtkIdType cellsIndicies[n_cells][n_nodes_per_cell];
        for (int i = 0; i < n_cells; i++) {
            for (int j = 0; j < n_nodes_per_cell; j++)
                cellsIndicies[i][j] = cells[n_nodes_per_cell * i + j];
        }

        // updating VTK grid

        grid->Allocate(n_cells);
        for (int i = 0; i < n_cells; i++)
            grid->InsertNextCell(VTK_QUAD, n_nodes_per_cell, cellsIndicies[i]);
        grid->SetPoints(points);
        grid->SetFieldData(pointFieldData);

        // passing to Catalyst

        if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
            vtkCPInputDataDescription *idd = dataDescription->GetInputDescriptionByName("input");
            idd->SetGrid(grid);
            Processor->CoProcess(dataDescription.GetPointer());
        }

        // writing VTU output

//        vtkSmartPointer <vtkXMLUnstructuredGridWriter> writer =
//                vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
//        string filename = "vtk-output_iteration-" + to_string(timeStep) + ".vtu";
//        writer->SetFileName(filename.c_str());
//        writer->SetInputData(grid);
//        writer->Write();
    }
} // namespace Adaptor

//Grid::Grid() {}
//
//void Grid::Initialize(const unsigned int numPoints[3],
//                      const double spacing[3]) {
//    // ...
//}
//
//size_t Grid::GetNumberOfPoints() { return this->Points.size() / 3; }
//
//size_t Grid::GetNumberOfCells() { return this->Cells.size() / 8; }
//
//double *Grid::GetPointsArray() {
//    if (this->Points.empty()) {
//        return NULL;
//    }
//    return &(this->Points[0]);
//}
//
//double *Grid::GetPoint(size_t pointID) {
//    if (pointID >= this->Points.size()) {
//        return NULL;
//    }
//    return &(this->Points[pointID * 3]);
//}
//
//unsigned int *Grid::GetCellPoints(size_t cellID) {
//    if (cellID >= this->Cells.size()) {
//        return NULL;
//    }
//    return &(this->Cells[cellID * 8]);
//}
//
//Attributes::Attributes() { this->GridPtr = NULL; }
//
//void Attributes::Initialize(Grid *grid) { this->GridPtr = grid; }
//
//void Attributes::UpdateFields(double time) {
//    // ...
//}
//
//double *Attributes::GetVelocityArray() {
//    if (this->Velocity.empty()) {
//        return NULL;
//    }
//    return &this->Velocity[0];
//}
//
//float *Attributes::GetPressureArray() {
//    if (this->Pressure.empty()) {
//        return NULL;
//    }
//    return &this->Pressure[0];
//}
