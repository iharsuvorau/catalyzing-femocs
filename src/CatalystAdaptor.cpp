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
#include "CatalystAdaptor.h"
#include "Femocs.h"

namespace CatalystAdaptor {
    vtkCPProcessor *Processor = NULL;
    vtkSmartPointer <vtkUnstructuredGrid> grid;

    void Initialize(const char *path) {
        printf("CatalystAdaptor::Initialize has been called\n");

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

        vtkNew <vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(path);
        Processor->AddPipeline(pipeline.GetPointer());
    }

    void Finalize() {
        printf("CatalystAdaptor::Finalize has been called\n");
        if (Processor) {
            Processor->Delete();
            Processor = NULL;
        }
    }

    void CoProcess(femocs::Femocs &project, const double time, const unsigned int timeStep, const bool lastTimeStep) {
        printf("CatalystAdaptor::CoProcess has been called\n");

        // Mesh data

        const double *nodes = NULL;
        const int n_coordinates = 3;
        const int n_nodes = project.export_data(&nodes, "nodes");
        printf("nodes: %.3d\n", n_nodes);
        // for (int i = 0; i < n_nodes; i++) {
        //   int I = n_coordinates * i;
        //   printf("%.3f %.3f %.3f\n", nodes[I], nodes[I + 1], nodes[I + 2]);
        // }

        const int *cells = NULL;
        const int n_nodes_per_cell = 4;
        const int n_cells = project.export_data(&cells, "quadrangles"); // export hexahedron
        printf("cells: %.3d\n", n_cells);
        // for (int i = 0; i < n_cells; i++) {
        //   int I = n_nodes_per_cell * i;
        //   printf("%d %d %d %d\n", cells[I], cells[I + 1], cells[I + 2],
        //          cells[I + 3]);
        // }

        // Field data

        double elfield_data[n_nodes] = {0};
        double temperature_data[n_nodes] = {0};
        project.export_data(elfield_data, n_nodes, "elfield");
        project.export_data(temperature_data, n_nodes, "temperature");

        vtkSmartPointer <vtkDoubleArray> elfield =
                vtkSmartPointer<vtkDoubleArray>::New();
        elfield->SetName("elfield");
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
        grid->Squeeze();

        // passing data to Catalyst

        vtkSmartPointer <vtkCPDataDescription> dataDescription = vtkSmartPointer<vtkCPDataDescription>::New();
        dataDescription->AddInput("input");
        dataDescription->SetTimeData(time, timeStep);

        if (lastTimeStep) {
            dataDescription->ForceOutputOn();
        }

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
}
