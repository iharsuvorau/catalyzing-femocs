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
    vtkCPProcessor *processor = NULL;
    vtkSmartPointer <vtkUnstructuredGrid> grid;
    char *cell_type_g;
    char *field_label_g;

    void Initialize(const char *path, const char *cell_type, const char *field_label) {
        printf("CatalystAdaptor::Initialize has been called\n");

        if (grid == NULL) {
            grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        } else {
            grid->Initialize();
        }

        if (processor == NULL) {
            processor = vtkCPProcessor::New();
            processor->Initialize();
        } else {
            processor->RemoveAllPipelines();
        }

        vtkNew <vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(path);
        processor->AddPipeline(pipeline.GetPointer());

        cell_type_g = cell_type;
        field_label_g = field_label;
    }

    void Finalize() {
        printf("CatalystAdaptor::Finalize has been called\n");
        if (processor) {
            processor->Delete();
            processor = NULL;
        }
    }

    void CoProcess(femocs::Femocs &project, const double time, const unsigned int timeStep, const bool lastTimeStep) {
        printf("CatalystAdaptor::CoProcess has been called\n");

        // Mesh data

        const double *nodes = NULL;
        const int n_coordinates = 3;
        const int n_nodes = project.export_data(&nodes, "nodes"); // TODO: make a parameter
        printf("nodes: %.3d\n", n_nodes);

        const int *cells = NULL;
        const int n_nodes_per_cell = 4;
        const int n_cells = project.export_data(&cells, cell_type_g); // export hexahedron, quadrangles
        printf("cells: %.3d\n", n_cells);

        // Field data

//        double elfield_data[n_nodes] = {0};
//        vtkSmartPointer <vtkDoubleArray> elfield =
//                vtkSmartPointer<vtkDoubleArray>::New();
//        elfield->SetName("elfield");
//        for (int i = 0; i < n_nodes; i++)
//            elfield->InsertNextValue(elfield_data[i]);

        double field_data[n_nodes] = {0};
        project.export_data(field_data, n_nodes, field_label_g);

        vtkSmartPointer <vtkDoubleArray> field_data =
                vtkSmartPointer<vtkDoubleArray>::New();
        field_data->SetName(field_label_g);
        for (int i = 0; i < n_nodes; i++)
            field_data->InsertNextValue(field_data[i]);

        vtkSmartPointer <vtkPointData> pointFieldData =
                vtkSmartPointer<vtkPointData>::New();
//        pointFieldData->AddArray(elfield);
        pointFieldData->AddArray(field_data);

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

        if (processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
            vtkCPInputDataDescription *idd = dataDescription->GetInputDescriptionByName("input");
            idd->SetGrid(grid);
            processor->CoProcess(dataDescription.GetPointer());
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
