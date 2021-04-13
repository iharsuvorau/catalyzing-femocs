#include <cstdlib>
#include <iterator>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertex.h>
#include <vtkMergeCells.h>
#include <vtkMergeFields.h>
#include "CatalystAdaptor.h"
#include "Femocs.h"

namespace CatalystAdaptor {
    vtkCPProcessor *processor = NULL;
    const char *cell_type_g;
    const char *field_label_g;
    int n_coordinates;
    int n_nodes_per_cell;
    VTKCellType vtk_cell_type;
    double **atoms_g;
    int natoms_g;

    void Initialize(const char *path, const char *cell_type, const char *field_label) {
        printf("CatalystAdaptor::Initialize has been called\n");

        if (processor == NULL) {
            processor = vtkCPProcessor::New();
            processor->Initialize();
        } else {
            processor->RemoveAllPipelines();
        }

        // path to the python visualization script
        vtkNew <vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(path);
        processor->AddPipeline(pipeline.GetPointer());

        // parsing arguments

        cell_type_g = cell_type; // TODO: make enum
        field_label_g = field_label; // TODO: make enum

        n_coordinates = 3;

        // https://vtk.org/doc/nightly/html/vtkCellType_8h.html
        if (strcmp(cell_type_g, "quadrangles") == 0) {
            n_nodes_per_cell = 4;
            vtk_cell_type = VTK_QUAD;
        } else if (strcmp(cell_type_g, "hexahedra") == 0) {
            n_nodes_per_cell = 8;
            vtk_cell_type = VTK_HEXAHEDRON;
        } else {
            std::cout << "cell type is undefined" << std::endl;
            exit(1);
        }
    }

    void Finalize() {
        printf("CatalystAdaptor::Finalize has been called\n");
        if (processor) {
            processor->Delete();
            processor = NULL;
        }
    }

    void ImportAtoms(double **atoms, int natoms) {
        atoms_g = atoms;
        natoms_g = natoms;
        std::cout << "CatalystAdaptor::ImportAtoms has been called" << std::endl
                  << "\tthere are " << natoms << " atoms" << std::endl
                  << "\tatoms_g[0] = " << *atoms_g[0] << std::endl
                  << "\tatoms_g[1] = " << *atoms_g[1] << std::endl
                  << "\tatoms_g[2] = " << *atoms_g[2] << std::endl
                  << "\tatoms_g[3] = " << *atoms_g[3] << std::endl
                  << "\tatoms_g[4] = " << *atoms_g[4] << std::endl
                  << "\tatoms_g[5] = " << *atoms_g[5] << std::endl;
    }

    void CoProcess(femocs::Femocs &project, const double time, const unsigned int timeStep, const bool lastTimeStep) {
        // creating data description
        vtkNew <vtkCPDataDescription> dataDescription;
        dataDescription->AddInput("input");
        dataDescription->SetTimeData(time, timeStep);
        if (lastTimeStep) {
            dataDescription->ForceOutputOn();
        }

        // determining if co-processing should be done
        if (processor->RequestDataDescription(dataDescription.GetPointer()) == 0) {
            return;
        }
        printf("CatalystAdaptor::CoProcess has been called\n");

        const double *nodes = NULL;
        const int n_nodes = project.export_data(&nodes, "nodes"); // TODO: rename n_nodes, it's n_atoms
        printf("nodes: %.3d\n", n_nodes);

//        const int *cells = NULL;
//        const int n_cells = project.export_data(&cells, cell_type_g);
//        printf("%s cells: %.3d\n", cell_type_g, n_cells);
//
//        // Field data
//
//        double temperature_data[n_nodes] = {0};
//        project.export_data(temperature_data, n_nodes, "temperature");
//        vtkNew <vtkDoubleArray> temperature;
//        temperature->SetName("temperature");
//        temperature->SetArray(temperature_data, n_nodes, 1);
//
//        // making VTK points and cells
//
//        double coordinates[n_nodes][n_coordinates];
//        for (int i = 0; i < n_nodes; i++) {
//            for (int j = 0; j < n_coordinates; j++)
//                coordinates[i][j] = nodes[n_coordinates * i + j];
//        }
//        vtkNew <vtkPoints> points;
//        for (int i = 0; i < n_nodes; i++)
//            points->InsertPoint(i, coordinates[i]);
//
//        vtkIdType cellsIndicies[n_cells][n_nodes_per_cell];
//        for (int i = 0; i < n_cells; i++) {
//            for (int j = 0; j < n_nodes_per_cell; j++) {
//                cellsIndicies[i][j] = cells[n_nodes_per_cell * i + j];
//            }
//        }
//
//        // creating VTK grid
//        vtkNew <vtkUnstructuredGrid> grid;
//        grid->Allocate(n_cells);
//        for (int i = 0; i < n_cells; i++)
//            grid->InsertNextCell(vtk_cell_type, n_nodes_per_cell, cellsIndicies[i]);
//        grid->SetPoints(points);
//        grid->GetPointData()->AddArray(temperature);
////            grid->GetCellData()->AddArray(temperature);
////            grid->GetCellData()->AddArray(elfield_norm);
//        grid->Squeeze();






//        // grid with mesh
//
//        const int *cells = NULL;
//        const int n_cells = project.export_data(&cells, cell_type_g);
//        printf("%s cells: %.3d\n", cell_type_g, n_cells);
//
//        double coordinates[n_nodes][n_coordinates];
//        for (int i = 0; i < n_nodes; i++) {
//            for (int j = 0; j < n_coordinates; j++)
//                coordinates[i][j] = nodes[n_coordinates * i + j];
//        }
//        vtkNew <vtkPoints> points;
//        for (int i = 0; i < n_nodes; i++)
//            points->InsertPoint(i, coordinates[i]);
//
//        vtkIdType cellsIndicies[n_cells][n_nodes_per_cell];
//        for (int i = 0; i < n_cells; i++) {
//            for (int j = 0; j < n_nodes_per_cell; j++) {
//                cellsIndicies[i][j] = cells[n_nodes_per_cell * i + j];
//            }
//        }
//
//        vtkNew <vtkUnstructuredGrid> mesh_grid;
//        mesh_grid->Allocate(n_cells);
//        for (int i = 0; i < n_cells; i++)
//            mesh_grid->InsertNextCell(vtk_cell_type, n_nodes_per_cell, cellsIndicies[i]);
//        mesh_grid->SetPoints(points);
//
//        // grid with atoms
//
//        vtkNew <vtkPoints> atom_points;
//        vtkNew <vtkCellArray> vertices;
//
//        for (int i = 0; i < natoms_g; i++) {
//            atom_points->InsertPoint(i, (*atoms_g)[i * 3], (*atoms_g)[i * 3 + 1], (*atoms_g)[i * 3 + 2]);
//
//            vertices->InsertNextCell(1);
//            vertices->InsertCellPoint(i);
//        }
//
//        double temperature_data[n_nodes] = {0};
//        project.export_data(temperature_data, n_nodes, "temperature");
//        vtkNew <vtkDoubleArray> temperature;
//        temperature->SetName("temperature");
//        temperature->SetArray(temperature_data, n_nodes, 1);
//
//        vtkNew <vtkUnstructuredGrid> atom_grid;
//        atom_grid->SetPoints(atom_points);
//        atom_grid->SetCells(VTK_VERTEX, vertices);
//        atom_grid->GetPointData()->AddArray(temperature);
//
//        // merge
//
//        vtkNew <vtkMergeCells> merge;
//        vtkNew <vtkUnstructuredGrid> merged_grid;
//
//        merge->SetUnstructuredGrid(merged_grid);
//        merge->SetTotalNumberOfCells(mesh_grid->GetNumberOfCells() + atom_grid->GetNumberOfCells());
//        merge->SetTotalNumberOfPoints(mesh_grid->GetNumberOfPoints() + atom_grid->GetNumberOfPoints());
//        merge->SetTotalNumberOfDataSets(2);
//        merge->MergeDataSet(mesh_grid);
//        merge->MergeDataSet(atom_grid);



        // grid with atoms

        vtkNew <vtkPoints> points;
        vtkNew <vtkCellArray> vertices;

        points->Allocate(natoms_g + n_nodes);

        for (int i = 0; i < natoms_g; i++) {
            points->InsertPoint(i, (*atoms_g)[i * 3], (*atoms_g)[i * 3 + 1], (*atoms_g)[i * 3 + 2]);

            vertices->InsertNextCell(1);
            vertices->InsertCellPoint(i);
        }

        double temperature_data[n_nodes] = {0};
        project.export_data(temperature_data, n_nodes, "temperature");
        double temp_data[n_nodes + natoms_g];
        for (int i = 0; i < n_nodes + natoms_g; i++) {
            if (i < n_nodes) {
                temp_data[i] = temperature_data[i];
            } else {
                temp_data[i] = 0;
            }
        }
        vtkNew <vtkDoubleArray> temperature;
        temperature->SetName("temperature");
        temperature->SetArray(temp_data, n_nodes + natoms_g, 1);

        // mesh

        const int *cells = NULL;
        const int n_cells = project.export_data(&cells, cell_type_g);

        vtkIdType mesh_point_ids[n_nodes];
        for (int i = 0; i < n_nodes; i++) {
            mesh_point_ids[i] = points->InsertNextPoint(nodes[i * 3], nodes[i * 3 + 1], nodes[i * 3 + 2]);
        }

        vtkIdType mesh_cell_ids[n_cells][n_nodes_per_cell];
        for (int i = 0; i < n_cells; i++) {
            for (int j = 0; j < n_nodes_per_cell; j++) {
                mesh_cell_ids[i][j] = cells[n_nodes_per_cell * i + j] + natoms_g;
            }
        }


//        double coordinates[n_nodes][n_coordinates];
//        for (int i = natoms_g; i < natoms_g + n_nodes; i++) {
//            for (int j = 0; j < n_coordinates; j++)
//                coordinates[i][j] = nodes[n_coordinates * i + j];
//        }
//        for (int i = natoms_g; i < natoms_g + n_nodes; i++)
//            points->InsertPoint(i, coordinates[i]);
//
//        vtkIdType cellsIndicies[n_cells][n_nodes_per_cell];
//        for (int i = natoms_g; i < natoms_g + n_cells; i++) {
//            for (int j = 0; j < n_nodes_per_cell; j++) {
//                cellsIndicies[i][j] = cells[n_nodes_per_cell * i + j];
//            }
//        }

        // grid

        vtkNew <vtkUnstructuredGrid> grid;
        grid->SetPoints(points);
        grid->SetCells(VTK_VERTEX, vertices);
        for (int i = 0; i < n_cells; i++)
            grid->InsertNextCell(vtk_cell_type, n_nodes_per_cell, mesh_cell_ids[i]);
        grid->GetPointData()->AddArray(temperature);


        // passing data to Catalyst
        dataDescription->GetInputDescriptionByName("input")->SetGrid(grid);
        processor->CoProcess(dataDescription.GetPointer());
    }
}
