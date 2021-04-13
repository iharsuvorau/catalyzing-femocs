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
        dataDescription->AddInput("MD");
        dataDescription->AddInput("FEM");
        dataDescription->SetTimeData(time, timeStep);
        if (lastTimeStep) {
            dataDescription->ForceOutputOn();
        }

        // determining if co-processing should be done
        if (processor->RequestDataDescription(dataDescription.GetPointer()) == 0) {
            return;
        }
        printf("CatalystAdaptor::CoProcess has been called\n");


        // Generating FEM grid using mesh data

        // preparing mesh points
        const double *nodes = NULL;
        const int n_nodes = project.export_data(&nodes, "nodes");
        std::cout << "n_nodes = " << n_nodes << ", natoms_g = " << natoms_g << std::endl;
        vtkNew <vtkPoints> mesh_points;
        mesh_points->Allocate(n_nodes);
        for (int i = 0; i < n_nodes; i++) {
            mesh_points->InsertNextPoint(nodes[i * 3], nodes[i * 3 + 1], nodes[i * 3 + 2]);
        }

        // preparing mesh cells
        const int *cells = NULL;
        const int n_cells = project.export_data(&cells, cell_type_g);
        vtkIdType mesh_cell_ids[n_cells][n_nodes_per_cell];
        for (int i = 0; i < n_cells; i++) {
            for (int j = 0; j < n_nodes_per_cell; j++) {
                mesh_cell_ids[i][j] = cells[n_nodes_per_cell * i + j];
            }
        }

        // making mesh grid
        vtkNew <vtkUnstructuredGrid> mesh_grid;
        mesh_grid->SetPoints(mesh_points);
        for (int i = 0; i < n_cells; i++)
            mesh_grid->InsertNextCell(vtk_cell_type, n_nodes_per_cell, mesh_cell_ids[i]);


        // Generating MD grid using atomistic data

        // preparing atoms' points and cells
        vtkNew <vtkPoints> atom_points;
        atom_points->Allocate(natoms_g);
        vtkNew <vtkCellArray> atom_cells;
        for (int i = 0; i < natoms_g; i++) {
            atom_points->InsertPoint(i, (*atoms_g)[i * 3], (*atoms_g)[i * 3 + 1], (*atoms_g)[i * 3 + 2]);
            atom_cells->InsertNextCell(1);
            atom_cells->InsertCellPoint(i);
        }

        // extracting field data for atoms
        double temperature_data[n_nodes] = {0};
        project.export_data(temperature_data, n_nodes, "temperature"); // NOTE: this field data is for atoms, not mesh nodes
        vtkNew <vtkDoubleArray> temperature;
        temperature->SetName("temperature");
        temperature->SetArray(temperature_data, n_nodes, 1);

        // making atomistic grid
        vtkNew <vtkUnstructuredGrid> atoms_grid;
        atoms_grid->SetPoints(atom_points);
        atoms_grid->SetCells(VTK_VERTEX, atom_cells);
        atoms_grid->GetPointData()->AddArray(temperature);


        // passing data to Catalyst
        dataDescription->GetInputDescriptionByName("MD")->SetGrid(atoms_grid);
        dataDescription->GetInputDescriptionByName("FEM")->SetGrid(mesh_grid);
        processor->CoProcess(dataDescription.GetPointer());
    }
}
