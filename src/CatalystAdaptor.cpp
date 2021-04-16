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
#include "CatalystAdaptor.h"
#include "Femocs.h"

namespace CatalystAdaptor {
    vtkCPProcessor *processor = NULL;

    // FEM-related globals
    const char *meshCellType;
    VTKCellType vtkMeshCellType;
    int numNodesPerCell;

    // MD-related globals
    double **atoms;
    int numAtoms;

    void Initialize(const char *path, const char *cellType) {
        std::cout << "CatalystAdaptor::Initialize has been called" << std::endl;

        if (processor == NULL) {
            processor = vtkCPProcessor::New();
            processor->Initialize();
        } else {
            processor->RemoveAllPipelines();
        }

        // registering python visualization script
        vtkNew <vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(path);
        processor->AddPipeline(pipeline.GetPointer());

        // parsing mesh-related arguments
        meshCellType = cellType;
        // https://vtk.org/doc/nightly/html/vtkCellType_8h.html
        if (strcmp(meshCellType, "quadrangles") == 0) {
            numNodesPerCell = 4;
            vtkMeshCellType = VTK_QUAD;
        } else if (strcmp(meshCellType, "hexahedra") == 0) {
            numNodesPerCell = 8;
            vtkMeshCellType = VTK_HEXAHEDRON;
        } else {
            std::cout << "cell type is undefined" << std::endl;
            exit(1);
        }
    }

    void Finalize() {
        std::cout << "CatalystAdaptor::Finalize has been called" << std::endl;
        if (processor) {
            processor->Delete();
            processor = NULL;
        }
    }

    void ImportAtoms(double **arr, int num) {
        std::cout << "CatalystAdaptor::ImportAtoms has been called" << std::endl;
        atoms = arr;
        numAtoms = num;
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
        std::cout << "CatalystAdaptor::CoProcess has been called" << std::endl;


        // Generating FEM grid using mesh data

        // preparing mesh points
        const double *nodes = NULL;
        const int numNodes = project.export_data(&nodes, "nodes");
        vtkNew <vtkPoints> meshPoints;
        meshPoints->Allocate(numNodes);
        for (int i = 0; i < numNodes; i++) {
            meshPoints->InsertNextPoint(nodes[i * 3], nodes[i * 3 + 1], nodes[i * 3 + 2]);
        }

        // preparing mesh cells
        const int *cells = NULL;
        const int numCells = project.export_data(&cells, meshCellType);
        vtkIdType meshCellIDs[numCells][numNodesPerCell];
        for (int i = 0; i < numCells; i++) {
            for (int j = 0; j < numNodesPerCell; j++) {
                meshCellIDs[i][j] = cells[numNodesPerCell * i + j];
            }
        }

        // making mesh grid
        vtkNew <vtkUnstructuredGrid> meshGrid;
        meshGrid->SetPoints(meshPoints);
        for (int i = 0; i < numCells; i++)
            meshGrid->InsertNextCell(vtkMeshCellType, numNodesPerCell, meshCellIDs[i]);


        // Generating MD grid using atomistic data

        // preparing atoms' points and cells
        vtkNew <vtkPoints> atomPoints;
        atomPoints->Allocate(numAtoms);
        vtkNew <vtkCellArray> atomCells;
        for (int i = 0; i < numAtoms; i++) {
            atomPoints->InsertPoint(i, (*atoms)[i * 3], (*atoms)[i * 3 + 1], (*atoms)[i * 3 + 2]);
            atomCells->InsertNextCell(1);
            atomCells->InsertCellPoint(i);
        }

        // extracting field data for atoms

        double temperatureData[numAtoms] = {0};
        project.export_data(temperatureData, numAtoms, "temperature");
        vtkNew <vtkDoubleArray> temperature;
        temperature->SetName("temperature");
        temperature->SetArray(temperatureData, numAtoms, 1);

        double elfieldNormData[numAtoms] = {0};
        project.export_data(elfieldNormData, numAtoms, "elfield_norm");
        vtkNew <vtkDoubleArray> elfieldNorm;
        elfieldNorm->SetName("electric field normals");
        elfieldNorm->SetArray(elfieldNormData, numAtoms, 1);

        double elfieldData[numAtoms * 3] = {0};
        project.export_data(elfieldData, numAtoms * 3, "elfield");
        vtkNew <vtkDoubleArray> elfield;
        elfield->SetName("electic field");
        elfield->SetNumberOfComponents(3);
        elfield->SetNumberOfTuples(numAtoms);
        for (int i = 0; i < numAtoms * 3; i++) {
            elfieldData[i] *= -1;
        }
        elfield->SetArray(elfieldData, numAtoms * 3, 1);

        // making atomistic grid
        vtkNew <vtkUnstructuredGrid> atomsGrid;
        atomsGrid->SetPoints(atomPoints);
        atomsGrid->SetCells(VTK_VERTEX, atomCells);
        atomsGrid->GetPointData()->AddArray(temperature);
        atomsGrid->GetPointData()->AddArray(elfieldNorm);
        atomsGrid->GetPointData()->SetVectors(elfield);


        // Passing data to Catalyst

        dataDescription->GetInputDescriptionByName("MD")->SetGrid(atomsGrid);
        dataDescription->GetInputDescriptionByName("FEM")->SetGrid(meshGrid);
        processor->CoProcess(dataDescription.GetPointer());
    }
}
