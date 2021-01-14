#include "Femocs.h"
#include "GeneralProject.h"
#include "Globals.h"
#include "Macros.h"
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkMeta.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace std;
using namespace femocs;

void print_progress(const string &message, const bool contition) {
  cout << message << ":  ";
  if (contition)
    cout << "passed" << endl;
  else
    cout << "failed" << endl;
}

void print_field(const string label, const int n_nodes, double *field_data) {
  for (int i = 0; i < n_nodes; i++) {
    if (field_data[i] == 0)
      continue;
    cout << label << "[" << i << "] : " << field_data[i] << endl;
  }
}

int main() {
  string filename = "in/md.in";

  femocs::Femocs project(filename);

  int success = 0;
  int n_iterations = 100;
  bool add_rnd_noise = true;

  for (int iter_i = 1; iter_i <= n_iterations; ++iter_i) {
    if (n_iterations > 1)
      cout << "\n> iteration " << iter_i << endl;

    success = project.import_atoms("", add_rnd_noise);
    success += project.run();

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

    vtkSmartPointer<vtkDoubleArray> elfield =
        vtkSmartPointer<vtkDoubleArray>::New();
    elfield->SetName("elfield");
    for (int i = 0; i < n_nodes; i++)
      elfield->InsertNextValue(elfield_data[i]);

    vtkSmartPointer<vtkDoubleArray> temperature =
        vtkSmartPointer<vtkDoubleArray>::New();
    temperature->SetName("temperature");
    for (int i = 0; i < n_nodes; i++)
      temperature->InsertNextValue(temperature_data[i]);

    vtkSmartPointer<vtkPointData> pointFieldData =
        vtkSmartPointer<vtkPointData>::New();
    pointFieldData->AddArray(elfield);
    pointFieldData->AddArray(temperature);

    // VTK mesh

    double coordinates[n_nodes][n_coordinates];
    for (int i = 0; i < n_nodes; i++) {
      for (int j = 0; j < n_coordinates; j++)
        coordinates[i][j] = nodes[n_coordinates * i + j];
    }
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < n_nodes; i++)
      points->InsertPoint(i, coordinates[i]);

    vtkIdType cellsIndicies[n_cells][n_nodes_per_cell];
    for (int i = 0; i < n_cells; i++) {
      for (int j = 0; j < n_nodes_per_cell; j++)
        cellsIndicies[i][j] = cells[n_nodes_per_cell * i + j];
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->Allocate(n_cells);
    for (int i = 0; i < n_cells; i++)
      grid->InsertNextCell(VTK_QUAD, n_nodes_per_cell, cellsIndicies[i]);
    grid->SetPoints(points);
    grid->SetFieldData(pointFieldData);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    string filename = "vtk-output_iteration-" + to_string(iter_i) + ".vtu";
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->Write();
  }

  print_progress("\n> full run of Femocs", success == 0);
  return 0;
}
