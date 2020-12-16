#include "AtomReader.h"
#include "Config.h"
#include "Femocs.h"
#include "GeneralProject.h"
#include "Macros.h"
#include "TetgenMesh.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vtkCellType.h>
#include <vtkMeta.h>
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

int main() {
  string filename = "in/md.in";

  femocs::Femocs project(filename);

  // double output_data[n_atoms * 3] = {0};
  int success = 0;
  int n_iterations = 1;
  bool add_rnd_noise = false;

  for (int iter_i = 1; iter_i <= n_iterations; ++iter_i) {
    if (n_iterations > 1)
      cout << "\n> iteration " << iter_i << endl;

    success = project.import_atoms("", add_rnd_noise);
    success += project.run();

    // Mesh data

    const double *nodes = NULL;
    const int n_coordinates = 3;
    const int n_nodes = project.export_data(&nodes, "nodes");
    for (int i = 0; i < n_nodes; i++) {
      int I = n_coordinates * i;
      printf("%.3f %.3f %.3f\n", nodes[I], nodes[I + 1], nodes[I + 2]);
    }

    const int *cells = NULL;
    const int n_nodes_per_cell = 4;
    const int n_cells = project.export_data(&cells, "quadrangles");
    for (int i = 0; i < n_cells; i++) {
      int I = n_nodes_per_cell * i;
      printf("%d %d %d %d\n", cells[I], cells[I + 1], cells[I + 2],
             cells[I + 3]);
    }

    // VTK mesh

    double coordinates[n_nodes][n_coordinates];
    for (int i = 0; i < n_nodes; i++) {
      for (int j = 0; j < n_coordinates; j++)
        coordinates[i][j] = nodes[n_coordinates * i + j];
    }
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < n_nodes; i++)
      points->InsertPoint(i, coordinates[i]);

    vtkIdType cellsIndecies[n_cells][n_nodes_per_cell];
    for (int i = 0; i < n_cells; i++) {
      for (int j = 0; j < n_nodes_per_cell; j++)
        cellsIndecies[i][j] = cells[n_nodes_per_cell * i + j];
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->Allocate(n_cells);
    for (int i = 0; i < n_cells; i++)
      grid->InsertNextCell(VTK_QUAD, n_nodes_per_cell, cellsIndecies[i]);
    grid->SetPoints(points);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    string filename = "vtk-output_iteration-" + to_string(iter_i) + ".vtu";
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->Write();

    // Field data

    // for (auto label : labels) {
    //   cout << "exporting " << label << endl;
    //   project.export_data(output_data, n_atoms, label);
    //   for (int j = 0; j < n_atoms * 3; j++) {
    //     cout << "\t" << label << "[" << j << "] : " << output_data[j] <<
    //     endl;
    //   }
    // }
  }

  print_progress("\n> full run of Femocs", success == 0);

  return 0;
}
