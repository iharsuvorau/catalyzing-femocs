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
#include <vtkUnstructuredGrid.h>

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

  for (int i = 1; i <= n_iterations; ++i) {
    if (n_iterations > 1)
      cout << "\n> iteration " << i << endl;

    success = project.import_atoms("", add_rnd_noise);
    success += project.run();

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

    // TetgenMesh *mesh;
    // vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();
    // export_vtk_unstructured_grid(mesh, grid);
    // grid->PrintSelf(cout, vtkIndent(2));
    // grid->Delete();

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
