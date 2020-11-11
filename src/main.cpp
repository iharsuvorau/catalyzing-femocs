#include "AtomReader.h"
#include "Config.h"
#include "Femocs.h"
#include "GeneralProject.h"
#include "Macros.h"
#include "TetgenMesh.h"
#include <stdlib.h>
#include <vector>
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

void read_ckx(const string &file_name, double *x, double *y, double *z) {
  ifstream in_file(file_name, ios::in);
  require(in_file.is_open(), "Did not find a file " + file_name);

  int n_atoms, type, id;
  string line;
  istringstream iss;

  getline(in_file, line); // Read number of atoms
  iss.clear();
  iss.str(line);
  iss >> n_atoms;

  getline(in_file, line); // Skip comments line

  id = -1;
  // keep storing values from the text file as long as data exists:
  while (++id < n_atoms && getline(in_file, line)) {
    iss.clear();
    iss.str(line);
    iss >> type >> x[id] >> y[id] >> z[id];
  }
}

void read_xyz(const string &file_name, double *x, double *y, double *z) {
  ifstream in_file(file_name, ios::in);
  require(in_file.is_open(), "Did not find a file " + file_name);

  int n_atoms, type, id;
  string elem, line;
  istringstream iss;

  getline(in_file, line); // Read number of atoms
  iss.clear();
  iss.str(line);
  iss >> n_atoms;

  getline(in_file, line); // Skip comments line

  id = -1;
  // keep storing values from the text file as long as data exists:
  while (++id < n_atoms && getline(in_file, line)) {
    iss.clear();
    iss.str(line);
    iss >> elem >> x[id] >> y[id] >> z[id] >> type;
  }
}

void read_atoms(const string &file_name, double *x, double *y, double *z) {
  string file_type = get_file_type(file_name);

  if (file_type == "xyz")
    read_xyz(file_name, x, y, z);
  else if (file_type == "ckx")
    read_ckx(file_name, x, y, z);
  else
    require(false, "Unsupported file type: " + file_type);
}

void read_n_atoms(const string &file_name, int &n_atoms) {
  ifstream in_file(file_name, ios::in);
  require(in_file.is_open(), "Did not find a file " + file_name);

  string line;
  istringstream iss;

  getline(in_file, line); // Read number of atoms
  iss.clear();
  iss.str(line);
  iss >> n_atoms;
}

void export_vtk_unstructured_grid(TetgenMesh *mesh, vtkUnstructuredGrid *grid) {
  const size_t n_nodes = mesh->hexs.get_n_nodes();
  const size_t n_cells = mesh->hexs.size();
  const size_t celltype = mesh->hexs.get_cell_type();

  // vtkUnstructuredGrid *grid = new vtkUnstructuredGrid::New();
  vtkPoints *points = vtkPoints::New();

  // size to allocate: numCells*(numPointsPerCell+1)
  grid->Allocate(n_cells * (celltype + 1));

  points->SetNumberOfPoints(n_nodes);
  for (size_t node = 0; node < n_nodes; node++) {
    auto n = mesh->hexs.get_node(node);
    points->SetPoint(node, n.x, n.y, n.z);
  }
  grid->SetPoints(points);
  points->Delete();

  for (size_t cl = 0; cl < n_cells; cl++) {
    auto cell = mesh->hexs.get_cell(cl);
    int cell_size = cell.size();
    vtkIdType ids[cell_size];
    for (int i = 0; i < cell_size; i++) {
      ids[i] = cell[i];
    }
    grid->InsertNextCell(celltype, cell_size, ids);
  }

  // TODO: add field data
  // TODO: is there a way to generalize the function more?
}

int main() {
  string filename = "in/md.in";

  // run(filename);

  femocs::Femocs project(filename);

  string cmd1 = "infile";
  string infile = "";
  int success = 0;

  success = project.parse_command(cmd1, infile);
  success = 0;
  print_progress("\n> reading " + cmd1, infile != "");

  int n_atoms = 0;

  if (infile != "")
    read_n_atoms(infile, n_atoms);

  double *x = (double *)malloc(n_atoms * sizeof(double));
  double *y = (double *)malloc(n_atoms * sizeof(double));
  double *z = (double *)malloc(n_atoms * sizeof(double));

  if (infile != "")
    read_atoms(infile, x, y, z);

  cout << "printing out x:" << endl;
  for (int i = 0; i < n_atoms; i++) {
    cout << "\t" << x[i] << endl;
  }

  int n_iterations = 1;
  bool add_rnd_noise = false;

  double output_data[n_atoms * 3] = {0};

  vector<string> labels = {
      // "vec",
      //                    "vec_norm",
      //                    "scalar",
      //  "potential",
      // "temperature",
      //  "rho",
      //  "rho_norm",
      //  "kin_energy",
      //  "pot_energy",
      //  "pair_potential",
      //  "parcas_force",
      //  "charge_force",
      //  "force_norm",
      //  "charge",
      //  "charge_density",
      //  "parcas_velocity",
      //  "velocity_norm",
      //  "area",
      //  "atom_type"

      // scalar
      //  "heat",
      // "elfield_norm",

      // 3D vector
      // "elfield",
      // "force",
      "velocity",
  };

  TetgenMesh *mesh;

  for (int i = 1; i <= n_iterations; ++i) {
    if (n_iterations > 1)
      cout << "\n> iteration " << i << endl;

    success = project.import_atoms(infile, add_rnd_noise);
    success += project.run();

    cout << "exporting a mesh" << endl;
    mesh = project.export_mesh();
    cout << "Mesh statistics from my program:" << endl
         << "#hexs=" << mesh->hexs.size() << ", #tets=" << mesh->tets.size()
         << ", #quads=" << mesh->quads.size() << ", #tris=" << mesh->tris.size()
         << ", #edges=" << mesh->edges.size()
         << ", #nodes=" << mesh->nodes.size() << endl;

    cout << "node: " << mesh->nodes[0] << endl;
    cout << "tet: " << mesh->tets[0] << endl;
    cout << "quad: " << mesh->quads[0] << endl;
    cout << "edge: " << mesh->edges[0] << endl;
    cout << "hexs: " << mesh->hexs[0] << endl;

    mesh->hexs.write("hexs_output.vtk");

    vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();
    export_vtk_unstructured_grid(mesh, grid);
    grid->PrintSelf(cout, vtkIndent(2));
    grid->Delete();

    for (auto label : labels) {
      cout << "exporting " << label << endl;
      project.export_data(output_data, n_atoms, label);
      for (int j = 0; j < n_atoms * 3; j++) {
        cout << "\t" << label << "[" << j << "] : " << output_data[j] << endl;
      }
    }
  }

  print_progress("\n> full run of Femocs", success == 0);

  free(x);
  free(y);
  free(z);

  return 0;
}