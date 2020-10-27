#include "Config.h"
#include "Femocs.h"
#include "Macros.h"
#include <stdlib.h>
#include <vector>

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

void config(string filename) {
  static bool first_call = true;
    bool fail;

  Config config;
  config.read_all(filename);
  MODES.MUTE = false;
  MODES.VERBOSE = true;

  // Pick correct flags for writing log file
  MODES.WRITELOG = conf.behaviour.n_write_log != 0;
  MODES.SHORTLOG = conf.behaviour.n_write_log < 0;


}

int main() {
  string filename = "in/md.in";
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

  ofstream output_mesh("output_mesh.vtk");

  for (int i = 1; i <= n_iterations; ++i) {
    if (n_iterations > 1)
      cout << "\n> iteration " << i << endl;

    success = project.import_atoms(infile, add_rnd_noise);
    success += project.run();

    for (auto label : labels) {
      cout << "exporting " << label << endl;
      project.export_data(output_data, n_atoms, label);
      for (int j = 0; j < n_atoms * 3; j++) {
        cout << "\t" << label << "[" << j << "] : " << output_data[j] << endl;
      }
    }

    project.export_mesh(output_mesh);
  }

  print_progress("\n> full run of Femocs", success == 0);

  free(x);
  free(y);
  free(z);

  return 0;
}