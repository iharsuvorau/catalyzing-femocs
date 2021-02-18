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
#include "CatalystAdaptor.h"

using namespace std;

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

int main(int argc, char *argv[]) {
    string filename = "in/md.in";

    femocs::Femocs project(filename);

    int success = 0;
    int n_iterations = 100;
    bool add_rnd_noise = true;

    CatalystAdaptor::Initialize(argc - 1, argv + 1);

    for (int iter_i = 1; iter_i <= n_iterations; ++iter_i) {
        if (n_iterations > 1)
            cout << "\n> iteration " << iter_i << endl;

        success = project.import_atoms("", add_rnd_noise);
        success += project.run();

        CatalystAdaptor::CoProcess(project, iter_i, iter_i, iter_i == n_iterations);
    }

    CatalystAdaptor::Finalize();

    print_progress("\n> full run of Femocs", success == 0);
    return 0;
}
