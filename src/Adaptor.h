#include <cstddef>
#include <vector>

class Grid {
public:
  Grid();
  void Initialize(const unsigned int numPoints[3], const double spacing[3]);
  size_t GetNumberOfPoints();
  size_t GetNumberOfCells();
  double *GetPointsArray();
  double *GetPoint(size_t pointID);
  unsigned int *GetCellPoints(size_t cellID);

private:
  std::vector<double> Points;
  std::vector<unsigned int> Cells;
};

class Attributes {
public:
  Attributes();
  void Initialize(Grid *grid);
  void UpdateFields(double time);
  double *GetVelocityArray();
  float *GetPressureArray();

private:
  std::vector<double> Velocity;
  std::vector<float> Pressure;
  Grid *GridPtr;
};

namespace Adaptor {
void Initializa(int numScripts, char *scripts[]);
void Finalize();
void CoProcess(Grid &grid, Attributes &attributes, double time,
               unsigned int timeStep, bool lastTimeStep);
} // namespace Adaptor
