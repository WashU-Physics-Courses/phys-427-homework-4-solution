#include "rk4.h"
#include <cmath>
#include <fstream>

int main() {
  // Define our right hand side
  auto dydx = [](double x, const std::vector<double> &y) {
    return std::vector<double>{std::cos(x)};
  };

  // Define our stop condition
  auto stop = [](double x, const std::vector<double> &y) { return x > 10; };

  // Integrate the equation until the stop condition
  runge_kutta_4 rk4(1);

  double h = 0.01;
  double x0 = 0.0;
  std::vector<double> y0 = {0.0};
  rk4.integrate(dydx, stop, h, x0, y0);

  // Output the results
  std::ofstream output_file("output1.csv");
  for (int i = 0; i < rk4.xs.size(); i++) {
    output_file << rk4.xs[i] << "," << rk4.result[i][0] << std::endl;
  }
  output_file.close();

  return 0;
}
