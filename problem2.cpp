#include "euler.h"
#include "rk4.h"
// #include "rk4_adaptive.h"
#include <cmath>
#include <fstream>

int main() {
  auto dydx = [](double x, const std::vector<double> &y) {
    std::vector<double> result(2);
    result[0] = y[1];
    result[1] = -std::sin(y[0]);
    return result;
  };

  auto stop = [](double x, const std::vector<double> &y) { return x > 24.511146 * 5; };

  euler_method euler;
  runge_kutta_4 rk(2);
  // rk4_adaptive rka(2, 1e-12, 1e-12);
  double h = 0.01;
  double th0 = 179.0 / 180.0 * M_PI;
  auto y = rk.integrate(dydx, stop, h, 0.0, {th0, 0.0});
  // y = rka.integrate(dydx, stop, h, 0.0, {th0, 0.0});
  y = euler.integrate(dydx, stop, h, 0.0, {th0, 0.0});

  // std::ofstream output_file("output-rka.csv");
  // for (int i = 0; i < rka.xs.size(); i++) {
  //   output_file << rka.xs[i] << "," << rka.result[i][0] << std::endl;
  // }
  // output_file.close();

  std::ofstream output_rk("output-rk.csv");
  for (int i = 0; i < rk.xs.size(); i++) {
    output_rk << rk.xs[i] << "," << rk.result[i][0] << std::endl;
  }
  output_rk.close();

  std::ofstream output_euler("output-euler.csv");
  for (int i = 0; i < euler.xs.size(); i++) {
    output_euler << euler.xs[i] << "," << euler.result[i][0] << std::endl;
  }
  output_euler.close();

  return 0;
}
