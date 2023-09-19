#pragma once

#include <vector>
#include <cmath>
#include <iostream>

class rk4_adaptive {
public:
  rk4_adaptive(int num_equations, double tolerance_abs, double tolerance_rel) {
    // constructor, initialize the data members
    n_eq = num_equations;
    atol = tolerance_abs;
    rtol = tolerance_rel;
    y_tmp.resize(n_eq);
    y_err.resize(n_eq);
    k1.resize(n_eq);
    k2.resize(n_eq);
    k3.resize(n_eq);
    k4.resize(n_eq);
  }

  template <typename F, typename StopCondition>
  std::vector<double> integrate(const F &f, const StopCondition &stop, double h,
                                double x0, const std::vector<double> &y0) {
    // clear the arrays so that we start fresh every time we call this function
    xs.clear();
    result.clear();

    double x = x0;
    std::vector<double> y = y0;
    double err = 0.0;

    // Always push the initial condition into the results
    xs.push_back(x);
    result.push_back(y);

    while (stop(x, y) == false) {
      y = step(f, h, x, y);
      err = error(y);
      // If err is fortuitously too small, set it to some lower bound
      err = std::max(err, 1.0e-10);

      // Accept the step if the scalar error is below 1, otherwise reject it and
      // do not move forward
      if (err < 1.0) {
        x += h;
        xs.push_back(x);
        result.push_back(y);
      }

      // Adjust h as needed

      // TODO: implement this part where we adjust h. Increase h if err is small
      // and decrease it when err is large. Do not let h go below hmin or above hmax

      // std::cout << "x = " << x << ", h = " << h << ", err = " << err << std::endl;
    }
    return y;
  }

  template <typename F>
  std::vector<double> step(const F& f, double h, double x, const std::vector<double> &y) {
    // Compute the next step in y, given x and y of the current step

    // TODO: Compute two estimates for the next y, one with step size h, the
    // other with two steps each of size h/2. Try to reuse the function you
    // implemented for rk4.h

    // TODO: Estimate y_err for each y in the vector using the difference
    // between y1 and y2
  }

  // TODO: You might want to copy your rk4 "step" function here, call it
  // step_rk4, so that you can reuse it to compute y1 and y2

  double error(const std::vector<double> &y) {
    // TODO: compute a scalar error from the error estimate y_err and return it.
  }

  int n_eq;
  double atol, rtol;
  double hmin = 1.0e-10;
  double hmax = 1.0;
  std::vector<double> k1, k2, k3, k4, y_tmp, y_err;
  std::vector<double> xs;
  std::vector<std::vector<double>> result;
};
