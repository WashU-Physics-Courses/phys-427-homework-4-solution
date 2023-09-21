#pragma once

#include <vector>

class runge_kutta_4 {
public:
  runge_kutta_4(int num_equations) {
    // constructor, initialize the data members
    n_eq = num_equations;
    y_tmp.resize(n_eq);
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
    while (stop(x, y) == false) {
      // append x and y to the result array
      xs.push_back(x);
      result.push_back(y);

      // compute the next step in y
      y = step(f, h, x, y);
      x += h;
    }
    return y;
  }

  template <typename F>
  std::vector<double> step(const F& f, double h, double x, const std::vector<double> &y) {
    // Compute the next step in y, given x and y of the current step

    // TODO: Finish this implementation, refer to the euler.h if you are
    // unsure where to start.
    std::vector<double> y_next = y;
    k1 = f(x, y);

    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + 0.5 * h * k1[i];
    }
    k2 = f(x + 0.5 * h, y_tmp);

    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + 0.5 * h * k2[i];
    }
    k3 = f(x + 0.5 * h, y_tmp);

    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * k3[i];
    }
    k4 = f(x + h, y_tmp);

    for (int i = 0; i < n_eq; i++) {
      y_next[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    return y_next;
  }

  int n_eq;
  std::vector<double> k1, k2, k3, k4, y_tmp;
  std::vector<double> xs;
  std::vector<std::vector<double>> result;
};
