# Physics 427 Homework 4

__Due 11:59pm Monday 9/25/2023__

## 1. Implementing RK4

In the repository, I have included an example implementation of the Euler's
method, `euler.h`. Take a look at the code, and try to understand it as much as
possible. I have included a skeleton file for a 4th-order Runge-Kutta solver,
`rk4.h`. You task is to fill out the function bodies with a TODO comment and finish the implementation.

To help you debug your code, I also included a simple test file `test1.cpp` to
run the solver through the simplest equation, $dy/dx = \cos x$. Its output,
`output1.csv`, can be checked against your expected solution (what is the
expected solution?) to see whether your RK4 is implemented correctly.

Note that RK4 requires some intermediate values during evaluation, and we can store them in the class. It is convenient to define `k1`, `k2`, `k3`, `k4`, and a temporary `y` vector to store the estimated midpoint positions. These should be initialized to the number of equations you are solving, so we resize them to the right length in the _constructor_ of the class. The constructor of a class is a function with the same name as the class name, and no return value. Its job is usually to initialize the various data members of the class:

``` c++
class runge_kutta_4 {
public:
  runge_kutta_4(int num_equations) {
    n_eq = num_equations;
    y_tmp.resize(n_eq);
    k1.resize(n_eq);
    k2.resize(n_eq);
    ...
  }
};
```

Commit your changes to the `rk4.h` file to the repository.

## 2. The Nonlinear Pendulum Revisited

For our first physics application, let's solve the complete equation of motion of the nonlinear pendulum we looked at in Homework 3. The equation of motion can be written as:

$$
\frac{d^2\theta}{dt^2} + \frac{g}{L}\sin\theta = 0
$$

If we define $x = t/\sqrt{L/g}$, then $x$ is dimensionless (see Lecture 4). The equation becomes:

$$
\frac{d^2\theta}{dx^2} = -\sin\theta
$$

This is a 2nd order ODE. To convert it to a system of 1st order ODEs, we will introduce a new unknown function $\omega = d\theta/dx$. Then, the two-equation system can be written as:

$$
\frac{d\theta}{dx} = \omega,\quad \frac{d\omega}{dx} = -\sin\theta
$$

Your task is to use both Euler's method and RK4 to solve the equations above, and compare the results. Choose your initial condition such that $x_0 = 0$ and $\theta_0 = \theta_\mathrm{max}$ (what is the initial $\omega$ then?). For maximum deviation from the regular pendulum, set $\theta_\mathrm{max} = 179^\circ$ (remember to convert it to radians in the code). Use a step size where these two methods give you visibly different results. Run for 5 periods. Remember you can use your program from last week to calculate the period. This is a place where $N=10$ Gauss quadrature might give you a not-so-accurate answer. Try adaptive Simpson and target a high precision like $10^{-10}$.

Finally, plot both curves of $\theta$ vs $x$ in Python and include the plot in the repo as `problem2.png`. Remember to add an $x$ label and $y$ label to the axes of the plot to indicate what they are! Also write your step size on the plot as the title. Add your source file to the repo as `problem2.cpp`, and your python file as `problem2.py` or `problem2.ipynb`. We would like someone with the code to be able to reproduce your plot. 

### Setting up a Python environment

If you already have a Python environment, then feel free to skip this section. Here are the ways you can set up Python on Mac, WSL, and Linux. I strongly suggest Python 3, since Python 2 is already deprecated. 

__On Mac__, you can install Python through `homebrew`:

``` sh
brew install python
```

Note that if you used `homebrew` to install `gcc` in HW0, and have been using `g++` since then, it is very likely that you are actually still using the default `g++` that comes with Mac. To force use the version installed by `homebrew`, you can call `g++-13` instead.

__On WSL__, you are effectively in a Ubuntu system. Install the Python package using:

``` sh
sudo apt-get install python3 python3-pip
```

__On Liunx__, use the package manager of your distribution. 

Once you have installed Python, install the required Python packages using `pip` (or `pip3` for ubuntu):

``` sh
pip install numpy matplotlib jupyter
```

These packages should get you started. You may want to also install the Jupyter plugin for VS Code, so that you can run notebooks directly in VS Code. Install other packages as we go if need be.

## 3. RK4 with Adaptive Step Size

Now let's modify your RK4 implementation to include the simplest adaptive step size control. Starting from a skeleton class I defined in the file `rk4_adaptive.h`. You will need to implement the parts of the code that are marked with "TODO:". 

The central idea of the adaptive RK4 is that, instead of simply doing a Runge-Kutta step at every loop, carry out two steps, one with step size $h$ and the other with two ministeps half the step size, $h/2$. The difference between these two results gives you an estimate for the error on $y_i$.

$$
\varepsilon_i = |y_2 - y_1|_i
$$

Rescale your error to respect to some target absolute tolerance `atol` and relative tolerance `rtol`, using the following definition:

$$
\mathrm{err} = \sqrt{\frac{1}{N_\mathrm{eq}} \sum_{i=0}^{N_\mathrm{eq} - 1}\left(\frac{\epsilon_i}{\mathrm{atol} + |y_i|\mathrm{rtol}}\right)^2}
$$

You can think of this rescaled error as a ratio between the estimated error and the target error computed using our tolerance. If `err > 1` then the estimated error is larger than our target, so we reject the step and run again with a smaller step size. Otherwise, accept the step, and increase the step size.

Commit your changes to `rk4_adaptive.h` to the repository. Use it to solve problem 2 again and produce a plot `problem3.png`. In this plot, draw a vertical dashed line at every output point $x$. This should give you an idea when the method takes larger steps and when it takes small steps.

## 4. (Optional, no extra credit) Dormand-Prince RK45 Adaptive Method

This problem is optional and will not count towards the grade of this homework (even for grad students). Do this if you want a challenge, or if you want a (relatively) simple and customizable ODE solver for the rest of your career.

Implement the Dormand-Prince version of the adaptive Runge-Kutta method we covered in class. It should be very similar to the adaptive RK4 you implemented in Problem 3. You can start with the same skeleton file `rk4_adaptive.h`. There are just many more coefficients that you'll need to type in by hand (and check their correctness).

Name your implementation `rk45_dormand_prince.h`, and commit it to the repository. 


