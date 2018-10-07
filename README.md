# Numerical methods in physics 2016-2017/1 @ ELTE

All of the codes were wrote in `C language`.

During the course, 4 projects had to be written, which were the following:

### 1. Vector-Matrix Multiplication

**Task description issued during the course:**
*Full task*: Write an algorithm, which reads a matrix and a vector from two separate files. The matrix's size should not be pre-initialized, but the program should find it out of the number of columns and rows in the file. There is not much bound to the file format, the main thing is to have a text file, which somehow "looks like" a matrix. Prepare different test files for the program to test the correct operation. **Please note that the files can only be read once and only sequentially and the files can be arbitrarily large.**


### 2. Gauss-Jordan Elimination

**Task description issued during the course:**
*Minimum task*: Write a program that scans a NxN-sized matrix (from text, space, or line-separated files), with real elements and outputs the inverse of the matrix as a result. As an algorithm, use the Gaussian elimination with row and column exchange. The program also should detect if the matrix to be inverted is singular!

*Full task*: After implementing the Gaussian elimination, use the singular value decomposition function of the LAPACK packet to determine the inverse of the matrix. Compare the result of the SVD algorithm with Gaussian elimination for some matrices! Test the speed of the two algorithms for matrices of different sizes! Let's try SVD for singular and quasi-singular matrices!


### 3. Polynomial Fitting

**Task description issued during the course:**
*Minimum task*: Write a program, using the linear equation system solving code elaborated in the previous task, which performs a general N variable polynomial fitting, using the linear least squares method. The program should wait for the number of independent variables, the order of the polynomial to be copied, and the file name of the input data, as command line parameters. The program should display the fitting parameters on the standard output, and should write a new text file, which cointains the original measurement values, and the estimates provided by the matched polynomial for the datapoints of each original input files. Make a plot of the fitted polynomial on the datapoints!
Optionally, create your own data files to test the algorithm, or just run the program for the data files below (which contain 5 variables).
If you make the minimum task, you don't need to take into account the "mixed-members" of multi-variable polynomials (eg. xy, xy2, etc...).

*Optional (not graded) task*: Compute non-linear function matching using the MCMC method, the Metropolis-Hastings algorithm. To test the program, create your own data files. Try Gaussian fit!


### 4. Solving Differential Equation

**Task description issued during the course:**
*Minimum task*: Write a program that uses the explicit Euler method to solve the equation of motion of the Earth-Moon system. Make a plot of the obtained results: show the Moon's track and the total energy of the system as a function of time! Check if there is any kind difference from the analytical solution after many periods! (Eg. is the conservation of energy violated or not?) 
When working out the problem, work in plane, with geocentric coordinates (Earth is fixed at 0.0 points).

Given constants:
Mass of Earth: $5.9736 \times 10^{24} kg$
Mass of Moon: $7.349 \times 10^{22} kg$
Moon's distance in Apogee: $405.500 km$
Moon's speed in Apogee: $964 \frac{m}{s}$
Moon's distance in Perigee : 363.300 km
Moon's speed in the Perigee : $1076 \frac{m}{s}$
Gravitational constant: $6.67384 \times 10^{-11} m^{3} kg^{-1} s^{-2}$

*Full task*: Integrate the equation of motion with a simple, and then an adaptive step-length-controlled RK4/classical Runge–Kutta method. The program should be completely general, no constraint for the variables to be integrated. Use functionpointers! Try the program to solve the equations of the Lorenz equation or the double pendulum. Show the results on a plot!
