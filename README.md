# Project Overview
This program solves time dependent Schrodinger equation within a box for two different scenarios zero and harmonic_well potentials 
using the Cranck-Nicolson method.


# Compilation Instructions
Compilation is all done using the makefile in the repository. Type `make` into your command line to compile the files.

Initially sets schrodinger as the executable name.
- Then creates the types object file.
- Then creates the linear_algebra object file using the types object file.
- Then creates the hamiltonian object file using the types object file.
- Then creates the quantum object file using the linear_algebra, hamiltonian and types object files.
- Then creates the read_write object file using the quantum and types object files.
- Then creates the main object file using the types, read_write, ode_solver and quantum object files.


# Usage Instructions 
Once you have compiled everything execute the program using ```./schrodinger arg```
The program can be excuted using just the executable ```./schrodinger``` (which will return results for default initial conditions) or it can take 
namelist files as optional arguments. Just make sure that you keep the namelist files in the directory you are working in. Make sure you set 
the `potential_scenario` to either `zero` or `harmonic_well` since those are the only two scenarios set for this program.

The namelist files need to be structured like:

```
&integration
  length = 5.0,
  n_points = 100,
  n_steps = 100,
  delta_t = 0.05
/
&wave_function
	width = 0.5,
  center = 0.0,
  potential_scenario = 'zero'
/
&oscillator
	k_oscillator = 0.0
/
&output
	time_file = 'zero_time_results.dat',
	density_file = 'zero_density_results.dat'
/

```

Make sure to include an extra line at the end of the file (after the last hash).
Once created you will be able to run the Jupyter notebook the data files produced. 


# Expected Behavior
Once you have executed it will perform the calculations, it will then writes the results of the program into two .dat files.

### Probability Density file
The first output file should have `n_steps` + 1 columns. 'x_points',....

- x_points: This column will contain the x position of the first object at each step.
- every column after: will contain the probabilty density at each point for every evolution of the wave function.

### Expectations file
If the `potential_scenario` is set to `zero`.
The second output file should have 5 columns. 'time', 'norm', 'position', 'width', 'analytic_width'.

- time: This column will contain the time at each step.
- norm: This column will contain the normalization for the probability densities at each time step.
- position: This column will contain the numerical position at each time step.
- width: This column will contain the numerical width at each time step.
- analytic_width: This column will contain the analytical width at each time step for comparison to the numerical result.

If the `potential_scenario` is set to `harmonic_well`.
The second output file should have 5 columns. 'time', 'norm', 'position', 'width', 'analytic_position'.

- time: This column will contain the time at each step.
- norm: This column will contain the normalization for the probability densities at each time step.
- position: This column will contain the numerical position at each time step.
- width: This column will contain the numerical width at each time step.
- analytic_position: This column will contain the analytical position at each time step for comparison to the numerical result.

