# Robofly_simulator
This repo is for devloping a 3D simulator for flapping wing robots on python. The controller is originaly written by Dr. Sawyer B. Fuller, coverted to python by Daksh Dhingra.

# Installation
Required libraries for running the controller

- numpy
- control
- csv
- matplotlib
- scipy

# Files 
- 'controller.py' - class with lateral and altitude controller methods
- 'robot.py' - class defining the dynamics of roboFly. Contains methods for calculating drag force and drag torque calculation.
- 'simulator.py' - Simulator class containing methods for plot, save data and ODE solver.
- 'test.py' - main fuction to define initial and desired state and call ODE solver.
- 'load_and_plot.py' - script used to plot the saved data.


# Refrences
Mellinger, D., Michael, N., & Kumar, V. (2012). Trajectory generation and control for precise aggressive maneuvers with quadrotors. The International Journal of Robotics Research, 31(5), 664–674.

McGhee, Robert, Eric Robert Bachmann, and Michael J. Zyda. "Rigid body dynamics, inertial reference frames, and graphics coordinate systems: A resolution of conflicting conventions and terminology." (2000). 
