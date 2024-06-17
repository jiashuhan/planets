# planets
### A Newtonian N-body simulation 
<img src="./results/videos/inner-solar-3d.gif" width="600" height="375"/>
<img src="./results/videos/inner-solar-2d.gif" width="600" height="300"/>

To start the simulation, run

```
python3 run.py
```

Add an object to the list by entering its __name__, __mass__ and the __state vectors__ (position and velocity). Alternatively, one can enter the Keplerian orbital elements and convert them into the state vectors. In this case, the __parent body__ must be specified (must be either 'N/A' or one of the objects already added to the simulation). The orbital parameters of the eight planets (plus Pluto and two comets) in the Solar System (epoch = 2451545.0 (2000-Jan-1.5)) are saved in the ```parameters``` folder and can be loaded into the program. Once all the objects are entered, set the __duration__ (in days) and __time resolution__ (steps per day) of the simulation before starting the run. The simulation integrates using a fixed-stepsize Dormand-Prince method by default, but other options are available:

- Symplectic Euler
- Velocity Verlet
- Leapfrog
- 4th order Runge-Kutta (RK4)
- Dormand-Prince (RK5(4)) with fixed stepsize (Default)
- Dormand-Prince (RK5(4)) with adaptive stepsize (```scipy.integrate.RK45```)

In the case of the adaptive Dormand-Prince method, a __maximum local position error__ ("Tolerance", in meters) should be specified instead of the time resolution.

For the solar system, which is sufficiently stable, using the least accurate symplectic Euler method with a time resolution of 10 steps per day should give reasonable accuracy. The program automatically saves the simulation snapshots in the ```results``` folder. To calculate the orbital elements (_e_, _a_, _i_, $\Omega$, $\bar{\omega}$, _L_) or (_e_, _a_, _i_, $\Omega$, $\omega$, _M_) of a particular object at the end of the simulation, run

```
python3 kepler.py simulation_output object_name central_body_name
```

The final state of the simulation run is plotted and can be found in ```results/plots```. To create videos from the snapshots, run

```
python3 visualization.py simulation_output plot_range_in_m colormap kepler_overlay
```

This will generate videos of the simulation in ```results/videos``` with the specified plot range (in meters) and color map. If the parameter ```kepler_overlay``` is set to 1, the Keplerian orbits specified by the initial orbital elements will be plotted in dashed lines alongside the actual orbits. Only orbits of objects directly orbiting the central object (i.e. no moons) will be plotted.

<img src="./results/videos/ksp.gif" width="600" height="300"/>

In addition, there is the option to generate a number of objects with random masses, positions and velocities within a specified radius, which can then be added to the system:

<img src="./results/videos/random.gif" width="600" height="375"/>
