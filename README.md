# Extreme_Env_Nav_Suite

**Author:** Heli Corgi  
**Purpose:** Learning/demo toolkit for navigation & environment adaptation in C++

**Note:** This project was developed with the assistance of an AI coding assistant to organize and expand the demo.

## Overview

Extreme_Env_Nav_Suite is a simplified C++ toolkit for educational purposes, demonstrating:

- **NavEngine_L1**: 3D L1 navigation engine (simplified, adaptive control demo)  
- **NavEstimator_SRUKF**: Skeleton of Square-root Unscented Kalman Filter for state estimation  
- **EnvAdaptive_Core**: PCM-based thermal model for environmental adaptation  

**Important:** This code is for learning/demo purposes only. It is not intended for real navigation or mission-critical applications.

## Features

- 3D navigation simulation with position and velocity in X, Y, Z axes
- SRUKF-based state estimation of the navigation engine
- Parallel PCM thermal model simulation
- Demonstrates sigma points, Cholesky updates, and simple adaptive control principles

## Build Instructions

### Dependencies

- C++17 compiler
- [Eigen3](https://eigen.tuxfamily.org/) (required)

### Using CMake

```bash
git clone https://github.com/HeliCorgi/X-Env-Nav-Suite-3D-L1-Navigation-SRUKF-Estimation-Lab.git
cd X-Env-Nav-Suite-3D-L1-Navigation-SRUKF-Estimation-Lab
mkdir build && cd build
cmake ..
cmake --build .
./nav_suite_run   # or nav_suite_run.exe on Windows
```

## Example Output

```
[t=0.00s] Raw: pos=(0.000,0.000,0.000) vel=(0.000,0.000,0.000) scale=1.000
     Est: pos=(0.001,0.000,0.000) vel=(0.000,0.000,0.000)
     PCM=100% Temp=20.024℃
...
```

## What You Can Learn

- How a simple L1 adaptive navigation engine works in 3D
- Basic usage of Unscented Kalman Filters for state estimation
- Eigen3 matrix operations and sigma point manipulation
- Interaction between a navigation engine and environmental constraints (thermal/PCM model)

## License

MIT

---

💡 **Key Points**

- main.cpp is ready to copy/paste and run
- NavEngine_L1 is fully 3D compatible
- SRUKF provides foundational concepts for state estimation learning

## References

- Atmosphere model:  
  U.S. Standard Atmosphere, 1976 (NASA/NOAA/USAF).  
  NASA TM X-74335 / NOAA-S/T-76-1562 (1976).  
  https://ntrs.nasa.gov/citations/19770009539

- Convective heating formula:  
  Sutton, K. and Graves, R. A., Jr., "A General Stagnation-Point Convective-Heating Equation for Arbitrary Gas Mixtures," NASA Technical Report R-376 (November 1971).

- Gravity model:  
  World Geodetic System 1984 (WGS84), NIMA TR8350.2.

- Coriolis term: Standard geophysical approximation (2ω × v).
