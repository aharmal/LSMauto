# LSMauto

**LSMauto (Lattice Spring Model automation)** is a python code that simulates fracture of single edge notched beam in 3-point bending using Lattice Spring Modeling from one command line. Several functions, which are usually used separately, are chained, making the simulation process easy to use. The objective of this code is to make this simulation method more accessible and easier to use for the materials science community. The simulation is calibrated for our samples (geopolymer composites, see reference below). Feel free to explore the functions used in this code and taylor it to your need.

![figure 1](https://user-images.githubusercontent.com/99771175/229556965-b8fd9a3d-34df-4529-ac82-106575c75bbc.jpg)

Please cite work related to this code (Find our latest work in our scholar profile: https://scholar.google.com/citations?user=OZMm7ogAAAAJ&hl=en&oi=ao):
> Harmal, A., Khouchani, O., El-Korchi, T., Tao, M., & Walker, H. W. (2023). Bioinspired brick-and-mortar geopolymer composites with ultra-high toughness. Cement and Concrete Composites, 104944.


Please also cite work related to Img2Particle code used in this code.
> Chiang Y, Chiu T-W and Chang S-W (2022) ImageMech: From Image to Particle Spring Network for Mechanical Characterization. *Front. Mater.* 8:803875. doi: 10.3389/fmats.2021.803875

## requirement

- Operating system: Below directions and examples are for Linux. If you are using Windows, you can install Linux on Windows with WSL (https://learn.microsoft.com/en-us/windows/wsl/install). You can also run LSMauto on your Windows command line, provided you follow equivalent directions for Windows.

- Simulation software: Install a pre-built LAMMPS distribution (for Ubuntu for example)

```sudo apt-get install lammps```

Alternatively, for more experienced users (if needed) you can build LAMMPS in your system using make or cmake (https://docs.lammps.org/Install.html). In this case, please make sure to install MOLECULE and MC (Monte Carlo) packages.

- Python packages: Install Python in your terminal and install the following python packages: os, math, PIL, numpy, matplotlib, datetime, argparse.

## How to run 

- Download LSMauto and save it into a directory (pathtoLSMauto for example)

- Navigate to pathtoLSMauto by changing directory

```cd ../pathtoLSMauto```

- Execute LSMauto to run the simulation

```python3 ./LSMauto.py```

The previous command will use the standard values that we set up for the simulation variables. You can change the parser standard values.
It is also possible to give different values to the variables using the command line as follows:

```python3 ./LSMauto.py -ar 4 -ov 2 ```

The previous command will change the aspect ratio to 1:4 and the overlap to 0.2

- Once the simulation is done, you can open the dump files using the open source software OVITO (https://www.ovito.org/about/).
You can also open the text file ss.txt to see displacement-load results at each timestep, and use it to plot strain-stress graphs.
The dump files give the position (x,y,z), stresses (σxx, σyy, σxy), and forces (fx, fy) for all the datapoints of the simulated system. This gives you the ability to compute other properties we did not cover in our post-processing.

- If you would like to take advantage of parallel computing, please change the line in LSMauto (in case you have 10 cores for example): 

```os.system("lmp -in script.in")```

into: 

```os.system("mpiexec -np 10 lmp -in script.in")```
