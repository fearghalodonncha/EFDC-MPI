# EFDC-MPI
This repo includes an extended version of the widely used Environmental Fluid Dynamics Code (EFDC) (https://www.epa.gov/exposure-assessment-models/efdc), a state of the art hydro-environmental modelling suite to simulate aquatic systems in one, two, and three dimensions. It serves as the hydrodynamic core of IBM Research asset DeepCurrent (http://www.research.ibm.com/labs/ireland/research_areas/deep_current.html), focusing on the prediction of environmental conditions in coastal oceans, rivers and lakes.

Extensions include capabilities to run in parallel using MPI, netCDF file I/O and incorporation of modules to simulate impeded flows.

Details on the parallel development of the code are provided here
[Computers & Geoscience](https://www.sciencedirect.com/science/article/pii/S009830041300304X) [[pdf]](https://www.researchgate.net/publication/259509004_Parallelization_study_of_a_three-dimensional_environmental_flow_model)
and applications here:
1) [Galway Bay](https://www.sciencedirect.com/science/article/pii/S0924796314002346) [[pdf]](https://www.researchgate.net/publication/268207331_Characterizing_observed_circulation_patterns_within_a_bay_using_HF_radar_and_numerical_model_simulations)
2) [Impeded aquaculture flows](https://www.tandfonline.com/doi/abs/10.1080/00221686.2015.1093036) [[pdf]](https://www.researchgate.net/profile/Fearghal_Odonncha/publication/283438617_Parameterizing_suspended_canopy_effects_in_a_three-dimensional_hydrodynamic_model/links/5a2133b1aca27229a06eb4b0/Parameterizing-suspended-canopy-effects-in-a-three-dimensional-hydrodynamic-model.pdf)
and 
3) [investigating marine renewable energy in cobscook Bay](https://www.sciencedirect.com/science/article/pii/S0960148116308898) [[pdf]](https://www.researchgate.net/publication/309306783_Modelling_study_of_the_effects_of_suspended_aquaculture_installations_on_tidal_stream_generation_in_Cobscook_Bay)

# Quickstart Guide
To build the EFDC model, clone this repository and take a peek at `QUICKSTART` (installs netCDF and MPI dependencies). You can directly execute it, given you are running a recent Ubuntu/Debian (and have sudo installed).

    $ git clone https://github.com/fearghalodonncha/EFDC-MPI.git
    $ cd Src/
    $ sudo ./Quickstart
    $ make

If netCDF or MPI libraries are installed in non-standard location, edit makefile to point to the installation path (NCDIR & MPICHDIR respectively).

If you are using a different distribution use your package manager to install all dependencies available. A list of dependencies can be viewed in the Dependencies section of this README.

The repo contains two sample model applications in the `SampleModels/` directory:

1) a simple harbour channel model
2) a real-world model of Cobscook Bay, ME, USA
3) a real-world model of Chesapeake Bay and
4) a simple harbour channel model with data assimilation

To run either of these examples, copy the `EFDC` executable from the `Src/` directory to the directory containing the `*.INP` files and run the code using ./EFDC (for serial examples)

# Dependencies
The dependencies for the EFDC model are NetCDF fortran (details on the installation are provided here https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html along with the associated dependencies) and mpi (e.g. openmpi (https://www.open-mpi.org/) or mpich (https://www.mpich.org/)). The `Quickstart` file included in the repo installs these on Ubuntu systems.

If one wishes to use the data assimilation libraries included in this repo then blas linear algebra libraries must be installed (e.g. http://www.openblas.net/) and edit makefile to point to installation path.

Once installed, the makefile included in the Src directory must be edited to point to the path of the MPI and netCDF libraries (`MPICHDIR` and `NCDIR` and `BLAS_PATH` respectively).

# Running examples
Sample input files are included in the `SampleModels/` directory. EFDC input files are of form `*.INP`
To run these examples in serial, copy the `EFDC` executable to the directory containing the `*.INP` files and run the code using ./EFDC

To run in parallel (using MPI), execute of form (across 4 compute cores):

$ mpirun -np 4 ./EFDC

Configurations for parallel simulations are defined in `LORP.INP` which describes the size of each sub-domain for distributed computing. This file is not specified by the user but rather a load balancing algorithm. The number of sub-domains specified in `LORP.INP` must equal the number passed to `mpirun`

A key consideration in parallel simulation of coastal ocean applications is distributing load across processors in a well balanced manner. This consists of ensuring that each sub-domain contains a relatively equal distribution of land/ocean cells (since land cells do not invoke any computational cost).

The repo contains a C++ load balancing algorithm that computes a rectilinear decomposition of a global domain into a number of subdomains in a locally optimal manner

To build the load balancing module, cd into the `Gorp` directory and compile using a C++ compiler:

  $ g++ Gorp.cpp GorpMain.cpp -o Gorp
    
Copy the Gorp executable into the same directory containing the input files and execute of form
  $ ./Gorp CELL.INP sc #

(Note: the CELL.INP header line must accurately specify the I and J extents of the domain as it is read by the load-balancing module).

Where # denotes  the number of subdomains to decompose the problem into (i.e. how many compute cores to distribute the problem over) while 'sc' is appended to the generated LORP file as an identifier. The generated files are of form `LORP_sc_#_v.INP`, `LORP_sc_#_h.INP` and `LORP_sc_#_r.INP` where h, v, and r denote decomposing the domain into balanced horizontal or vertical strips or into cartesian rectilinear domains respectively

Details on the parallel implementation and load balancing module are provided in:

O'Donncha, F, Ragnoli, E. and Suits, F. "Parallelisation study of a three-dimensional environmental flow model." Computers & Geosciences 64 (2014): 96-103.

# Details on the model
EFDC is a state-of-the-art hydrodynamic model that can be used to simulate aquatic systems in one, two, and three dimensions. It has evolved over the past decades to become one of the most widely used and technically defensible hydrodynamic models in the world. EFDC uses stretched or sigma vertical coordinates and Cartesian or curvilinear, orthogonal horizontal coordinates to represent the physical characteristics of a waterbody. It solves three-dimensional, vertically hydrostatic, free surface, turbulent averaged equations of motion for a variable-density fluid. Dynamically-coupled transport equations for turbulent kinetic energy, turbulent length scale, salinity and temperature are also solved. The EFDC model allows for drying and wetting in shallow areas by a mass conservation scheme. 

Extensions to this version of the code include MPI parallelism described in:

O'Donncha, F, Ragnoli, E. and Suits, F. "Parallelisation study of a three-dimensional environmental flow model." Computers & Geosciences 64 (2014): 96-103.

Inclusion of modules to represent aquaculture structures described in:

O'Donncha, F., Hartnett, M., & Plew, D. R. (2015). Parameterizing suspended canopy effects in a three-dimensional hydrodynamic model. Journal of Hydraulic Research, 53(6), 714-727.

And inclusion of modules to represent marine hydrokinetic structures described in:

O'Donncha, F., James, S. C., & Ragnoli, E. (2017). Modelling study of the effects of suspended aquaculture installations on tidal stream generation in Cobscook Bay. Renewable Energy, 102, 65-76.

# Vagrant example

To assist users in getting started with the model the repo includes a `Vagrantfile` that instantiates a Ubuntu machine and installs all required dependencies. To use the Vagrantfile, users must have Virtualbox and Vagrant installed

Install Virtualbox:
https://www.virtualbox.org/wiki/Downloads
Installing Vagrant:
Download from: https://www.vagrantup.com/downloads.html

After successfully installing vagrant

Run: vagrant plugin install vagrant-vbguest

For Windows - Run: vagrant plugin install vagrant-share --plugin-version 1.1.8

To bring up a Vagrant virtual machine, type:

    $ vagrant up
    $ vagrant ssh

Within the vagrant VM, users can easily build the model by simply invoking:

    $make

# Documentation

The repo contains a set of manuals describing the different components of the EFDC model. The fundamental hydrodynamic model input files and setup are described in the [EFDC User manual from the EPA](https://github.com/fearghalodonncha/EFDC-MPI/blob/master/Doc/Manuals/EFDC_User_Manual_EPA_Ver-101.pdf). Also included is a half day tutorial related to the EFDC MPI model and developing parallel model simulation presented at the [IEEE Oceans 2017 conference in Anchorage](http://www.oceans17mtsieeeanchorage.org/wp-content/uploads/2017/08/EFDC-MPI-IBM_tutorial_description_final.pdf). Further tutorial material related to the development of marine hydrokinetic applications is available from [Zenodo](https://doi.org/10.5281/zenodo.1292498).


  
