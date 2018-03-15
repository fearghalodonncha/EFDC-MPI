#!/bin/bash
# This file is only used by the Vagrantfile to configure installations
# If one would like additional software installed, one can edit this file
# (either to install from package manager or source)

# USER=efdc

# HOME=/home/${USER} 

echo "Provisioning virtual machine..."
sudo apt-get update

echo "Installing required packages"
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y > /dev/null
    sudo apt-get update -y > /dev/null
    sudo apt-get install automake apt-utils gcc gfortran openssh-server wget  -y > /dev/null
    sudo apt-get install python-dev python-numpy python-pip python-matplotlib python-tk -y > /dev/null
    sudo apt-get install git make m4 zlib1g-dev -y > /dev/null
    sudo apt-get install libopenmpi-dev openmpi-bin libhdf5-openmpi-dev openmpi-common -y >  /dev/null
    sudo apt-get install libnetcdf-dev libblas-dev liblapack-dev -y > /dev/null
    sudo apt-get install libfreetype6-dev pkg-config libpng12-dev -y > /dev/null
    export PATH=$PATH:/usr/bin/
    

# Add NetCDF & OpenBlas libraries
# 1) We need HDF
#Required for NetCDF integration; source is supported at http://www.hdfgroup.org/ftp/HDF5/current/src

HDF_VERSION="1.10.1"  
# Download build and install HDF5
    mkdir ~/temp
    cd ~/temp
    wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-${HDF_VERSION}.tar.bz2; 
    tar -xjvf hdf5-${HDF_VERSION}.tar.bz2; 
    cd hdf5-${HDF_VERSION}; 
    ./configure --enable-shared --prefix=/usr/local/hdf5; 
    make;  
    sudo make install; 
    cd ..;  
    rm -rf /hdf5-${HDF_VERSION} /hdf5-${HDF_VERSION}.tar.bz2; 


# 2)
#Build netcdf 
# First we need to build NetCDF C version 
# (http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)
# Download and install netcdf C

   NCD_VERSION="4.3.3.1"
   NCD_VERSION="4.6.0"
   wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-${NCD_VERSION}.tar.gz;  
   tar -xzvf netcdf-${NCD_VERSION}.tar.gz; 
   cd netcdf-${NCD_VERSION}; 
   ./configure --prefix=/usr/local/netcdf CC=gcc LDFLAGS=-L/usr/local/hdf5/lib CFLAGS=-I/usr/local/hdf5/include --disable-dap; 
   make ; 
   sudo make install; 
   cd .. ;
   rm -rf netcdf-${NCD_VERSION} netcdf-${NCD_VERSION}.tar.gz

# 3) Build NetCDF fortran version
# (http://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html)
# Download and install NetCDF fortran
     NCF_VERSION="4.4.2"
     NCF_VERSION="4.4.4"
     wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-${NCF_VERSION}.tar.gz; 
     tar -xzvf netcdf-fortran-${NCF_VERSION}.tar.gz;  
     cd netcdf-fortran-${NCF_VERSION}; 
    ./configure --prefix=/usr/local/netcdf  --disable-fortran-type-check CC=gcc FC=gfortran  LDFLAGS=-L/usr/local/netcdf/lib CFLAGS=-I/usr/local/netcdf/include; 
     make ;  
     sudo make install; 
     cd ..  
     rm -rf netcdf-fortran-${NCF_VERSION}  netcdf-fortran-${NCF_VERSION}.tar.gz

## Add netcdf path to vagrant machine 
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/netcdf/lib" >> /home/vagrant/.bashrc
echo "export PATH=$PATH:/usr/local/netcdf/bin" >> /home/vagrant/.bashrc

## Include NCDIR and HDF5DIR directories for python netCDF4 installation
echo "export HDF5_DIR=/usr/local/hdf5"  >> /home/vagrant/.bashrc
echo "export NETCDF4_DIR=/usr/local/netcdf"  >> /home/vagrant/.bashrc

# Update any python requirements
sudo pip install --upgrade pip 
sudo pip install --upgrade setuptools 
sudo pip install -U matplotlib 
sudo pip install -U numpy 
sudo pip install -U netCDF4 
 