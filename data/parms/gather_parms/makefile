#!/bin/bash

INC_NETCDF=/discover/nobackup/projects/lis/libs/netcdf/4.3.3.1_intel-14.0.3.174_sp3/include/
LIB_NETCDF=/discover/nobackup/projects/lis/libs/netcdf/4.3.3.1_intel-14.0.3.174_sp3/lib/
INC_HDF5=/discover/nobackup/projects/lis/libs/hdf5/1.8.14_intel-14.0.3.174_sp3/include/
LIB_HDF5=/discover/nobackup/projects/lis/libs/hdf5/1.8.14_intel-14.0.3.174_sp3/lib/


LDFLAGS = -g -L$(LIB_NETCDF) -L$(LIB_HDF5) 
FFLAGS = -I$(INC_NETCDF) -I$(INC_HDF5) #-I$(INC_HDF4)
FC = ifort

EXEFILE = main.exe
OBJS = gather_parameters.o
SRCS = gather_parameters.f90

all: $(EXEFILE) partclean

$(EXEFILE): $(OBJS) 
	$(FC) -o $(EXEFILE) $(OBJS) -L$(LIB_NETCDF) -lnetcdff -lnetcdf  #-L$(LIB_HDF5) -lhdf5_fortran -lhdf5_hl -lhdf5  
$(OBJS):  
	$(FC) -c -g $(FFLAGS) $(SRCS) 

partclean:
	/bin/rm $(OBJS)

clean: 
	/bin/rm $(OBJS) $(EXEFILE)  


# DO NOT DELETE THIS LINE -- make depend depends on it.
