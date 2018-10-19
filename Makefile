#### GNU compilers ####
FFLAGS = -O3
WARN  = -Wall
#BOUND = -fbounds-check
f2pyCC = gfortran
CC = mpifort
fexec = mpif90
LOG = > log.txt

#### NERSC ####
#FFLAGS = -O3 
#BOUND = 
#f2pyCC = intelem
#CC = ftn 
#fexec = ftn

BASE0 = ConstData Math pGrid2D fileIO
BASE1 = CompDom
BASE2 = Radiation
BASE3 = eBeam
BASE4 = Element
BASE5 = Undulator Drift Quad Propagator Dipole DiagnoElems
BASE6 = ImpactFEL
BASE  = $(BASE0) $(BASE1) $(BASE2) $(BASE3) $(BASE4) $(BASE5) $(BASE6)

SRC0 = $(addsuffix .f90,${BASE0})
SRC1 = $(addsuffix .f90,${BASE1})
SRC2 = $(addsuffix .f90,${BASE2})
SRC3 = $(addsuffix .f90,${BASE3})
SRC4 = $(addsuffix .f90,${BASE4})
SRC5 = $(addsuffix .f90,${BASE5})
SRC6 = $(addsuffix .f90,${BASE6})

OBJ0 = $(addsuffix .o,${BASE0})
OBJ1 = $(addsuffix .o,${BASE1})
OBJ2 = $(addsuffix .o,${BASE2})
OBJ3 = $(addsuffix .o,${BASE3})
OBJ4 = $(addsuffix .o,${BASE4})
OBJ5 = $(addsuffix .o,${BASE5})
OBJ6 = $(addsuffix .o,${BASE6})
OBJs = $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)

###### ImpactFEL library compile  #######
impactFEL.a : $(OBJ6)
	ar cr impactFEL.a $(OBJs) 

############### individual steps ###################
$(OBJ0) : 
	$(CC) -c $(BOUND) $(FFLAGS)       $(SRC0)

$(OBJ1) : $(OBJ0)
	$(CC) -c $(BOUND) $(FFLAGS) $(WARN) $(SRC1)

$(OBJ2) : $(OBJ1)
	$(CC) -c $(BOUND) $(FFLAGS) $(WARN) $(SRC2)

$(OBJ3) : $(OBJ2)
	$(CC) -c $(BOUND) $(FFLAGS)       $(SRC3)

$(OBJ4) : $(OBJ3)
	$(CC) -c $(BOUND) $(FFLAGS) $(WARN) $(SRC4)

$(OBJ5) : $(OBJ4)
	$(CC) -c $(BOUND) $(FFLAGS)       $(SRC5)

$(OBJ6) : $(OBJ5)
	$(CC) -c $(BOUND) $(FFLAGS) $(WARN) $(SRC6)


############### test_run ###################

test_run : impactFEL.a
	$(CC) $(BOUND) $(FFLAGS) $(WARN) test_run.f90 -o test_run impactFEL.a 

########## clean ############## 
.PHONY: clean
clean:
	rm -f *~ *.o *.mod impactFEL.a

.PHONY: cleanMod
cleanMod:
	rm -f -r *.mod

.PHONY: cleanO
cleanO:
	rm -f -r *.o


