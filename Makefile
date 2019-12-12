
INCLUDES=-I ~/Code/dist/Debug/include/nspr/ -I ~/Code/nss/lib/freebl/mpi -I ~/Code/nss/lib/util/
LIBS=-L ~/Code/dist/Debug/lib
MPI_FILES=~/Code/nss/lib/freebl/mpi/mpprime.c ~/Code/nss/lib/freebl/mpi/mpi.c ~/Code/nss/lib/freebl/mpi/mpmontg.c ~/Code/nss/lib/freebl/mpi/mplogic.c ~/Code/nss/lib/freebl/mpi/mpcpucache.c

all:
	gcc ROCA.c -o ROCA -lgmp -std=c99 $(INCLUDES) $(MPI_FILES)
