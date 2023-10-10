FC	= gfortran
FFLAGS  = -O

LD	= gfortran
LFLAGS  = -O

OBJ	= \
LJ.o coul.o diag.o eigsrt.o field.o force.o gauss.o gks-force.o setQ.o setup.o tqli.o tred2.o dstep.o

# =================================================================

clean:
	rm -f *.o 

# =================================================================

gks-force.x:	$(OBJ)
	$(LD) $(LFLAGS) $(OBJ) -o gks-force.x

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

.c.o:
	cc -c $<
