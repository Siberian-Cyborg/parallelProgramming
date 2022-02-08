CC=c++
MPICC=mpicxx

CXX_FLAGS=-O3 -Wall -Wextra -mavx
MPICXX_FLAG=-DOMPI_SKIP_MPICXX

STD=c++14

OPENMP=-fopenmp
OMPDEFINE=BUILDOMP
MPIDEFINE=BUILDMPI
HYBRIDDEFINE=BUILDHYBRID
SERIALDEFINE=BUILDSERIAL

VISUALDEFINE=VISUALIZATION

LIBS=-lm

utility.o: utility.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -c -o $@ $^

serial.o: serial.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -c -o $@ $^

visualserrial.o: serial.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -D$(VISUALDEFINE) -c -o $@ $^

ompkernel.o: ompbh.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) $(OPENMP) -c -o $@ $^

ompinterface.o: omp_interface.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) $(OPENMP) -c -o $@ $^
mpikernel.o: mpibh.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -c -o $@ $^
mpiinterface.o: mpi_interface.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -c -o $@ $^
hybridkernel.o: hybridbh.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(LIBS) $(OPENMP) -c -o $@ $^

hybridinterface.o: hybrid_interface.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(LIBS) $(OPENMP) -c -o $@ $^

mpibuildinterface.o: interface.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -D$(MPIDEFINE) -c -o $@ $^
ompbuildinterface.o: interface.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) $(OPENMP) -D$(OMPDEFINE) -c -o $@ $^

hybridbuildinterface.o: interface.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(LIBS) $(OPENMP) -D$(HYBRIDDEFINE) -c -o $@ $^

serialinterface.o: interface.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -D$(SERIALDEFINE) -c -o $@ $^

ompbh: utility.o serial.o ompkernel.o ompinterface.o ompbuildinterface.o main.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) $(OPENMP) -D$(OMPDEFINE) -o $@ $^

mpibh: utility.o serial.o mpikernel.o mpiinterface.o mpibuildinterface.o main.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(MPICXX_FLAG) $(LIBS) -D$(MPIDEFINE) -o $@ $^

hybridbh: utility.o serial.o hybridkernel.o hybridinterface.o hybridbuildinterface.o main.cpp
	$(MPICC) -std=$(STD) $(CXX_FLAGS) $(MPICXX_FLAG) $(LIBS) $(OPENMP) -D$(HYBRIDDEFINE) -o $@ $^

serialbh: utility.o serial.o serialinterface.o main.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -D$(SERIALDEFINE) -o $@ $^

visualizer.o: visualizer.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -c -o $@ $^

visualbh: utility.o visualizer.o visualserrial.o serialinterface.o main.cpp
	$(CC) -std=$(STD) $(CXX_FLAGS) $(LIBS) -D$(SERIALDEFINE) -o $@ $^
	rm -rf *.o

visualtest: visualbh
	./visualbh -n 32768 -t 16000

serialtest: serialbh
	./serialbh -n 1024 -t 1000

omptest: ompbh
	export OMP_NUM_THREADS=4
	./serialbh -n 1024 -t 1000

mpitest: mpibh
	mpiexec --oversubscribe -n 4 ./mpibh -n 1024 -t 1000

hybridtest: hybridbh
	export OMP_NUM_THREADS=2
	mpiexec --oversubscribe -n 2 ./hybridbh -n 1024 -t 1000

clean:
	rm -rf *.o serialbh ompbh mpibh hybridbh visualbh images/*