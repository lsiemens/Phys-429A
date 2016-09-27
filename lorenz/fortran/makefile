compiler = gfortran -Wall -Wextra -Wimplicit-interface -fmax-errors=1 -Wconversion -pedantic -ffree-line-length-none
debug = -g -fcheck=all -fbacktrace
normal = -Werror -O3 -march=native -ffast-math -funroll-loops

lyapunov.exe : find_lyapunov.f90 constants.o
	$(compiler) $(normal) constants.o find_lyapunov.f90 -o lyapunov.exe

constants.o : constants.f90
	$(compiler) $(normal) -c constants.f90 -o constants.o

run : lyapunov.exe
	./lyapunov.exe

new : clean run
	

clean :
	rm -f ./solution.txt
	rm -f ./constants.o
	rm -f ./constants.mod
	rm -f ./lyapunov.exe
