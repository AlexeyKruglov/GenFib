genfib: genfib.c convol.c
	gcc -O9 $^ -lfftw3 -o $@
