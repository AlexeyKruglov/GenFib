genfib: genfib.c convol.c
	gcc $^ -lfftw3 -o $@
