assn3: nBody.c submit.c 
	mpicc -Wall -m64 -DBUILD_64 -O3 -DNDEBUG  -o assn3 nBody.c submit.c -lm

clean: 
	rm *~
	rm assn3
