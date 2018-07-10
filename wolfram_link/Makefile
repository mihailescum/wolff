
CC = clang
WSPREP = /opt/Mathematica/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions/wsprep
CFLAGS = -g -Os -O3 -Wall -fno-strict-aliasing -Wstrict-overflow -Wno-missing-field-initializers -flto -fopenmp=libiomp5 -I/usr/lib/gcc/x86_64-linux-gnu/4.8/include/ -march=native -I/opt/Mathematica/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions/ -I../lib/
LDFLAGS = -L/opt/Mathematica/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions -lm -lpthread -lrt -lstdc++ -ldl -luuid -l WSTP64i4

convex.o: ../lib/convex.c
	@${CC} -c -o $@ $< ${CFLAGS}

convexminorant: convexminorant.tm convex.o
	@${WSPREP} $< -o $@.c
	@${CC} -o $@ convex.o $@.c ${CFLAGS} ${LDFLAGS} 

