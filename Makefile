CC = clang
CFLAGS = -O3 -march=native -ffast-math -Wall -Wextra
LDFLAGS = -static -lopenblas -lm -flto

dft.out: dft.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

run: dft.out
	@time ./dft.out
	
clean:
	rm -f *.out