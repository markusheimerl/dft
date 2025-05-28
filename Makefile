CC = clang
CFLAGS = -O3 -march=native -ffast-math -Wall -Wextra
LDFLAGS = -lopenblas -llapacke -lm -lcint -flto

hf.out: hf.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

run: hf.out
	@time ./hf.out
	
clean:
	rm -f *.out