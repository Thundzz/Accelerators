CC=gcc
CFLAGS= -W -Wall -O2 -std=c99
LDFLAGS= -W -Wall -O2 -lm
EXEC=main.out
BENCH=
TEST=

all: $(TEST) $(BENCH) $(EXEC)

test: $(TEST)

send:
	scp * formation:~/accelerateur/
get:
	scp formation:~/accelerateur/datafile .
plot:
	python plot.py

cuda:
	nvcc -deviceemu cuda_nbody.cu -o cuda_nbody

	
main.out: particule.o main.o
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)
	
clean:
	rm -rf *.o  *.[oe][0-9]* *~ *# *.out *.gif datafile* particles*.plt
