CC=gcc
CFLAGS= -W -Wall -O2 -std=c99
LDFLAGS= -W -Wall -O2 -lm
EXEC=main.out
BENCH=
TEST=

all: $(TEST) $(BENCH) $(EXEC)

test: $(TEST)

send:
	scp * formation:GPU/

plot:
	gnuplot plot.plt

	
main.out: particule.o main.o
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)
	
clean:
	rm -rf *.o  *.[oe][0-9]* *~ *# *.out *.gif datafile* particles*.plt
