FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: g.o spmat.o Bhat.o main.o checkError.o queue.o
	gcc main.o g.o spmat.o checkError.o Bhat.o queue.o -o cluster $(LIBS)
clean:
	rm -rf *.o cluster

g.o: g.c g.h checkError.h
	gcc $(FLAGS) -c g.c

spmat.o: spmat.c spmat.h checkError.h
	gcc $(FLAGS) -c spmat.c

Bhat.o: Bhat.c Bhat.h g.h checkError.h
	gcc $(FLAGS) -c Bhat.c

queue.o: queue.c queue.h g.h checkError.h
	gcc $(FLAGS) -c queue.c

checkError.o: checkError.c checkError.h
	gcc $(FLAGS) -c checkError.c

main.o: main.c spmat.h Bhat.h g.h checkError.h queue.h
	gcc $(FLAGS) -c main.c