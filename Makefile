all: a.out

clean:
	rm -rf *.o
	rm -rf *~

CC = mpicxx

FLAGS = -Wall 

a.out: main.o matrix.o utils.o sbc_matrix.o init_functions.o
	$(CC) $(FLAGS) sbc_matrix.o main.o matrix.o utils.o  init_functions.o -o a.out

main.o: main.cpp
	$(CC) $(FLAGS) -c main.cpp

matrix.o: matrix.cpp
	$(CC) $(FLAGS) -c matrix.cpp

utils.o: utils.cpp
	$(CC) $(FLAGS) -c utils.cpp

sbc_matrix.o: sbc_matrix.cpp
	$(CC) $(FLAGS) -c sbc_matrix.cpp

init_functions.o: init_functions.cpp
	$(CC) $(FLAGS) -c init_functions.cpp
