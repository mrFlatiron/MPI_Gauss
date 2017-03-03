all: a.out

clean:
	rm -rf *.o
	rm -rf *~

CC = mpicxx

BUILD_DIR = ./build
BIN_DIR = ./bin
SRC_DIR = ./src
OBJS = main.o matrix.o utils.o sbc_matrix.o init_functions.o

FLAGS = -Wall -O3 --fast-math 

a.out: $(OBJS)
	cd $(BUILD_DIR); $(CC) $(FLAGS) sbc_matrix.o main.o matrix.o utils.o  init_functions.o -o ../$(BIN_DIR)/a.out

main.o: $(SRC_DIR)/main.cpp
	$(CC) $(FLAGS) -c  $(SRC_DIR)/main.cpp -o $(BUILD_DIR)/main.o

matrix.o: $(SRC_DIR)/matrix.cpp
	$(CC) $(FLAGS) -c  $(SRC_DIR)/matrix.cpp -o $(BUILD_DIR)/matrix.o

utils.o: $(SRC_DIR)/utils.cpp
	$(CC) $(FLAGS) -c  $(SRC_DIR)/utils.cpp -o $(BUILD_DIR)/utils.o

sbc_matrix.o: $(SRC_DIR)/sbc_matrix.cpp
	$(CC) $(FLAGS) -c  $(SRC_DIR)/sbc_matrix.cpp -o $(BUILD_DIR)/sbc_matrix.o

init_functions.o: $(SRC_DIR)/init_functions.cpp
	$(CC) $(FLAGS) -c $(SRC_DIR)/init_functions.cpp -o $(BUILD_DIR)/init_functions.o

