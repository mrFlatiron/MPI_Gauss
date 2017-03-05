EXE = a.out
BUILD_DIR = ./build
BIN_DIR = ./bin
SRC_DIR = ./src
LIB_DIR = ./src/lib
LIBS = matrix utils init_functions sbc
LIBS_DIR = $(patsubst %, $(LIB_DIR)/%, $(LIBS))
OBJS = matrix.o sbc_matrix.o init_functions.o utils.o main.o
INC = -I./include

FLAGS = -Wall -O3 --fast-math -fsanitize=address

all: a.out

clean:
	@cd $(BUILD_DIR); rm -f *.o; cd ..;
	@cd $(BIN_DIR); rm -f $(EXE);
	@cd ..

CC = mpicxx


a.out: $(OBJS)
	cd $(BUILD_DIR); $(CC) $(FLAGS) $(OBJS) -o ../$(BIN_DIR)/a.out

main.o: $(SRC_DIR)/main.cpp
	$(CC) $(FLAGS) $(INC) -c  $(SRC_DIR)/main.cpp -o $(BUILD_DIR)/main.o

matrix.o: $(SRC_DIR)/lib/matrix/matrix.cpp
	$(CC) $(FLAGS) $(INC) -c  $(SRC_DIR)/lib/matrix/matrix.cpp -o $(BUILD_DIR)/matrix.o

utils.o: $(SRC_DIR)/lib/utils/utils.cpp
	$(CC) $(FLAGS) $(INC) -c  $(SRC_DIR)/lib/utils/utils.cpp -o $(BUILD_DIR)/utils.o

sbc_matrix.o: $(SRC_DIR)/lib/sbc/sbc_matrix.cpp
	$(CC) $(FLAGS) $(INC) -c  $(SRC_DIR)/lib/sbc/sbc_matrix.cpp -o $(BUILD_DIR)/sbc_matrix.o

init_functions.o: $(SRC_DIR)/lib/init_functions/init_functions.cpp
	$(CC) $(FLAGS) $(INC) -c $(SRC_DIR)/lib/init_functions/init_functions.cpp -o $(BUILD_DIR)/init_functions.o

