.PHONY: all clean

LIB_DIR = Libraries
COLLOCATION_DIR = Collocation
STAR_DIR = Star

LIB_SRC := $(wildcard $(LIB_DIR)/*.cpp)
LIB_BINS := $(patsubst $(LIB_DIR)/%.cpp, obj/libs/%.o, $(LIB_SRC))

COLLOCATION_SRC := $(wildcard $(COLLOCATION_DIR)/*.cpp)
COLLOCATION_BINS := $(patsubst $(COLLOCATION_DIR)/%.cpp, obj/colloc/%.o, $(COLLOCATION_SRC))

STAR_SRC := $(wildcard $(STAR_DIR)/*.cpp)
STAR_BINS := $(patsubst $(STAR_DIR)/%.cpp, obj/star/%.o, $(STAR_SRC))

HEAD_FOLDERS = $(LIB_DIR) $(COLLOCATION_DIR) $(STAR_DIR)
HEAD_LINK = $(patsubst %, -I%, $(HEAD_FOLDERS))

CC = g++ # compiler to use
CFLAGS = -c -O3 -fopenmp -Iobj $(HEAD_LINK) # flags to use at the compliation
LINKERFLAG = -lm -lfftw3 -fopenmp
BINARIES = $(LIB_BINS) $(COLLOCATION_BINS) $(STAR_BINS)

all: prog

obj:
	@echo "Making directory for binaries..."
	mkdir obj
	mkdir obj/libs
	mkdir obj/colloc
	mkdir obj/star
	xxd -i ./$(STAR_DIR)/par_solitonic.dat ./obj/par_solitonic.h 
	xxd -i ./$(STAR_DIR)/par_polynomial.dat ./obj/par_polynomial.h 

obj/libs/%.o: $(LIB_DIR)/%.cpp obj
	@echo "Compiling library binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/colloc/%.o: $(COLLOCATION_DIR)/%.cpp obj
	@echo "Compiling collocation binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/star/%.o: $(STAR_DIR)/%.cpp obj
	@echo "Compiling star binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/solver.o: solver.cpp obj
	@echo "Compiling the solver..."
	${CC} ${CFLAGS} -o obj/solver.o solver.cpp 

prog: $(BINARIES) obj/solver.o
	@echo "Linking the $(PROGRAM)..."
	${CC} -o multidomain $(BINARIES) obj/solver.o $(LINKERFLAG)

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -vf	multidomain
	rm -vf	extrapolator