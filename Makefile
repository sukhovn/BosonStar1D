.PHONY: all clean

LIB_DIR = Libraries
COLLOCATION_DIR = Collocation
STAR_DIR = Star
EXTRA_DIR = Extrapolation

LIB_SRC := $(wildcard $(LIB_DIR)/*.cpp)
LIB_BINS := $(patsubst $(LIB_DIR)/%.cpp, obj/libs/%.o, $(LIB_SRC))

COLLOCATION_SRC := $(wildcard $(COLLOCATION_DIR)/*.cpp)
COLLOCATION_BINS := $(patsubst $(COLLOCATION_DIR)/%.cpp, obj/colloc/%.o, $(COLLOCATION_SRC))

STAR_SRC := $(wildcard $(STAR_DIR)/*.cpp)
STAR_BINS := $(patsubst $(STAR_DIR)/%.cpp, obj/star/%.o, $(STAR_SRC))

EXTRA_SRC := $(wildcard $(EXTRA_DIR)/*.cpp)
EXTRA_BINS := $(patsubst $(EXTRA_DIR)/%.cpp, obj/extra/%.o, $(EXTRA_SRC))

HEAD_FOLDERS = $(LIB_DIR) $(COLLOCATION_DIR) $(STAR_DIR) $(EXTRA_DIR)
HEAD_LINK = $(patsubst %, -I%, $(HEAD_FOLDERS))

CC = g++ # compiler to use
CFLAGS = -c -g -O3 -fopenmp -Iobj $(HEAD_LINK) # flags to use at the compliation
LINKERFLAG = -lm -lfftw3 -fopenmp
BINARIES = $(LIB_BINS) $(COLLOCATION_BINS) $(STAR_BINS) $(EXTRA_BINS)

all: prog

obj:
	@echo "Making directory for binaries..."
	mkdir obj
	mkdir obj/libs
	mkdir obj/colloc
	mkdir obj/star
	mkdir obj/extra
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

obj/extra/%.o: $(EXTRA_DIR)/%.cpp obj
	@echo "Compiling extrapolation binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/extrapolator.o: extrapolator.cpp obj
	@echo "Compiling the extrapolator..."
	${CC} ${CFLAGS} -o obj/extrapolator.o extrapolator.cpp

obj/solver.o: solver.cpp obj
	@echo "Compiling the solver..."
	${CC} ${CFLAGS} -o obj/solver.o solver.cpp 

prog: $(BINARIES) obj/solver.o obj/extrapolator.o
	@echo "Linking the $(PROGRAM)..."
	${CC} -o multidomain $(BINARIES) obj/solver.o $(LINKERFLAG)
	${CC} -o extrapolator $(BINARIES) obj/extrapolator.o $(LINKERFLAG)

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -vf	multidomain
	rm -vf	extrapolator