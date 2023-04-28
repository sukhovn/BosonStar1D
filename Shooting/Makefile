.PHONY: all clean

CC = g++ # compiler to use
CFLAGS = -c -g -O3 -I. -IRoutines #-Wall -Werror -Wextra -pedantic# flags to use at the compliation
LINKERFLAG = -lm

all: prog

prog: obj/star_shooter.o
	${CC} -o shooter $(LIB_BINS) obj/star_shooter.o ${LINKERFLAG}

obj:
	@echo "Making obj directory..."
	mkdir obj

obj/star_shooter.o: obj star_shooter.cpp
	@echo "Compiling star_shooter.cpp..."
	${CC} ${CFLAGS} -o obj/star_shooter.o star_shooter.cpp

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -rvf Data
	rm -vf shooter
