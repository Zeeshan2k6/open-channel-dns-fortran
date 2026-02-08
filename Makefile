FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -std=f2008

SRC = main.f90
OBJ = $(SRC:.f90=.o)

.PHONY: all clean run

all: output open_channel_dns

output:
	mkdir -p output

open_channel_dns: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

run: open_channel_dns
	./open_channel_dns

clean:
	rm -f $(OBJ) open_channel_dns