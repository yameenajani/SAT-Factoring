MROOT := $(shell pwd)
export MROOT

default:
	make -C core rs
	mkdir -p tests/solvers/
	cp -f core/maplesat_static tests/solvers/maplesat_cb

clean:
	make -C core clean
