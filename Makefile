# Top level makefile for ce-gipaw

all:    build

build:
	@echo "Building ce-gipaw..."
	make -C src

clean:
	@echo "Cleaning ce-gipaw..."
	make -C src clean
	-/bin/rm -f bin/gipaw.x

distclean:
	make -C src distclean
	-/bin/rm -f config.log config.status configure makedeps.sh
	-/bin/rm -Rf autom4te.cache

