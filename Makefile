# Top level makefile for qe-gipaw

all:    build

build:
	@echo "Building qe-gipaw..."
	$(MAKE) -C src

clean:
	@echo "Cleaning qe-gipaw..."
	$(MAKE) -C src clean
	-/bin/rm -f bin/gipaw.x

distclean:
	$(MAKE) -C src distclean
	-/bin/rm -f config.log config.status configure makedeps.sh
	-/bin/rm -Rf autom4te.cache

