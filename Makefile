# Top level makefile for qe-gipaw

all:    build

build:
	@echo "Building qe-gipaw..."
	$(MAKE) -C src

clean:
	@echo "Cleaning qe-gipaw..."
	if test -s src/Makefile ; then ( $(MAKE) -C src clean ); fi
	-/bin/rm -f bin/gipaw.x

distclean:
	$(MAKE) -C src distclean
	-/bin/rm -f config.log config.status makedeps.sh
	-/bin/rm -Rf autom4te.cache

