SHELL = /bin/sh

include ../Pablo/makeinclude

DIRS = src

all:
	@for dir in $(DIRS); do\
		echo "******** Making $$dir ********";\
		(cd $$dir; make) || break;\
	done

clean:
	@for dir in $(DIRS); do\
		echo "******** Cleaning $$dir ********";\
		(cd $$dir; make clean) || break;\
	done
