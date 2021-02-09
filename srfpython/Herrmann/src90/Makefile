# Make file, will be executed by python setuptools

.DEFAULT_GOAL := all
.PHONY : clean compile test all

clean:
	bash clean.sh;

compile:
	bash compile.sh;

test:
	bash ./test.sh

all: clean compile test
