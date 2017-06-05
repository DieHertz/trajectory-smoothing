PREFIX := /usr/local/Cellar
INCLUDE_DIRS := $(PREFIX)/nlopt/2.4.2_2/include \
				$(PREFIX)/gsl/2.3/include
LIB_DIRS := $(PREFIX)/nlopt/2.4.2_2/lib \
			$(PREFIX)/gsl/2.3/lib

default: quintic.cpp
	g++ quintic.cpp -o quintic $(addprefix -I,$(INCLUDE_DIRS)) $(addprefix -L,$(LIB_DIRS)) -std=c++1z -lm -lnlopt -lgsl -O3 -Ofast

cubic: main.cpp
	g++ main.cpp -o cubic $(addprefix -I,$(INCLUDE_DIRS)) $(addprefix -L,$(LIB_DIRS)) -std=c++1z -lm -lnlopt -lgsl -O3 -Ofast