root_dir= $(shell chdir)
src_dir= $(shell  dir  src /D /b)
CPP_C = g++
bin_dir=$(root_dir)\bin
obj_dir=$(root_dir)\obj
SOURCES := $(patsubst %, $(root_dir)/%,$(wildcard src/*.cpp src/*/*.cpp src/*/*/*.cpp))
CPP_FLAGS=-g
vpath %.o $(bin_dir)
vpath %.exe $(root_dir)
export CPP_C root_dir src_dir bin_dir obj_dir CPP_FLAGS SOURCES
run: $(src_dir) C_Main

	
$(src_dir): 
	make -C src/$@

C_Main:
	make -C main
	make -C obj
.PHONEY: MAIN

py:
	make -C python

all: $(src_dir) C_Main py
