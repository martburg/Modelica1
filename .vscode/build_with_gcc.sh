#!/usr/bin/bash

SRC_UNIX=$(cygpath -u "$1")

gcc -O2 -Wall -shared "$SRC_UNIX" -o ../Resources/Library/spring_solver.dll -llapack -lblas -lm
