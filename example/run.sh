#!/usr/bin/env bash

if [ ! -e ../bin/jinv ]; then
	make clean -C ..
	make -C ..
	make install -C ..
fi

../bin/jinv -f synth/input_mag.in -g synth/input_grv.in -l 0:-1 -s settings.par -v 2>&1 | tee log

