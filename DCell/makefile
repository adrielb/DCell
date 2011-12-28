export MODULES=Common LevelSetMethod IIM FluidField DWorld

rebuild:
	@clear
	@make -C Common cleanLib > /dev/null
	@echo "PETSC_ARCH = " ${PETSC_ARCH}
	@START_TIME=`date +%s`; \
	for d in ${MODULES}; do \
		echo Making $$d; \
		${MAKE} -C $$d clean makeLib > /dev/null; \
		echo ; echo; \
	done; \
	END_TIME=`date +%s.%N`; \
	echo "BUILD TIME: " `echo "$$END_TIME - $$START_TIME" | bc` " sec"\

build:
	@for d in ${MODULES}; do \
		make -C $$d makeLib; \
	done;

opt:
	export PETSC_ARCH=gcc-opt; make rebuild



shared: rebuild
	@echo "DCell.so"
	${MAKE} -C Common dcellshared

