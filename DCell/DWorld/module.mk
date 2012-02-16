# DWorld Module

LOCOBJS := DWorld.o DCell.o DCellsArray.o Euler.o RK2.o BFGS.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${LOCOBJS}: DWorld/DWorld.h

LOCDEPS := DWorld FluidField IIM LevelSetMethod Common
${eval ${call test-library,${subdirectory},dworld,${LOCDEPS}, 2 }}
${eval ${call test-library,${subdirectory},bfgs,${LOCDEPS}, 1 }}

#-ksp_monitor
RUNOPTS =  \
-ksp_atol 1e-2 -ksp_rtol 1e-3 -ksp_max_it 100 \
-fieldsplit_p_ksp_max_it 4 \
-fieldsplit_v_fieldsplit_0_ksp_type preonly \
-fieldsplit_v_fieldsplit_1_ksp_type preonly \
-fieldsplit_v_fieldsplit_0_pc_type cholesky \
-fieldsplit_v_fieldsplit_1_pc_type cholesky \
-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \
-viewer_binary_skip_info \
-info ${PETSC_TMP}/info.log
${eval ${call test-library,${subdirectory},implicit,${LOCDEPS}, 1, ${RUNOPTS} }}
