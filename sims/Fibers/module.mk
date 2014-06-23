# Fiber simulation

RUNOPTS :=                  \
-timax 10000                \
-CFL 0.1                    \
-dtmax 1                    \
-ksp_monitor                \
-viewer_binary_skip_info    \
-info ${PETSC_TMP}/info.log \
-on_error_abort

#-log_summary                                                             \

${eval ${call simulation,Fibers,8,${RUNOPTS}}}
