# Fiber simulation

RUNOPTS :=                  \
-timax 2                    \
-fiber_bend 0.70            \
-fiber_dt 0.5               \
-fluidDrag 10               \
-max_verts 100              \
-num_fibers_per_proc 50     \
-ksp_monitor                \
-viewer_binary_skip_info    \
-info ${PETSC_TMP}/info.log \
-on_error_abort

#-log_summary                                                             \

${eval ${call simulation,Fibers,8,${RUNOPTS}}}
