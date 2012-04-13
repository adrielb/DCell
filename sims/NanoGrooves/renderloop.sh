
export SGE_TASK_ID = 0
while true; do
  sleep 0.1
  echo $SGE_TASK_ID 
  math < viz.m
  export SGE_TASK_ID=$[${SGE_TASK_ID} + 1 ]
done

#  math < viz.m    
