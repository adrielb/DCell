

for i in `seq 1 5000`;
do
  export SGE_TASK_ID=$i
  echo "SGE_TASK_ID = " $SGE_TASK_ID
  math < viz.m
done    

