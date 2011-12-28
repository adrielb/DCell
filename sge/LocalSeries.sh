# USAGE
# cluster-fork LocalSeries.sh filelist.txt

#set bash sesssion
# -e exit on error
# -x display command
# -u error on unset variable

set -e -x -u

trap 'echo $SGE_TASK_ID >> $LOG' ERR

#limit the amount of virtual memory consumed to ~1GB
ulimit -v 1000000

LOCAL=$TMP/../temp
mkdir $LOCAL
$SCP $IP:$REMOTE_DIR $LOCAL


rm -rf $LOCAL
