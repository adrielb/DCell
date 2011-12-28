set -v

# returns device number
grep aac /proc/devices

# make sure device number is the same one given to mknod
cd /dev
mknod afa0 c 253 0
ls -l afa0
