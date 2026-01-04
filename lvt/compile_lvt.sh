source /etc/profile.d/modules.sh  

mpdule purge
module load intel21/compiler-21
module load intel21/mpich-3.4.1

./compile -j 10
