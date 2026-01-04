source /etc/profile.d/modules.sh

# 모듈 로드
module purge
module load intel21/compiler-21
module load intel21/mpich-3.4.1

# 아키텍처 지정 (LIS에서 인식하는 값)
#export LIS_ARCH=Linux_intel_mpich
export LIS_ARCH=linux_ifc
#i 컴파일러
export LIS_FC=mpifort
export LIS_CC=mpicc
#export LIS_FC=ifort
#export LIS_CC=icc

# 라이브러리 경로 (LDT와 동일하게 최신 사용)
export LIS_JASPER=/usr/local/openjpeg/2.5.0_intel21
export LIS_NETCDF=/usr/local/netcdf/4.9.2_intel21_parallel
export LIS_HDF5=/usr/local/hdf5/1.14.6_parallel
export LIS_HDF4=/usr/local/hdf4/4.3.0_intel21
export LIS_HDFEOS=/usr/local/hdfeos2/3.0_intel21
export LIS_ECCODES=/usr/local
export LIS_MODESMF=/usr/local/esmf/8.7.0_intel21_mpich/mod/modO/Linux.intel.64.mpich.default
export LIS_LIBESMF=/usr/local/esmf/8.7.0_intel21_mpich/lib/libO/Linux.intel.64.mpich.default
export LIS_LAPACK=/usr/local/intel/oneapi/mkl/latest

# LD_LIBRARY_PATH 설정
export LD_LIBRARY_PATH=\
${LIS_HDF4}/lib:${LIS_HDF5}/lib:${LIS_NETCDF}/lib:\
${LIS_HDFEOS}/lib:${LIS_ECCODES}/lib64:${LIS_JASPER}/lib64:\
${LIS_LIBESMF}:${LIS_MODESMF}:${LD_LIBRARY_PATH}

# 컴파일 시작
./configure
#./compile -j 10
