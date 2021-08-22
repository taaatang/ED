#pragma once

#include <iostream>
#include <mpi.h>
#include <omp.h>
#include "global/constant.hpp"

inline void ompInfo(int workerID){
    int ompThreadsNum;
    #pragma omp parallel
    {
        #pragma omp master
        ompThreadsNum = omp_get_num_threads();
    }
    if (workerID == MPI_MASTER) std::cout << "openMP turned on with " << ompThreadsNum << " threads" << std::endl;
}

inline void mpiInfo(int& workerID, int& workerNum){
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID == MPI_MASTER) std::cout << "Total MPI Workers: " << workerNum << std::endl;
}

inline void init(int &workerID, int &workerNum, bool &isMaster) {
    MPI_Init(NULL, NULL);
    mpiInfo(workerID, workerNum);
    ompInfo(workerID);
    isMaster = (workerID == MPI_MASTER);
}

/* MPI send */
inline void MPI_Send(int* buf, int count, int dest){
    MPI_Send(buf,count,MPI_INT,dest,0,MPI_COMM_WORLD);
}
/* MPI recv */
inline void MPI_Recv(int* buf, int count, int source){
    MPI_Recv(buf,count,MPI_INT,source,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
/* MPI all gather */
inline void MPI_Allgather(int *vec, int *vec_tot, idx_t nloc) {
    MPI_Allgather(vec, nloc, MPI_INT, vec_tot, nloc, MPI_INT, MPI_COMM_WORLD);
}
inline void MPI_Allgather(long long *vec, long long *vec_tot, idx_t nloc){
    MPI_Allgather(vec, nloc, MPI_LONG_LONG, vec_tot, nloc, MPI_LONG_LONG, MPI_COMM_WORLD);
}
inline void MPI_Allgather(unsigned long long *vec, unsigned long long *vec_tot, idx_t nloc){
    MPI_Allgather(vec, nloc, MPI_UNSIGNED_LONG_LONG, vec_tot, nloc, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
}
inline void MPI_Allgather(double *vec, double *vec_tot, idx_t nloc) {
    MPI_Allgather(vec, nloc, MPI_DOUBLE, vec_tot, nloc, MPI_DOUBLE, MPI_COMM_WORLD);
}
inline void MPI_Allgather(cdouble *vec, cdouble *vec_tot, idx_t nloc) {
    MPI_Allgather(vec, nloc, MPI_DOUBLE_COMPLEX, vec_tot, nloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
}

/* MPI reduce */
inline void MPI_Reduce(int *vec_part, int *vec_tot, int nloc, int pidx) {
    MPI_Reduce(vec_part, vec_tot, nloc, MPI_INT, MPI_SUM, pidx, MPI_COMM_WORLD);
}
inline void MPI_Reduce(double *vec_part, double *vec_tot, int nloc, int pidx) {
    MPI_Reduce(vec_part, vec_tot, nloc, MPI_DOUBLE, MPI_SUM, pidx, MPI_COMM_WORLD);
}
inline void MPI_Reduce(cdouble *vec_part, cdouble *vec_tot, int nloc, int pidx) {
    MPI_Reduce(vec_part, vec_tot, nloc, MPI_DOUBLE_COMPLEX, MPI_SUM, pidx, MPI_COMM_WORLD);
}

/* MPI Allreduce */
inline void MPI_Allreduce(int *part, int* sum, int count){
    MPI_Allreduce(part, sum, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}
inline void MPI_Allreduce(long long *part, long long *sum, int count){
    MPI_Allreduce(part, sum, count, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
}
inline void MPI_Allreduce(unsigned long long *part, unsigned long long *sum, int count){
    MPI_Allreduce(part, sum, count, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
}
inline void MPI_Allreduce(double *part, double *sum, int count){
    MPI_Allreduce(part, sum, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
} 
inline void MPI_Allreduce(cdouble *part, cdouble *sum, int count){
    MPI_Allreduce(part, sum, count, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
} 

/* MPI all to all */
inline void MPI_Alltoall(long long *sendBuff, long long *recBuff,int count){
    MPI_Alltoall(sendBuff,count,MPI_LONG_LONG,recBuff,count,MPI_LONG_LONG,MPI_COMM_WORLD);
}
inline void MPI_Alltoall(unsigned long long *sendBuff, unsigned long long *recBuff,int count){
    MPI_Alltoall(sendBuff,count,MPI_UNSIGNED_LONG_LONG,recBuff,count,MPI_UNSIGNED_LONG_LONG,MPI_COMM_WORLD);
}
inline void MPI_Alltoall(double *sendBuff, double *recBuff,int count){
    MPI_Alltoall(sendBuff,count,MPI_DOUBLE,recBuff,count,MPI_DOUBLE,MPI_COMM_WORLD);
}
inline void MPI_Alltoall(cdouble *sendBuff, cdouble *recBuff,int count){
    MPI_Alltoall(sendBuff,count,MPI_DOUBLE_COMPLEX,recBuff,count,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
}

/* MPI all to all v */
inline void MPI_Alltoallv(long long *sendBuff,int *sendCounts,int *sdispls,long long *recBuff, int *recCounts, int *rdispls){
    MPI_Alltoallv(sendBuff,sendCounts,sdispls,MPI_LONG_LONG,recBuff,recCounts, rdispls,MPI_LONG_LONG,MPI_COMM_WORLD);
}
inline void MPI_Alltoallv(unsigned long long *sendBuff,int *sendCounts,int *sdispls,unsigned long long *recBuff, int *recCounts, int *rdispls){
    MPI_Alltoallv(sendBuff,sendCounts,sdispls,MPI_UNSIGNED_LONG_LONG,recBuff,recCounts, rdispls,MPI_UNSIGNED_LONG_LONG,MPI_COMM_WORLD);
}
inline void MPI_Alltoallv(double *sendBuff,int *sendCounts,int *sdispls,double *recBuff, int *recCounts, int *rdispls){
    MPI_Alltoallv(sendBuff,sendCounts,sdispls,MPI_DOUBLE,recBuff,recCounts, rdispls,MPI_DOUBLE,MPI_COMM_WORLD);
}
inline void MPI_Alltoallv(cdouble *sendBuff,int *sendCounts,int *sdispls,cdouble *recBuff, int *recCounts, int *rdispls){
    MPI_Alltoallv(sendBuff,sendCounts,sdispls,MPI_DOUBLE_COMPLEX,recBuff,recCounts, rdispls,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
}