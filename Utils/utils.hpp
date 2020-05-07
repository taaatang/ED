//
//  utils.hpp
//  ED
//
//  Created by tatang on 10/27/19.
//  Copyright © 2019 tatang. All rights reserved.
//

/* utils_hpp */
#ifndef utils_hpp
#define utils_hpp

#include "../Global/globalPara.hpp"

#include <stdio.h>
#include <random>
#include <iostream>
#include <fstream>
#include <list>
#include <cmath>
#include <complex>
#include <mpi.h>
#ifdef OMP_
    #include <omp.h>
#endif

/*
    ****************
    * Display Info *
    ****************
*/

inline void OMP_Info(int workerID){
    #ifdef OMP_
        int ompThreadsNum;
        #pragma omp parallel
        {
            #pragma omp master
            ompThreadsNum = omp_get_num_threads();
        }
        if (workerID==MPI_MASTER) std::cout<<"openMP turned on with "<<ompThreadsNum<<" threads"<<std::endl;
    #else
        if (workerID==MPI_MASTER) std::cout<<"openMP turned off"<<std::endl;
    #endif
}

inline void mpi_info(int& workerID, int& workerNum){
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID==MPI_MASTER) std::cout<<"Total MPI Workers:"<<workerNum<<std::endl;
    OMP_Info(workerID);
}

void work_load(ind_int size, int workerID, int workerNum, ind_int& idxStart, ind_int& idxEnd){
    assert((workerNum>0) && (workerID<workerNum));
    ind_int nlocmax = (size + workerNum - 1)/workerNum;
    idxStart = workerID * nlocmax;
    idxEnd = (idxStart + nlocmax)<size?(idxStart + nlocmax):size;
}

void errExit(std::string msg);

inline void assert_msg(bool condition, std::string msg){
    if(!condition){
        std::cout<<msg<<std::endl;
        exit(1);
    }
}

/*
    ***************
    * Combinatory *
    ***************
*/
unsigned long long pow(unsigned long long base, unsigned long long power);

long long pow(long long base, long long power);

template <class Type>
Type factorial(Type n) {
    if (n < 0){std::cerr<<"n must be a non-negative integer!"; exit(EXIT_FAILURE);}
    Type f = 1;
    while (n > 0){
        f *= n;
        n--;
    }
    return f;
}

template <class Type>
Type combination(Type n, Type k){
    if (n < k or n < 1) {std::cerr<<"n can't be smaller than k and both n,k should be positive integer!"<<std::endl; exit(EXIT_FAILURE);}
    if (n == k or k==0) return 1;
    Type n1 = k<(n-k)?k:n-k;
    Type* v1 = new Type[n1];
    Type* v2 = new Type[n1];
    for (Type i = 0; i < n1; i++){
        v1[i] = i + 1;
        v2[i] = n - n1 + i + 1;
    }
    for (Type i = n1; i > 0; i--){
        for (Type j = 0; j < n1; j++){
            if (v2[j] % v1[i-1] == 0){
                v2[j] /= v1[i-1];
                v1[i-1] = 1;
                break;
            }
        }
    }
    Type result1 = 1;
    Type result2 = 1;
    for (Type i = 0; i < n1; i++){
        result1 *= v1[i];
        result2 *= v2[i];
    }
    delete[] v1;
    delete[] v2;
    return result2/result1;
}

// calculate (-1)^# of permutations -> fermion sign
inline int seqSign(std::vector<int>& seq){
    int counter = 0;
    for(auto it = seq.begin(); it != seq.end(); it++){
        for (auto itp = it + 1; itp != seq.end(); itp++){
            if (*it > *itp) counter++;
        }
    }
    if (counter%2==0) return 1;
    return -1;
}

// Gaussian Pulse
inline double GaussPulse(double t, double amp, double sigma, double freq){
    return amp*std::exp(-t*t/sigma/sigma)*std::cos(freq*t);
}
/*
    ***********
    * Algebra *
    ***********
*/
template <class T>
inline void randInit(std::vector<T>& vec, int range=100){
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,range);
    T ranget = range;
    for (auto it = vec.begin(); it!=vec.end();it++) *it = T(dist(rng))/ranget;
}
inline void vecXAdd(int mul1, const int* input1, int mul2, const int* input2, int* output, int size){
    for (int i = 0; i < size; i++){
        output[i] = mul1 * input1[i] + mul2 * input2[i];
    }
}

inline void vecXAdd(double mul1, const double* input1, double mul2, const double* input2, double* output, int size){
    for (int i = 0; i < size; i++){
        output[i] = mul1 * input1[i] + mul2 * input2[i];
    }
}

inline void vecXAdd(double mul1, const double* input1, double mul2, const double* input2, double mul3, const double* input3, double* output, int size){
    for (int i = 0; i < size; i++){
        output[i] = mul1 * input1[i] + mul2 * input2[i] + mul3 * input3[i];
    }
}


template <class T, class U>
inline void vecInit(T *vec, U size, T val){
#ifdef OMP_
    #pragma omp parallel for
#endif
    for (U i = 0; i < size; i++) vec[i] = val;
}

template <class T>
inline void veqv(T* vecOut, T* vecIn, ind_int size){
    #ifdef OMP_
    #pragma omp parallel for
    #endif
    for (ind_int i = 0; i < size; i++) vecOut[i] = vecIn[i];
}

template <class T>
inline void veqatv(T* vecOut, T a, T* vecIn, ind_int size){
    #ifdef OMP_
    #pragma omp parallel for
    #endif
    for (ind_int i = 0; i < size; i++) vecOut[i] = a * vecIn[i];
}

template <class T>
inline void saxpy(T* x, T a, T* y, ind_int size){
    #ifdef OMP_
    #pragma omp parallel for
    #endif
    for (ind_int i = 0; i < size; i++) x[i] += a * y[i];
}

template <class T>
inline void vdeva(T* x, T a, ind_int size){
    #ifdef OMP_
    #pragma omp parallel for
    #endif
    for (ind_int i = 0; i < size; i++) x[i] /= a;
}

// vector dot prodoct
template <class T1, class T2>
inline void vConjDotv(T1* vconj, T2* v, cdouble* result, ind_int size){
    cdouble partResult = 0.0;
#ifdef OMP_
    cdouble tmp = 0.0;
    double partResultReal = 0.0, partResultImag = 0.0;
    #pragma omp parallel for private(tmp) reduction(+:partResultReal,partResultImag)
#endif
    for (ind_int i = 0; i < size; i++){
#ifdef OMP_
        tmp = std::conj(vconj[i]) * v[i];
        partResultReal += std::real(tmp);
        partResultImag += std::imag(tmp);
#else
        partResult += std::conj(vconj[i]) * v[i];
#endif
    }
#ifdef OMP_
    partResult = partResultReal + CPLX_I * partResultImag;
#endif
    MPI_Allreduce(&partResult, result, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
}

/* MPI send */
inline void MPI_Send(int* buf, int count, int dest){
    MPI_Send(buf,count,MPI_INT,dest,MPI_ANY_TAG,MPI_COMM_WORLD);
}
/* MPI recv */
inline void MPI_Recv(int* buf, int count, int source){
    MPI_Recv(buf,count,MPI_INT,source,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
/* MPI all gather */
inline void MPI_Allgather(int *vec, int *vec_tot, ind_int nloc) {
    MPI_Allgather(vec, nloc, MPI_INT, vec_tot, nloc, MPI_INT, MPI_COMM_WORLD);
}
inline void MPI_Allgather(long long *vec, long long *vec_tot, ind_int nloc){
    MPI_Allgather(vec, nloc, MPI_LONG_LONG, vec_tot, nloc, MPI_LONG_LONG, MPI_COMM_WORLD);
}
inline void MPI_Allgather(unsigned long long *vec, unsigned long long *vec_tot, ind_int nloc){
    MPI_Allgather(vec, nloc, MPI_UNSIGNED_LONG_LONG, vec_tot, nloc, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
}
inline void MPI_Allgather(double *vec, double *vec_tot, ind_int nloc) {
    MPI_Allgather(vec, nloc, MPI_DOUBLE, vec_tot, nloc, MPI_DOUBLE, MPI_COMM_WORLD);
}
inline void MPI_Allgather(cdouble *vec, cdouble *vec_tot, ind_int nloc) {
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
/*
    *******
    * I/O *
    *******
*/

template <class T>
inline void save(T *d_pt, int size, std::ofstream *f_pt, std::string filename, bool is_app=false){
    if(is_app)f_pt->open(filename, std::ios::binary|std::ios::app);
    else f_pt->open(filename, std::ios::binary);
    if (f_pt->is_open()){
        f_pt->write(reinterpret_cast<char*>(d_pt), size * sizeof(T));
        f_pt->close();
        if(is_app)std::cout<<"Data appended to "<<filename<<std::endl;
        else std::cout<<"Data wrote to "<<filename<<std::endl;
    }else{
        std::cout<<filename<<" failed to open!"<<std::endl;
        exit(1);
    }
}

template <class T>
inline void save(T *d_pt, ind_int size, std::ofstream *f_pt, std::string filename, bool is_app=false){
    if(is_app)f_pt->open(filename, std::ios::binary|std::ios::app);
    else f_pt->open(filename, std::ios::binary);
    if (f_pt->is_open()){
        f_pt->write(reinterpret_cast<char*>(d_pt), size * sizeof(T));
        f_pt->close();
        if(is_app)std::cout<<"Data appended to "<<filename<<std::endl;
        else std::cout<<"Data wrote to "<<filename<<std::endl;
    }else{
        std::cout<<filename<<" failed to open!"<<std::endl;
        exit(1);
    }
}



template <class T>
inline void read(T *d_pt, int size, std::ifstream *f_pt, std::string filename){
    f_pt->open(filename, std::ios::binary);
    if (f_pt->is_open()){
        f_pt->read(reinterpret_cast<char*>(d_pt), size * sizeof(T));
        f_pt->close();
        // std::cout<<"Data loaded from "<<filename<<std::endl;
    }else{
        std::cout<<filename<<" failed to open!"<<std::endl;
        exit(1);
    }
}

template <class T>
inline void read(T *d_pt, ind_int size, std::ifstream *f_pt, std::string filename){
    f_pt->open(filename, std::ios::binary);
    if (f_pt->is_open()){
        f_pt->read(reinterpret_cast<char*>(d_pt), size * sizeof(T));
        f_pt->close();
        // std::cout<<"Data loaded from "<<filename<<std::endl;
    }else{
        std::cout<<filename<<" failed to open!"<<std::endl;
        exit(1);
    }
}

template <class T>
inline void read(std::vector<T> *d_pt, std::string filename){
    std::ifstream infile;
    infile.open(filename, std::ios::in|std::ios::binary);
    if (infile.is_open()){
        infile.seekg(0, infile.end);
        ind_int size = infile.tellg();
        infile.seekg(0, infile.beg);
        assert((size % sizeof(T))==0);
        d_pt->resize(size/sizeof(T));
        infile.read(reinterpret_cast<char*>(d_pt->data()), size);
        infile.close();
        // std::cout<<"Data loaded from "<<filename<<std::endl;
    }else{
        std::cout<<filename<<" failed to open!"<<std::endl;
        exit(1);
    }
}

template <class T>
inline void read(std::vector<T> *d_pt, std::string filename, int workerID, int workerNum){
    std::ifstream infile;
    infile.open(filename, std::ios::in|std::ios::binary);
    if (infile.is_open()){
        int el_size = sizeof(T);
        infile.seekg(0, infile.end);
        ind_int size = infile.tellg();
        infile.seekg(0, infile.beg);
        assert((size % el_size)==0);
        ind_int dim = size/el_size;
        ind_int nlocmax = (dim + workerNum - 1)/workerNum;
        ind_int startRow = workerID * nlocmax;
        ind_int endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
        ind_int nloc = endRow - startRow;
        infile.seekg(startRow*el_size, infile.beg);
        d_pt->resize(nloc);
        infile.read(reinterpret_cast<char*>(d_pt->data()), nloc*el_size);
        infile.close();
        // std::cout<<"Data loaded from "<<filename<<std::endl;
    }else{
        std::cout<<filename<<" failed to open!"<<std::endl;
        exit(1);
    }
}
/*
    ******************
    * Bit Operations *
    * ****************
*/
inline ind_int bitMask(int pos){
    return ind_int{1}<<pos;
}
inline bool bitTest(ind_int n, int pos){
    return (n>>pos) & ind_int{1};
}
inline void bitSet(ind_int& n, int pos){
    n |= bitMask(pos);
}
inline void bitFlip(ind_int& n, int pos){
    n ^= bitMask(pos);
}
inline ind_int bitCount(ind_int& n, const VecI& idxs){
    ind_int sum = 0;
    for(auto it=idxs.begin(); it!=idxs.end(); it++) sum += n>>(*it) & ind_int{1};
    return sum; 
}
inline void bitPrint(ind_int n, int range){
    for(int pos=range-1;pos>=0;pos--)std::cout<<((n>>pos) & ind_int{1});
    std::cout<<std::endl;
}
// From Yao's code
template<class UnsignedType>
inline UnsignedType nextLexicographicalNumber(UnsignedType x) {
    UnsignedType t = (x | (x - 1)) + 1; //find pivot position and set it to 1
    return t | ((((t & -t) / (x & -x)) >> 1) - 1); //reverse bits after the pivot
}

// push data to an unordered map
inline void MapPush(MAP* map_pt, ind_int key, dataType val){
    auto it = map_pt->find(key);
    if (it == map_pt->end()){
        (*map_pt)[key] = val;
    }else{
        it->second += val;
    }
}

#endif // utils_hpp
