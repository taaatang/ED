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

#include "Global/globalPara.hpp"

#include <stdio.h>
#include <random>
#include <iostream>
#include <type_traits>

#ifdef CPP_17
#include <filesystem>
#else
#include <boost/filesystem.hpp>
#endif //CPP_17

#include <iomanip>
#include <fstream>
#include <list>
#include <cmath>
#include <complex>
#include <mpi.h>
#ifdef OMP_
    #include <omp.h>
#endif

void work_load(idx_t size, int workerID, int workerNum, idx_t& idxStart, idx_t& idxEnd);

/*
    ****************
    * Display Info *
    ****************
*/

void OMP_Info(int workerID);

void mpi_info(int& workerID, int& workerNum);

void exit_msg(std::string msg);

void assert_msg(bool condition, std::string msg);

/******
 * IO *
 ******/

void printModel(LATTICE_MODEL model);

std::string tostr(double val, int digit = 2);

std::string tostr(int val);

std::ostream& operator<<(std::ostream& os, LATTICE_MODEL model);
std::ostream& operator<<(std::ostream& os, ORBITAL orb);
std::ostream& operator<<(std::ostream& os, LINK_TYPE linkt);
std::ostream& operator<<(std::ostream& os, PointGroup p);
std::ostream& operator<<(std::ostream& os, SPIN s);
std::ostream& operator<<(std::ostream& os, LADDER t);

std::ostream& operator<<(std::ostream& os, const VecD& vec);
std::ostream& operator<<(std::ostream& os, const VecI& vec);

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
    for (Type i = 0; i < n1; ++i){
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
    for (Type i = 0; i < n1; ++i){
        result1 *= v1[i];
        result2 *= v2[i];
    }
    delete[] v1;
    delete[] v2;
    return result2/result1;
}

// calculate (-1)^# of permutations -> fermion sign
inline double seqSign(std::vector<int>& seq){
    int counter = 0;
    for(auto it = seq.begin(); it != seq.end(); it++){
        for (auto itp = it + 1; itp != seq.end(); itp++){
            if (*it > *itp) counter++;
        }
    }
    if (counter%2==0) return 1.0;
    return -1.0;
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
    for (int i = 0; i < size; ++i){
        output[i] = mul1 * input1[i] + mul2 * input2[i];
    }
}

inline void vecXAdd(double mul1, const double* input1, double mul2, const double* input2, double* output, int size){
    for (int i = 0; i < size; ++i){
        output[i] = mul1 * input1[i] + mul2 * input2[i];
    }
}

inline void vecXAdd(double mul1, const double* input1, double mul2, const double* input2, double mul3, const double* input3, double* output, int size){
    for (int i = 0; i < size; ++i){
        output[i] = mul1 * input1[i] + mul2 * input2[i] + mul3 * input3[i];
    }
}


template <class T, class U>
inline void vecInit(T *vec, U size, T val){
    #pragma omp parallel for
    for (U i = 0; i < size; ++i) vec[i] = val;
}

template <class T>
inline void veqv(T* vecOut, T* vecIn, idx_t size){
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) vecOut[i] = vecIn[i];
}

template <class T>
inline void veqatv(T* vecOut, T a, T* vecIn, idx_t size){
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) vecOut[i] = a * vecIn[i];
}

template <class T>
inline void saxpy(T* x, T a, T* y, idx_t size){
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) x[i] += a * y[i];
}

template <class T>
inline void vdeva(T* x, T a, idx_t size){
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) x[i] /= a;
}

inline double vdotv(VecD v1, VecD v2){
    assert_msg(v1.size()==v2.size(),"utils::vdotv, v1.size() != v2.size().");
    double result = 0.0;
    #pragma omp parallel for reduction(+:result)
    for(size_t i = 0; i < v1.size(); ++i){
        result += v1[i]*v2[i];
    }
    return result;
}

inline void normalize(VecD& v) {
    vdeva(v.data(),std::sqrt(vdotv(v,v)),v.size());
}

// vector dot prodoct
template <class T1, class T2>
inline void vConjDotv(T1* vconj, T2* v, cdouble* result, idx_t size){
    cdouble partResult = 0.0;
#ifdef OMP_
    cdouble tmp = 0.0;
    double partResultReal = 0.0, partResultImag = 0.0;
    #pragma omp parallel for private(tmp) reduction(+:partResultReal,partResultImag)
#endif
    for (idx_t i = 0; i < size; ++i){
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

/*******
 * MPI *
 *******/

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
/*
    *******
    * I/O *
    *******
*/

#ifdef CPP_17
inline void mkdir_fs(std::string dir){
    std::filesystem::path p(dir);
    std::filesystem::create_directories(p);
    bool succeed = std::filesystem::is_directory(p);
    assert_msg(succeed, dir + " failed to creat!");
}
#else
inline void mkdir_fs(std::string dir){
    boost::filesystem::path p(dir);
    boost::filesystem::create_directories(p);
	bool succeed = boost::filesystem::is_directory(p);
    assert_msg(succeed, dir + " failed to creat!");
}
#endif //CPP_17

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
inline void save(T *d_pt, idx_t size, std::ofstream *f_pt, std::string filename, bool is_app=false){
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
void save(T *d_pt, int size, std::string filename, bool is_app=false){
    std::ofstream os;
    save<T>(d_pt, size, &os, filename, is_app);
}

template <class T>
void save(T *d_pt, idx_t size, std::string filename, bool is_app=false){
    std::ofstream os;
    save<T>(d_pt, size, filename, is_app);
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
inline void read(T *d_pt, idx_t size, std::ifstream *f_pt, std::string filename){
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
        idx_t size = infile.tellg();
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
        idx_t size = infile.tellg();
        infile.seekg(0, infile.beg);
        assert((size % el_size)==0);
        idx_t dim = size/el_size;
        idx_t nlocmax = (dim + workerNum - 1)/workerNum;
        idx_t startRow = workerID * nlocmax;
        idx_t endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
        idx_t nloc = endRow - startRow;
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

template <class T>
inline void infile(std::vector<T*> para, std::string filename){
    std::ifstream file(filename);
    if(file.is_open()){
        for(auto it = para.begin(); it!= para.end(); it++){
            if(!file.eof()){
                file>>*(*it);
            }else{
                std::cout<<filename<<" contains parameters less than expected:"<<para.size()<<"\n";
            }   
        }
    }else{
        std::cout<<filename<<" failed to open!\n";
        exit(1);
    }
    file.close();
}
/*
    ******************
    * Bit Operations *
    * ****************
*/
inline idx_t bitMask(int pos){
    return idx_t{1}<<pos;
}
inline bool bitTest(idx_t n, int pos){
    return (n>>pos) & idx_t{1};
}
inline void bitSet(idx_t& n, int pos){
    n |= bitMask(pos);
}
inline void bitFlip(idx_t& n, int pos){
    n ^= bitMask(pos);
}
inline idx_t bitCount(idx_t& n, const VecI& idxs){
    idx_t sum = 0;
    for(auto it=idxs.begin(); it!=idxs.end(); it++) sum += n>>(*it) & idx_t{1};
    return sum; 
}
inline void bitPrint(idx_t n, int range){
    for(int pos=range-1;pos>=0;pos--)std::cout<<((n>>pos) & idx_t{1});
    std::cout<<std::endl;
}
// From Yao's code
template<class UnsignedType>
inline UnsignedType nextLexicographicalNumber(UnsignedType x) {
    if(x==0)return x+1;
    UnsignedType t = (x | (x - 1)) + 1; //find pivot position and set it to 1
    return t | ((((t & -t) / (x & -x)) >> 1) - 1); //reverse bits after the pivot
}

// push data to an unordered map
template <typename T>
inline void MapPush(MAP<T>* map_pt, idx_t key, T val){
    auto it = map_pt->find(key);
    if constexpr (std::is_same<cdouble, T>::value) {
        if (it == map_pt->end()) {
            (*map_pt)[key] = std::conj(val);
        } else {
            it->second += std::conj(val);
        }
    } else if constexpr (std::is_same<double, T>::value) {
        if (it == map_pt->end()) {
            (*map_pt)[key] = val;
        } else {
            it->second += val;
        }
    } else {
        std::cout<<"MapPush only defined for double/cdouble!\n";
        exit(1);
    }
}

#endif // utils_hpp
