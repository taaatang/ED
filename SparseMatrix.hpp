//
// SparseMatrix.hpp
// ED
//
// Created by tatang on 10/26/19
//

#ifndef SparseMatrix_hpp
#define SparseMatrix_hpp

#include "globalPara.hpp"
#include "utils.hpp"
#include "mpi.h"
#ifdef OMP_
    #include <omp.h>
#endif

template <class T>
class BaseMatrix{
public:
    int workerID;
    int workerNum;
    ind_int dim;
    ind_int nlocmax;
    ind_int ntot;
    ind_int nloc;
    ind_int startRow;
    ind_int endRow;

    BaseMatrix(ind_int totDim = 0):dim(totDim){
        MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
        MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
        nlocmax = (dim + workerNum - 1)/workerNum;
        ntot = nlocmax * workerNum;
        startRow = workerID * nlocmax;
        endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
        nloc = endRow - startRow;
    }
    ~BaseMatrix(){}

    ind_int get_dim() const {return dim;};
    int get_workerID() const {return workerID;};
    // ntot = nlocmax * workerNum != totDim
    ind_int get_ntot() const {return ntot;};
    ind_int get_nlocmax() const {return nlocmax;};
    // totDim = sum(nloc)
    ind_int get_nloc() const {return nloc;};

    virtual void MxV(T *vecIn, T *vecOut) = 0;
};
/*
    ****************************
    * Sparse Matrix Base Class *
    ****************************
*/

/*
    Sparse Matrix in CSR form (ROW_PARTITION)
    can be easily changed to CSC form for symmetric or hermitian matrix (COL_PARTITION)
    Use the global constant PARTITION in globalPara.hpp to choose the format
    The format will affect how MxV is carried out.
    The buffer vector for MxV in COL_PARTITION setting is scalable with MPI task number.
*/
template <class T>
class SparseMatrix: public BaseMatrix<T>{
/*
  Version 3
*/
protected:
    MATRIX_PARTITION partition;
    int dmNum; // number of diagonal matrix
    int spmNum; // number of sparse matrix in csr format
    std::vector<std::vector<T>> diagValList;
    std::vector<std::vector<T>> valList;
    std::vector<std::vector<ind_int>> colList, rowInitList;
    std::vector<ind_int> counter; // counter[i] is the current number of non-zero elements of i-th sparse matrix
    std::vector<T> diagParameters;
    std::vector<T> parameters;
    static std::vector<T> vecBuf;
    static bool is_vecBuf;
public:
    SparseMatrix(ind_int totDim_ = 0, int spmNum_ = 0, int dmNum_ = 0, MATRIX_PARTITION partition_ = MATRIX_PARTITION::ROW_PARTITION):\
        BaseMatrix<T>(totDim_),partition(partition_),spmNum(spmNum_),dmNum(dmNum_),\
        counter(spmNum_,0),valList(spmNum_),colList(spmNum_),rowInitList(spmNum_),diagValList(dmNum_){
            parameters.resize(spmNum_);
            diagParameters.resize(dmNum_); 
            // for(int i = 0; i < spmNum; i++) rowInitList.at(i).push_back(0);
        };
    ~SparseMatrix(){};

    ind_int nzCount() const {ind_int count=0; for(int i=0;i<spmNum;i++)count+=valList.at(i).size(); return count;}
    
    void reserve(std::vector<T> *pt, ind_int sizePerRow) {pt->reserve(sizePerRow * get_nloc());}
    void reserve(std::vector<ind_int> *pt, ind_int sizePerRow){pt->reserve(sizePerRow * get_nloc());}
    void clear(){
        for (int matID = 0; matID < dmNum; matID++) diagValList[matID].clear();
        for (int matID = 0; matID < spmNum; matID++){valList[matID].clear();colList[matID].clear();rowInitList[matID].clear();}
    }
    void pushRow(std::unordered_map<ind_int,T>* rowMap, int matID=0, MATRIX_PARTITION partition = MATRIX_PARTITION::ROW_PARTITION){
        counter.at(matID) += rowMap->size();
        rowInitList[matID].push_back(counter[matID]);
        for (auto it = rowMap->begin(); it != rowMap->end(); it++){
            colList[matID].push_back(it->first);
            switch(PARTITION){
                case ROW_PARTITION:{
                    valList[matID].push_back(std::conj(it->second));
                    break;
                }
                case COL_PARTITION:{
                    valList[matID].push_back(it->second);
                    break;
                }
            }
        }
    }
    void pushDiag(T val, int matID=0){diagVal.at(matID).push_back(val);}

    void setBuf(){
        switch(partition){
        case ROW_PARTITION:
            vecBuf.resize(BaseMatrix<T>::ntot);
            break;
        case COL_PARTITION:
            vecBuf.resize(BaseMatrix<T>::nlocmax);
            break;
        default:break;
        }
        is_vecBuf = true;
    }
    void clearBuf(){vecBuf.clear(); is_vecBuf=false;}
    void MxV(T *vecIn, T *vecOut);
};

// static members
// buffer for matrix vector multiplication
template <class T>
std::vector<T> SparseMatrix<T>::vecBuf;

// false means resources for vecBuf is not allocated
template <class T>
bool SparseMatrix<T>::is_vecBuf = false; 
/*
    *****************************
    * Sparse Matrix Source Code *
    *****************************
*/
/*
 Version 2. M = a1*M1 + a2*M2 + ...
 */
template <class T>
void SparseMatrix<T>::MxV(T *vecIn, T *vecOut){
    if(!is_vecBuf) setBuf();
    switch(partition){
        // row partition
        case MATRIX_PARTITION::ROW_PARTITION:{
            vecAllGather(vecIn,vecBuf.data(),BaseMatrix<T>::nlocmax);
            // initialize vecOut
            #pragma omp parallel for
            for (ind_int i = 0; i < BaseMatrix<T>::nloc; i++) vecOut[i] = 0.0;
            // constant diagVal
            if (!diagValList.at(0).empty()){
                #pragma omp parallel for
                for (ind_int i = 0; i < BaseMatrix<T>::nloc; i++) vecOut[i] += diagValList[0][i]*vecIn[i];
            }
            // time dependent diagVal
            for (int j = 1; j < dmNum; j++){
                #pragma omp parallel for
                for (ind_int i = 0; i < BaseMatrix<T>::nloc; i++) vecOut[i] += diagParameters[j]*diagValList[j][i]*vecIn[i];
            }
            // constant sparse matrix
            if (!valist.at(0).empty()){
                #pragma omp parallel for
                for (ind_int i = 0; i < BaseMatrix<T>::nloc; i++){
                    for (ind_int j = rowInitList[0].at(i); j < rowInitList[0].at(i+1); j++){
                        vecOut[i] += valList[matID].at(j) * vecBuf[colList[matID].at(j)];
                    }
                }
            }
            // time dependent sparse matrix
            for (int matID = 1; matID < spmNum; matID++){
                #pragma omp parallel for
                for (ind_int i = 0; i < BaseMatrix<T>::nloc; i++){
                    T tmp = 0.0;
                    for (ind_int j = rowInitList[matID].at(i); j < rowInitList[matID].at(i+1); j++){
                        tmp += valList[matID].at(j) * vecBuf[colList[matID].at(j)];
                    }
                    vecOut[i] += parameters.at(matID) * tmp;
                }
            }
            break;
        }
        // column partition need modification !!! 
        case MATRIX_PARTITION::COL_PARTITION:{
            // vecBufTmp = new(std::nothrow) T[nlocmax]; assert(vecBufTmp != NULL);
            for (int id = 0; id < BaseMatrix<T>::workerNum; id++){
                // initialization
                #pragma omp parallel for
                for (ind_int i = 0; i < BaseMatrix<T>::nlocmax; i++){
                    vecBuf[i] = 0.0;
                }
                ind_int rowStart = id * BaseMatrix<T>::nlocmax;
                ind_int rowEnd = (rowStart + BaseMatrix<T>::nlocmax)<BaseMatrix<T>::dim?(rowStart + BaseMatrix<T>::nlocmax):BaseMatrix<T>::dim;

                for (int matID = 0; matID < spmNum; matID++){
                    #pragma omp parallel for
                    for (ind_int col = 0; col < BaseMatrix<T>::nloc; col++){
                        ind_int row;
                        T tmp;
                        #ifdef COMPLEX_OP
                            double *r, *i;
                        #endif
                        for (ind_int rowid = rowInitList[matID].at(col); rowid < rowInitList[matID].at(col+1); rowid++){
                            if ((colList[matID].at(rowid) >= rowStart) && (colList[matID].at(rowid) < rowEnd)){
                                row = colList[matID].at(rowid) - rowStart;
                                tmp = parameters.at(matID)*valList[matID].at(rowid) * vecIn[col];
                                #ifdef COMPLEX_OP
                                    r = &reinterpret_cast<double(&)[2]>(vecBuf[row])[0];
                                    i = &reinterpret_cast<double(&)[2]>(vecBuf[row])[1];
                                    #pragma omp atomic
                                    *r += std::real(tmp);
                                    #pragma omp atomic
                                    *i += std::imag(tmp);
                                #else
                                    #pragma omp atomic
                                    vecBuf[row] += tmp;
                                #endif
                            }
                        }
                    }
                    
                }
                MPI_Barrier(MPI_COMM_WORLD);
                vecReduce(vecBuf.data(), vecOut, rowEnd - rowStart, id);
            }
            break;
        }
        default:{
            std::cout<<"Matrix Vector Multiplication not performed! check the matrix partition format. Only defined for ROW_PARTITION and COL_PARTITION."<<std::endl;
            exit(1);
            break;
        } 
    }
}

#endif
