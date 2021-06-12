//
// SparseMatrix.hpp
// ED
//
// Created by tatang on 10/26/19
//

#ifndef SparseMatrix_hpp
#define SparseMatrix_hpp

#include <type_traits>
#include <climits>
#include <omp.h>

#include "Global/globalPara.hpp"
#include "Basis/Basis.hpp"
#include "Utils/mpiwrap.hpp"
#include "Algebra/algebra.hpp"
#include "Utils/timer.hpp"

template <class T>
class BaseMatrix{
public:
    BaseMatrix(idx_t totRowDim = 0, idx_t totColDim = 0);
    ~BaseMatrix( ) { }

    int getWorkerID( ) const { return workerID; }
    
    void setDim(idx_t totDim = 0);
    idx_t getDim( ) const { return dim; }
    idx_t getntot( ) const { return ntot; }
    idx_t getnlocmax( ) const { return nlocmax; }
    idx_t getnloc( ) const { return nloc; }

    void setColDim(idx_t totDim = 0);
    idx_t getColDim( ) const { return dim_col; }
    idx_t getColntot( ) const { return ntot_col; }
    idx_t getColnlocmax( ) const { return nlocmax_col; }
    idx_t getColnloc( ) const { return nloc_col; }

    // Matrix Vector Multiplication Interface
    virtual void MxV(T *vecIn, T *vecOut) = 0;
    virtual void getDiag(std::vector<T> &diag) {}
    // project to symmetry subspace labeled by val
    virtual void project(double val, T* vec) { };

    int workerID;
    int workerNum;
    // row major. row dimension
    // tot row dim of all workers
    // totDim = sum(nloc)
    idx_t dim;
    idx_t nlocmax;
    // ntot = nlocmax * workerNum != totDim
    idx_t ntot;
    idx_t nloc;
    idx_t startRow;
    idx_t endRow;

    idx_t dim_col;
    idx_t nlocmax_col;
    idx_t ntot_col;
    idx_t nloc_col;
    idx_t startCol;
    idx_t endCol;
};

template <class T>
class Identity: public BaseMatrix<T> {
public:
    Identity(idx_t totDim): BaseMatrix<T>(totDim, totDim) { };
    void MxV(T* vecIn, T* vecOut) {
        copy(BaseMatrix<T>::nloc, vecIn, vecOut);
    }
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
  Version 4
*/
public:
    SparseMatrix( ) { }
    SparseMatrix(Basis *Bi_, Basis *Bf_, int spmNum_ = 0, int dmNum_ = 0);
    ~SparseMatrix( );

    // local matrix nonzero elements count
    idx_t nzCount( ) const;    

    // reserve memory for diagonal part
    void reserveDiag( ) { for(int i = 0; i < dmNum; ++i) diagValList.reserve(this->getnloc()); }

    // reserve memory for off diagonal part
    void reserve(idx_t sizePerRow, int matID = 0);

    void created( ) { isMatrixFree = false; }

    // clear memory of sparse matrixes
    void clear();

    void setSpmNum(int num);

    void setDmNum(int num);

    void pushDiag(T val, int matID=0) { diagValList.at(matID).push_back(val); }

    void putDiag(T val, idx_t idx, int matID = 0);

    // check if matrix vector multiplication buff is properly set
    bool isMVBuf( );

    // set buffer for MxV
    void setBuf(idx_t size) { vecBuf.resize(size); }

    void setBuf( );

    void clearBuf( ) { vecBuf.clear(); }

    void setVal(int matID, T val) { parameters.at(matID) = val; }

    void MxV(T *vecIn, T *vecOut);

    void getDiag(std::vector<T> &diag);

    T vMv(T *vecL, T *vecR);

    virtual void print(std::string info, std::ostream& os = std::cout) const;

    // create one row. store sparse matrixes data in corresponding rowMaps. (repI, val) for distributed basis. (idx,val) otherwise
    virtual void row(idx_t rowID, std::vector<MAP<T>>& rowMaps){};


#ifdef DISTRIBUTED_BASIS

    // initialize sendBuff/recvBuff for distributed basis MPI_Alltoall communication
    void setMpiBuff(idx_t idx_val);

    // construct sparse matrix in parallel. each thread create #rowPerThread.
    virtual void construct(int rowCount=50, int rowPerIt=1000);

#else // DISTRIBUTED_BASIS

    void pushRow(std::unordered_map<idx_t,T>* rowMap, int matID = 0);

    void build( ) {
        for (int i = 0; i < spmNum; ++i) {
            MKL::create(A.at(i), BaseMatrix<T>::nloc, Bi->getSubDim(), rowInitList.at(i), colList.at(i), valList.at(i));
        }
        isMatrixFree = false;
    }

    virtual void construct(int rowPerThread = 1);

#endif // DISTRIBUTED_BASIS

protected:

    MATRIX_PARTITION partition{ROW_PARTITION};
    // <Bf|Op|Bi>.
    Basis *Bi{nullptr}, *Bf{nullptr};
    int dmNum{0}; // number of diagonal matrix
    int spmNum{0}; // number of sparse matrix in csr format
    bool isMatrixFree{true};

    std::vector<T> diagParameters; // diagParameters[i]*diagValList[i]. parameter for i-th diagonal matrix
    std::vector<T> parameters; // parameters[i]*valList[i]. parameter for i-th sparse matrix

    static std::vector<T> vecBuf; // common buffer for matrix vector multiplication

#ifdef DISTRIBUTED_STATE
    int blockNum;
    std::vector<std::vector<T>> diagValList;
    //[matrix_id, row_block_id, values]
    std::vector<std::vector<std::vector<T>>> valList;
    std::vector<std::vector<std::vector<idx_t>>> colList, rowInitList;

    std::vector<std::vector<idx_t>> idxSendBuff, idxRecvBuff;
    std::vector<std::vector<T>> valSendBuff, valRecvBuff;

    std::vector<std::vector<idx_t>> counter; // counter[i,j] is the current number of non-zero elements of i-th sparse matrix's j-th block
    // [matrix_id,row_block_id]
    std::vector<std::vector<sparse_matrix_t>> A;

#else

    std::vector<std::vector<T>> diagValList;
    std::vector<std::vector<T>> valList;
    std::vector<std::vector<idx_t>> colList, rowInitList;
    std::vector<idx_t> counter; // counter[i] is the current number of non-zero elements of i-th sparse matrix
    std::vector<sparse_matrix_t> A;

#endif

};

template <class T>
BaseMatrix<T>::BaseMatrix(idx_t totRowDim, idx_t totColDim) {
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    setDim(totRowDim);
    setColDim(totColDim);
}

template <class T>
void BaseMatrix<T>::setDim(idx_t totDim) {
    dim = totDim;
    nlocmax = (dim + workerNum - 1)/workerNum;
    ntot = nlocmax * workerNum;
    startRow = workerID * nlocmax;
    endRow = (startRow + nlocmax) < dim ? (startRow + nlocmax) : dim;
    nloc = endRow - startRow;
}

template <class T>
void BaseMatrix<T>::setColDim(idx_t totDim) {
    dim_col = totDim;
    nlocmax_col = (dim_col + workerNum - 1)/workerNum;
    ntot_col = nlocmax_col * workerNum;
    startCol = workerID * nlocmax_col;
    endCol = (startCol + nlocmax_col) < dim_col ? (startCol + nlocmax_col) : dim_col;
    nloc_col = endCol - startCol;
}

// buffer for matrix vector multiplication
template <class T>
std::vector<T> SparseMatrix<T>::vecBuf;

/*
    *****************************
    * Sparse Matrix Source Code *
    *****************************
*/

template <class T>
SparseMatrix<T>::SparseMatrix(Basis *Bi_, Basis *Bf_, int spmNum_ , int dmNum_):BaseMatrix<T>(Bf_->getSubDim(), Bi_->getSubDim()), \
Bi(Bi_), Bf(Bf_) {

    setSpmNum(spmNum_);
    setDmNum(dmNum_);

    #ifdef DISTRIBUTED_BASIS

        partition =  MATRIX_PARTITION::COL_PARTITION;

    #else

        partition = MATRIX_PARTITION::ROW_PARTITION;   

    #endif
}

template <class T>
SparseMatrix<T>::~SparseMatrix( ) {
    for(int i = 0; i < spmNum; ++i) {
        #ifdef DISTRIBUTED_BASIS

            for(int b = 0; b < BaseMatrix<T>::workerNum;b++) MKL::destroy(A.at(i).at(b));

        #else

            MKL::destroy(A.at(i));

        #endif
    }
}

template <class T>
idx_t SparseMatrix<T>::nzCount( ) const {
    idx_t count = 0;

    #ifdef DISTRIBUTED_BASIS 

        for(int i=0;i<spmNum;++i)for(int b=0;b<BaseMatrix<T>::workerNum;b++)count+=valList.at(i).at(b).size(); 

    #else

        for(int i = 0; i < spmNum; ++i) count += valList.at(i).size();

    #endif

    for(int i = 0; i < dmNum; ++i) count += diagValList.at(i).size();
    return count;
}

template <class T>
void SparseMatrix<T>::reserve(idx_t sizePerRow, int matID) {
    #ifdef DISTRIBUTED_BASIS

        for(int b=0; b<blockNum; b++){
            valList.at(matID).at(b).reserve(sizePerRow * BaseMatrix<T>::nloc); 
            colList.at(matID).at(b).reserve(sizePerRow * BaseMatrix<T>::nloc);
            rowInitList.at(matID).at(b).reserve(1+BaseMatrix<T>::nloc);
        }

    #else

        valList.at(matID).reserve(sizePerRow * BaseMatrix<T>::nloc); 
        colList.at(matID).reserve(sizePerRow * BaseMatrix<T>::nloc);
        rowInitList.at(matID).reserve(1+BaseMatrix<T>::nloc);

    #endif
}

template <class T>
void SparseMatrix<T>::clear( ) {
    for (int matID = 0; matID < dmNum; ++matID) diagValList[matID].clear();

    #ifdef DISTRIBUTED_BASIS

        for (int matID = 0; matID < spmNum; ++matID) {
            for(int b = 0; b < blockNum; ++b) {
                counter.at(matID).at(b) = 0;
                valList.at(matID).at(b).clear();
                colList.at(matID).at(b).clear();
                rowInitList.at(matID).at(b).clear();
            }
        }

    #else

        for (int matID = 0; matID < spmNum; ++matID) {
            counter[matID] = 0;
            valList[matID].clear();
            colList[matID].clear();
            rowInitList[matID].clear();
        }
    
    #endif
    isMatrixFree = true;
}

template <class T>
void SparseMatrix<T>::putDiag(T val, idx_t idx, int matID){

    #ifdef DISTRIBUTED_BASIS

        diagValList.at(matID).at(idx) = val;

    #else

        diagValList.at(matID).at(idx-BaseMatrix<T>::startRow) = val;

    #endif

}

template <class T>
void SparseMatrix<T>::setSpmNum(int num) {
    spmNum = num; 
    parameters = std::vector<T>(spmNum,1.0);
    counter.resize(spmNum);
    valList.resize(spmNum);
    colList.resize(spmNum);
    rowInitList.resize(spmNum);
    A.resize(spmNum);

    for (auto &c : counter) {
        c = 0;
    }

    for (auto &row : rowInitList) {
        row.push_back(0);
    }

    #ifdef DISTRIBUTED_STATE

        blockNum = BaseMatrix<T>::workerNum;
        idxSendBuff.resize(num);
        idxRecvBuff.resize(num);
        valSendBuff.resize(num);
        valRecvBuff.resize(num);
        for(int matID=0; matID<spmNum; matID++) {
            counter.at(matID).resize(blockNum,0);
            rowInitList.at(matID).resize(blockNum);
            colList.at(matID).resize(blockNum);
            valList.at(matID).resize(blockNum);
            A.at(matID).resize(blockNum);
        }

    #endif
}

template <class T>
void SparseMatrix<T>::setDmNum(int num) {
    dmNum = num; 
    diagParameters = std::vector<T>(num,1.0);
    diagValList.resize(num);
    for(int i = 0; i < num; ++i) diagValList[i].resize(BaseMatrix<T>::nloc); 
}

#ifndef DISTRIBUTED_BASIS

    template <class T>
    void SparseMatrix<T>::pushRow(std::unordered_map<idx_t,T>* rowMap, int matID){
        counter.at(matID) += rowMap->size();
        rowInitList[matID].push_back(counter[matID]);
        for (auto it = rowMap->begin(); it != rowMap->end(); it++){
            colList[matID].push_back(it->first);
            if constexpr (std::is_same<cdouble, T>::value) {
                valList[matID].push_back(std::conj(it->second));
            } else if constexpr (std::is_same<double, T>::value) {
                valList[matID].push_back(it->second); 
            } else {
                std::cout<<"Sparse matrix only defined for couble/cdouble! \n";
                exit(1);
            }
        }
    }

#endif

template <class T>
bool SparseMatrix<T>::isMVBuf( ) {    
    // nlocmax_col = (Bi->getSubDim() + BaseMatrix<T>::workerNum - 1)/BaseMatrix<T>::workerNum;
    // ntot_col = nlocmax_col * BaseMatrix<T>::workerNum;

    #ifdef DISTRIBUTED_BASIS

        return vecBuf.size() == this->getColnlocmax();

    #else
 
        return vecBuf.size() == this->getColntot();

    #endif
}

template <class T>
void SparseMatrix<T>::setBuf( ) {
    // nlocmax_col = (Bi->getSubDim() + BaseMatrix<T>::workerNum - 1) / BaseMatrix<T>::workerNum;
    // ntot_col = nlocmax_col * BaseMatrix<T>::workerNum;
    #ifdef DISTRIBUTED_BASIS
        vecBuf.resize(this->getColnlocmax());
    #else
        vecBuf.resize(this->getColntot());
    #endif
}

#ifdef DISTRIBUTED_BASIS
template <class T>
void SparseMatrix<T>::setMpiBuff(idx_t idx_val){
    for(auto vecit = idxSendBuff.begin(); vecit != idxSendBuff.end(); vecit++){
        for(auto it = (*vecit).begin(); it != (*vecit).end(); it++) *it = idx_val;
    }
    for(auto vecit = idxRecvBuff.begin(); vecit != idxRecvBuff.end(); vecit++){
        for(auto it = (*vecit).begin(); it != (*vecit).end(); it++) *it = idx_val;
    }
    for(auto vecit = valSendBuff.begin(); vecit != valSendBuff.end(); vecit++){
        for(auto it = (*vecit).begin(); it != (*vecit).end(); it++) *it = 0.0;
    }
    for(auto vecit = valRecvBuff.begin(); vecit != valRecvBuff.end(); vecit++){
        for(auto it = (*vecit).begin(); it != (*vecit).end(); it++) *it = 0.0;
    }
}
#endif

#ifdef DISTRIBUTED_BASIS
    template <class T>
    void SparseMatrix<T>::construct(int rowCount, int rowPerIt){
        /*
            sendBuff layout:
            [matrix_1, matrix_2, ....]
            matrix_i is represented by a single vector_i
            matrix_i:
            [block_1,block_2,...]
            block_i is of size sendCount
            block_i:
            [row_1,row_2,...]
            row_i is of size rowCount
        */
        Timer timer;
        Timer tInit, tMPI;
        setDim(Bf->getSubDim());
        // do MPI_Alltoall communication after #rowPerIt row iterations 
        int sendCount = rowCount * rowPerIt;
        int buffSize = blockNum * sendCount;
        for(int matID=0; matID<spmNum; matID++){
            idxSendBuff.at(matID).resize(buffSize);
            idxRecvBuff.at(matID).resize(buffSize);
            valSendBuff.at(matID).resize(buffSize);
            valRecvBuff.at(matID).resize(buffSize);
        }
        
        clear();

        idx_t endRowFlag = ULLONG_MAX;

        std::vector<std::vector<int>> ridxBlockStart(blockNum);
        for(int bid = 0; bid < blockNum; bid++){
            ridxBlockStart.at(bid).resize(rowPerIt);
            for(int row=0;row<rowPerIt;row++){ridxBlockStart.at(bid).at(row)=bid*sendCount+row*rowCount;}
        }

        // std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<", block num:"<<idxSendBuff.size()<<", sendsize:"<<idxSendBuff[0].size()<<std::endl;
        for(int i = 0; i < dmNum; ++i) diagValList.at(i).resize(BaseMatrix<T>::nloc);
        for(int i = 0; i < spmNum; ++i) for(int b=0;b<blockNum;b++)rowInitList.at(i).at(b).push_back(counter.at(i).at(b));
        for (idx_t rowID = 0; rowID < BaseMatrix<T>::nloc; rowID+=rowPerIt){
            timer.tik();
            // std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<", Start row:"<<rowID<<std::endl;
            tInit.tik();
            setMpiBuff(endRowFlag); // set mpi sent/recv buff
            tInit.tok();
            idx_t iterStart = rowID;
            idx_t iterEnd = (rowID+rowPerIt)<BaseMatrix<T>::nloc ? (rowID+rowPerIt):BaseMatrix<T>::nloc;
            #pragma omp parallel for
            for(int iter = 0; iter < (iterEnd-iterStart); iter++){
                std::vector<MAP<T>> rowMaps(spmNum);
                // std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<", Start iter:"<<iterStart+iter<<std::endl;
                row(iterStart+iter,rowMaps);
                // std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<", after iter:"<<iterStart+iter<<std::endl;
                for(int matID=0; matID < spmNum; matID++){
                    int bid;
                    std::vector<int> rowBlockCount(blockNum, 0);
                    for(auto it = rowMaps.at(matID).begin(); it != rowMaps.at(matID).end(); it++){
                        if(Bi->getBid(it->first,bid)){
                            if(rowBlockCount.at(bid)<rowCount){
                                int curidx = ridxBlockStart.at(bid).at(iter)+rowBlockCount.at(bid);
                                idxSendBuff.at(matID).at(curidx) = it->first;
                                valSendBuff.at(matID).at(curidx) = it->second;
                                rowBlockCount.at(bid)++;
                            }
                            else{
                                std::cout<<"row elements count exceeds sendBuff size:"<<rowCount<<std::endl;
                                exit(1);
                            }
                        }
                    }
                }
            }
            // std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<", Start mpi all to all"<<std::endl;
            // mpi all to all 
            tMPI.tik();
            for(int matID=0; matID<spmNum; matID++){
                MPI_Alltoall(idxSendBuff.at(matID).data(), idxRecvBuff.at(matID).data(), sendCount);
                MPI_Alltoall(valSendBuff.at(matID).data(), valRecvBuff.at(matID).data(), sendCount);
            }
            tMPI.tok();
            // std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<" Start filter row:"<<rowID<<std::endl;
            // filter recv buff and push data to sparse matrix
            for(int matID=0; matID<spmNum; matID++){
                #pragma omp parallel for
                for(int bid=0; bid<blockNum; bid++){
                    for(int row=0; row<rowPerIt; row++){
                        int count_tmp = 0;
                        for(int idx=ridxBlockStart.at(bid).at(row);idx<ridxBlockStart.at(bid).at(row)+rowCount;idx++){
                            if(idxRecvBuff.at(matID).at(idx)==endRowFlag) break;
                            idx_t colID;
                            if(Bi->search(idxRecvBuff.at(matID).at(idx),colID)){
                                colList.at(matID).at(bid).push_back(colID);
                                valList.at(matID).at(bid).push_back(std::conj(valRecvBuff.at(matID).at(idx)/Bi->getNorm(colID)));
                                count_tmp++;
                            }
                        }
                        counter.at(matID).at(bid) += count_tmp;
                        rowInitList.at(matID).at(bid).push_back(counter.at(matID).at(bid));
                    }
                }
            }
            timer.tok();
            if(BaseMatrix<T>::workerID==MPI_MASTER)std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<"row:"<<rowID<<", init:"<<tInit.elapse()<<", MPI:"<<tMPI.elapse()<<", tot:"<<timer.elapse()<<"ms\n";
            // std::cout<<"workerID:"<<BaseMatrix<T>::workerID<<" finished filter row:"<<rowID<<std::endl;
        }
        for (int matID = 0; matID < spmNum; matID++) for(int bid=0; bid < blockNum; bid++){
            idx_t rowNum = bid<(blockNum-1) ? BaseMatrix<T>::nlocmax:(BaseMatrix<T>::dim - BaseMatrix<T>::nlocmax * (blockNum-1));
            MKL::create(A.at(matID).at(bid), rowNum, BaseMatrix<T>::nloc, rowInitList.at(matID).at(bid), colList.at(matID).at(bid), valList.at(matID).at(bid));
        } 
        isMatrixFree = false;
    }
#else
    template <class T>
    void SparseMatrix<T>::construct(int rowPerThread){
        clear();
        this->setDim(Bf->getSubDim());
        this->setColDim(Bi->getSubDim());
        int threadNum;
        #pragma omp parallel
        {
            #pragma omp master
            threadNum = omp_get_num_threads();
        }
        int rowPerIt = rowPerThread * threadNum;
        // for each iteration, each thread construct spmNum*rowPerThread rowMap
        std::vector<std::vector<MAP<T>>> rowMapList(rowPerIt);
        for (int i = 0; i < rowPerIt; ++i) {rowMapList[i].resize(spmNum);}

        // the starting index for each thread to copy data from rowMap to colList and valList
        std::vector<std::vector<idx_t>> startList(threadNum);
        for (int i = 0; i < threadNum; ++i) startList[i].resize(spmNum);
        // initialize containers
        for(int i = 0; i < dmNum; ++i) diagValList.at(i).resize(BaseMatrix<T>::nloc);
        for(int i = 0; i < spmNum; ++i) rowInitList.at(i).push_back(counter.at(i));
        for (idx_t rowID = BaseMatrix<T>::startRow; rowID < BaseMatrix<T>::endRow; rowID+=rowPerIt){
            #pragma omp parallel shared(rowMapList, startList, counter)
            {   
                assert(omp_get_num_threads()==threadNum);
                int threadID = omp_get_thread_num();
                int thStart = threadID * rowPerThread;
                idx_t thRowStart = rowID + thStart;
                idx_t thRowEnd = (thRowStart+rowPerThread) < BaseMatrix<T>::endRow ? (thRowStart+rowPerThread):BaseMatrix<T>::endRow;
                // clear map before construction
                for (int i = thStart; i < thStart+rowPerThread; ++i) for(int j=0;j<spmNum;j++){rowMapList[i][j].clear();}
                for(idx_t thRowID = thRowStart; thRowID<thRowEnd; thRowID++) row(thRowID, rowMapList.at(thRowID - rowID));
                #pragma omp barrier
                #pragma omp master
                {
                    for (int i = 0; i < threadNum; ++i){
                        for (int ii = 0; ii < spmNum; ++ii) startList[i][ii] = counter[ii];
                        for (int j = 0; j < rowPerThread; j++){
                            for(int jj = 0; jj < spmNum; jj++){
                                counter[jj] += rowMapList[i*rowPerThread+j][jj].size();
                                rowInitList[jj].push_back(counter[jj]);
                            }
                        }
                    }
                    for (int i = 0; i < spmNum; ++i){
                        colList[i].resize(counter[i]);
                        valList[i].resize(counter[i]);
                    }   
                }
                #pragma omp barrier
                std::vector<int> count_tmp(spmNum);
                for (int i = thStart; i < thStart + rowPerThread; ++i){
                    for (int j = 0; j < spmNum; j++){
                        for (auto it = rowMapList[i][j].begin(); it != rowMapList[i][j].end(); it++){
                            colList[j].at(startList[threadID][j]+count_tmp[j]) = it->first;
                            switch(PARTITION){
                                case ROW_PARTITION:{
                                    if constexpr (std::is_same<cdouble, T>::value) {
                                        valList[j].at(startList[threadID][j]+count_tmp[j]) = std::conj(it->second);
                                    } else if constexpr (std::is_same<double, T>::value) {
                                        valList[j].at(startList[threadID][j]+count_tmp[j]) = it->second;
                                    } else {
                                        std::cout<<"Sparse matrix only defined for double/cdouble !\n";
                                        exit(1);
                                    }
                                    break;
                                }
                                case COL_PARTITION:{
                                    valList[j].at(startList[threadID][j]+count_tmp[j]) = it->second;
                                    break;
                                }
                            }
                            count_tmp[j]++;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < spmNum; ++i) MKL::create(A.at(i), BaseMatrix<T>::nloc, Bi->getSubDim(), rowInitList.at(i), colList.at(i), valList.at(i));
        isMatrixFree = false;
    }
#endif

/*
 Version 2. M = a1*M1 + a2*M2 + ...
 */
template <class T>
void SparseMatrix<T>::MxV(T *vecIn, T *vecOut) {
    if(!isMVBuf()) setBuf();
    
    if(!isMatrixFree){
       
        #ifdef DISTRIBUTED_BASIS

            for (int id = 0; id < BaseMatrix<T>::workerNum; id++){
                // initialization
                #pragma omp parallel for
                for (idx_t i = 0; i < BaseMatrix<T>::nlocmax; ++i){
                    vecBuf[i] = 0.0;
                }
                idx_t rowStart = id * BaseMatrix<T>::nlocmax;
                idx_t rowEnd = (rowStart + BaseMatrix<T>::nlocmax)<BaseMatrix<T>::dim?(rowStart + BaseMatrix<T>::nlocmax):BaseMatrix<T>::dim;

                for (int matID = 0; matID < spmNum; matID++){
                    MKL::MxV(A.at(matID).at(id),vecIn,vecBuf.data(),parameters.at(matID));  
                }
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Reduce(vecBuf.data(), vecOut, rowEnd - rowStart, id);
            }

        #else

            #pragma omp parallel for
            for (idx_t i = 0; i < BaseMatrix<T>::nloc; ++i) vecOut[i] = 0.0;
            MPI_Allgather(vecIn,vecBuf.data(),this->getColnlocmax());
            for (int i = 0; i < spmNum; ++i) {
                MKL::MxV(A.at(i),vecBuf.data(),vecOut,parameters.at(i));
            }
        #endif

        // diagonal part
        if (dmNum>0){
            // constant diagVal
            if (!diagValList.at(0).empty()){
                #pragma omp parallel for
                for (idx_t i = 0; i < BaseMatrix<T>::nloc; ++i) vecOut[i] += diagValList[0][i]*vecIn[i];
            }
            // time dependent diagVal
            for (int j = 1; j < dmNum; j++){
                #pragma omp parallel for
                for (idx_t i = 0; i < BaseMatrix<T>::nloc; ++i) vecOut[i] += diagParameters[j]*diagValList[j][i]*vecIn[i];
            }
        }
    } else {
        MPI_Allgather(vecIn,vecBuf.data(), this->getColnlocmax());
        #pragma omp parallel for
        for (idx_t i = 0; i < BaseMatrix<T>::nloc; ++i) vecOut[i] = 0.0;
        #pragma omp parallel for
        for(idx_t rowID = BaseMatrix<T>::startRow; rowID < BaseMatrix<T>::endRow; ++rowID){
            idx_t rowID_loc = rowID - BaseMatrix<T>::startRow;
            std::vector<MAP<T>> rowMaps; rowMaps.resize(spmNum);
            row(rowID, rowMaps);
            int matID = 0;
            // off diagonal
            for(const auto& rowMap:rowMaps){
                for(const auto& kv:rowMap){
                    vecOut[rowID_loc] += parameters.at(matID) * kv.second * vecBuf.at(kv.first);
                }
                matID++;
            }
            // diagonal
            for (int j = 1; j < dmNum; j++){
                vecOut[rowID_loc] += diagParameters[j]*diagValList[j][rowID_loc]*vecIn[rowID_loc];
            }
        }
    }        
}

template<class T>
void SparseMatrix<T>::getDiag(std::vector<T> &diag) {
    diag = std::vector<T>(this->getnloc(), 0.0);
    for (int i = 0; i < dmNum; ++i) {
        #pragma omp parallel for
        for (idx_t j = 0; j < this->getnloc(); ++j) {
            diag[j] += diagValList[i][j];
        }
    }    
    for (int i = 0; i < spmNum; ++i) {
        idx_t rowid = BaseMatrix<T>::startRow;
        #pragma omp parallel for
        for (idx_t j = 0; j < this->getnloc(); ++j) {
            for (idx_t loc = rowInitList[i][j]; loc < rowInitList[i][j+1]; ++loc) {
                if (colList[i][loc] == rowid) {
                    diag[j] += valList[i][loc];
                }
            }
            ++rowid;
        }
    }
}

template <class T>
T SparseMatrix<T>::vMv(T *vecL, T *vecR){
    std::vector<T> vecTmp(BaseMatrix<T>::nlocmax);
    MxV(vecR, vecTmp.data());
    return mpiDot<T>(vecL, vecTmp.data(), BaseMatrix<T>::nloc);
}

template <class T>
void SparseMatrix<T>::print(std::string info, std::ostream& os) const {
    os<<info<<":\n";
    os<<"matrix local/tot dim: "<<this->getnloc()<<" / "<<this->getDim()<<".\n";
    os<<"non-zero elements count: "<<this->nzCount()<<".\n";
}

#endif // SparseMatrix_hpp
