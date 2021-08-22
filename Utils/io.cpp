#include "utils/io.hpp"
#include "global/constant.hpp"

/********
 * INFO *
 ********/

void OMP_Info(int workerID){
    int ompThreadsNum;
    #pragma omp parallel
    {
        #pragma omp master
        ompThreadsNum = omp_get_num_threads();
    }
    if (workerID==MPI_MASTER) std::cout<<"openMP turned on with "<<ompThreadsNum<<" threads"<<std::endl;
}

void mpi_info(int& workerID, int& workerNum){
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID==MPI_MASTER) std::cout<<"Total MPI Workers:"<<workerNum<<std::endl;
    OMP_Info(workerID);
}


/******
 * IO *
 ******/

std::string tostr(double val, int digit){
    std::ostringstream strTmp;
    strTmp<<std::fixed<<std::setprecision(digit)<<val;
    return strTmp.str();
}

std::string tostr(int val){
    return std::to_string(val);
}

void tolower(std::string &str) {
    for (auto &c : str) {
        c = tolower(c);
    }
}

void toupper(std::string &str) {
    for (auto &c : str) {
        c = toupper(c);
    }
}

void printLine(int n, char c) {
    std::cout << std::string(n, c) << '\n';
}

std::ostream& operator<<(std::ostream& os, const VecD& vec) {
    os<<"[";
    for (auto val:vec) {
        os<<" "<<val;
    }
    os<<" ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<cdouble>& vec) {
    os<<"[";
    for (auto val:vec) {
        os << " " << (std::abs(std::real(val)) > 1e-12 ? std::real(val) : 0) << "+" << (std::abs(std::imag(val)) > 1e-12 ? std::imag(val) : 0) << "i";
    }
    os<<" ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const VecI& vec) {
    os<<"[";
    for (auto val:vec) {
        os<<" "<<val;
    }
    os<<" ]";
    return os;
}