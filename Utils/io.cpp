#include "io.hpp"

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

void exit_msg(std::string msg){
    std::cerr<<msg<<std::endl;
    exit(EXIT_FAILURE);
}

void assert_msg(bool condition, std::string msg){
    if(!condition){
        std::cout<<msg<<std::endl;
        exit(1);
    }
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


std::ostream& operator<<(std::ostream& os, LINK_TYPE linkt) {
    os<<"LINK TYPE:";
    switch (linkt) {
        case LINK_TYPE::SUPER_EXCHANGE_J:
            os << "SUPER_EXCHANGE_J";
            break;
        case LINK_TYPE::CHIRAL_K:
            os << "CHIRAL_K";
            break;
        case LINK_TYPE::HOPPING_T:
            os << "HOPPING_T";
            break;
        case LINK_TYPE::CHARGE_TRANSFER_V:
            os << "CHARGE_TRANSFER_V";
            break;
        case LINK_TYPE::HUBBARD_U:
            os << "HUBBARD_U";
            break;
        case LINK_TYPE::EXCHANGE_J:
            os << "EXCHANGE_J";
            break;
        case LINK_TYPE::PAIR_HOPPING_J:
            os << "PAIR_HOPPING_J";
            break;
        case LINK_TYPE::PHONON_W0:
            os << "PHONON_W0";
            break;
        case LINK_TYPE::NCHARGE_SITE_PHONON:
            os << "NCHARGE_SITE_PHONON";
            break;
        case LINK_TYPE::NCHARGE_BOND_PHONON:
            os << "NCHARGE_BOND_PHONON";
            break;
        case LINK_TYPE::HOPPING_BOND_PHONON:
            os << "HOPPING_BOND_PHONON";
            break;
        default:
            os<<"not defined!\n";
            exit(1);
            break;
    }
    return os;
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