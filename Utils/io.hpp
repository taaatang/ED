#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <omp.h>
// for filesystem. if not c++ 17, use boost library
#define CPP_17_FS

#ifdef CPP_17_FS
    #include <filesystem>
#else
    #include <boost/filesystem.hpp>
#endif //CPP_17_FS

#include "global/globalType.hpp"
#include "global/constant.hpp"
#include "utils/mpiwrap.hpp"

#define LOCATION(cond) if (cond) {std::cout << "line number: "<< __LINE__ << "\nfrom func: " << __func__ << "\nfrom file: "<< __FILE__ << '\n';}

/**
 * @brief print option
 * SILENT: don't print
 * BRIEF: print short info
 * VERBOSE: print detailed info
 * 
 */
enum class PRINT {SILENT, BRIEF, VERBOSE};
/*
    ****************
    * Display Info *
    ****************
*/

void OMP_Info(int workerID);

void mpi_info(int& workerID, int& workerNum);

inline void init(int &workerID, int &workerNum, bool &isMaster) {
    MPI_Init(NULL, NULL);
    mpi_info(workerID, workerNum);
    isMaster = (workerID == 0);
}

void exit_msg(std::string msg);

void assert_msg(bool condition, std::string msg);

void printModel(LATTICE_MODEL model);

std::string tostr(double val, int digit = 2);

std::string tostr(int val);

void tolower(std::string &str);

void toupper(std::string &str);

void printLine(int n = 50, char c = '*');

std::ostream& operator<<(std::ostream& os, LATTICE_MODEL model);
std::ostream& operator<<(std::ostream& os, LINK_TYPE linkt);
std::ostream& operator<<(std::ostream& os, SPIN s);
std::ostream& operator<<(std::ostream& os, LADDER t);

std::ostream& operator<<(std::ostream& os, const VecD& vec);
std::ostream& operator<<(std::ostream& os, const VecI& vec);
std::ostream& operator<<(std::ostream& os, const std::vector<cdouble>& vec);

#ifdef CPP_17_FS
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
#endif //CPP_17_FS

template <class T>
inline void save(T *d_pt, int size, std::ofstream *f_pt, std::string filename, bool is_app = false, bool print = true){
    if(is_app)f_pt->open(filename, std::ios::binary|std::ios::app);
    else f_pt->open(filename, std::ios::binary);
    if (f_pt->is_open()){
        f_pt->write(reinterpret_cast<char*>(d_pt), size * sizeof(T));
        f_pt->close();
        if (print) {
            if (is_app) {
                std::cout << "Data appended to " << filename << std::endl;
            } else {
                std::cout << "Data wrote to " << filename << std::endl;
            }
        }    
    } else {
        std::cout << filename << " failed to open!" << std::endl;
        exit(1);
    }
}

template <class T>
inline void save(T *d_pt, idx_t size, std::ofstream *f_pt, std::string filename, bool is_app = false, bool print = true){
    if(is_app)f_pt->open(filename, std::ios::binary|std::ios::app);
    else f_pt->open(filename, std::ios::binary);
    if (f_pt->is_open()){
        f_pt->write(reinterpret_cast<char*>(d_pt), size * sizeof(T));
        f_pt->close();
        if (print) {
            if (is_app) {
                std::cout << "Data appended to " << filename << std::endl;
            } else {
                std::cout << "Data wrote to " << filename << std::endl;
            }
        }    
    }else{
        std::cout << filename << " failed to open!" << std::endl;
        exit(1);
    }
}

template <class T>
void save(T *d_pt, int size, std::string filename, bool is_app = false, bool print = true) {
    std::ofstream os;
    save<T>(d_pt, size, &os, filename, is_app, print);
}

template <class T>
void save(T *d_pt, idx_t size, std::string filename, bool is_app = false, bool print = true) {
    std::ofstream os;
    save<T>(d_pt, size, filename, is_app, print);
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
