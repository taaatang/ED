#pragma once

#include <cstdio>
#include <fstream>
#include <concepts>
#include <filesystem>
#include "global/typeAlias.hpp"
#include "utils/runtimeCheck.hpp"

inline void mkdir_fs(const std::string &dir){
    std::filesystem::path p(dir);
    std::filesystem::create_directories(p);
    bool succeed = std::filesystem::is_directory(p);
    assert_msg(succeed, dir + " failed to creat!");
}

struct RedirectStdOut {
    RedirectStdOut(const std::string &fileName) {
        assert_msg(freopen(fileName.c_str(), "w", stdout) != nullptr, "failed to redirect stdout to " + fileName);
    }
    ~RedirectStdOut() { 
        fclose(stdout);
    }
};

template <typename data_t, std::integral index_t>
inline void save(const data_t *d_pt, index_t size, std::string filename, bool is_app = false, bool print = true){
    std::ofstream os;
    if (is_app) {
        os.open(filename, std::ios::binary | std::ios::app);
    } else {
        os.open(filename, std::ios::binary);
    }
    if (os.is_open()){
        os.write(reinterpret_cast<const char*>(d_pt), size * sizeof(data_t));
        os.close();
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

template <typename data_t, std::integral index_t>
inline void read(data_t *d_pt, index_t size, std::string filename){
    std::ifstream is;
    is.open(filename, std::ios::binary);
    if (is.is_open()) {
        is.read(reinterpret_cast<char*>(d_pt), size * sizeof(data_t));
        is.close();
        // std::cout << "Data loaded from " << filename << std::endl;
    } else {
        std::cout << filename << " failed to open!" << std::endl;
        exit(1);
    }
}

template <class T>
inline void read(std::vector<T> *d_pt, std::string filename) {
    std::ifstream infile;
    infile.open(filename, std::ios::in | std::ios::binary);
    if (infile.is_open()) {
        infile.seekg(0, infile.end);
        idx_t size = infile.tellg();
        infile.seekg(0, infile.beg);
        assert((size % sizeof(T))==0);
        d_pt->resize(size/sizeof(T));
        infile.read(reinterpret_cast<char*>(d_pt->data()), size);
        infile.close();
        // std::cout<<"Data loaded from "<<filename<<std::endl;
    } else {
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