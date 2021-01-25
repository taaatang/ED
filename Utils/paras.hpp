#ifndef __PARAS_H__
#define __PARAS_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <unordered_map>
#include <utility>
#include <type_traits>

template<typename T>
using map = std::unordered_map<std::string,T>;

class Parameters{
public:
    Parameters(){};
    ~Parameters(){};

    bool read(const std::string& filename);
    void print(std::ostream& os);
private:
    map<int> mapi;
    map<double> mapd;
    map<std::string> maps;
    map<std::vector<std::string>> mapvecs;
    map<std::vector<double>> mapvecd;
    map<std::vector<std::vector<double>>> maparrd;
};

inline void removeComment(std::string& s, char delim){
    std::istringstream ins (s);
    std::string tmp;
    std::getline(ins,tmp,delim);
    s = std::move(tmp);
}

template<typename T>
std::vector<T> readVec(std::string& s, char sep=','){
    std::istringstream ins(s);
    std::vector<T> result;
    std::string s_val;
    while(std::getline(ins,s_val,sep)){
        T val;std::istringstream ins_val(s_val);ins_val>>val;
        result.push_back(val);
    }
    return result;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& myvec){
    if(!std::is_same<T,std::vector<double>>::value)os<<"\n";
    if(std::is_floating_point<T>::value)os<<std::setprecision(2)<<std::fixed;
    for(int i=0;i<myvec.size()-1;i++){
        os<<myvec.at(i);
        if(!std::is_same<T,std::vector<double>>::value)os<<", ";
    }
    if(!myvec.empty())os<<myvec.back();
    if(std::is_floating_point<T>::value)os.unsetf(std::ios::floatfield); 
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const map<T>& mymap){
    for(const auto& pair:mymap){
        os<<pair.first<<":"<<pair.second<<"\n";
    }
    return os;
}
#endif // __PARAS_H__