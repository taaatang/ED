#ifndef __PARAS_H__
#define __PARAS_H__

#include <iostream>
#include <memory>
#include <optional>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include <unordered_map>
#include <utility>
#include <type_traits>
#include <assert.h>

#include "Global/globalPara.hpp"
#include "Geometry/Geometry.hpp"
#include "Basis/Basis.hpp"
#include "Operator/OperatorsBase.hpp"
#include "Operator/Operators.hpp"
#include "Operator/links.hpp"

template<typename T>
using map = std::map<std::string,T>;

enum class REQUIRED {YES, NO};

template<typename T>
bool check(const map<T> &dic, std::string key) {
    return (dic.find(key) != dic.end());
}

/**
 * @brief Parameter class. Handling the input
 * 
 */
class Parameters{
    friend class Path;
    template <typename>
    friend struct System;
    friend void setlatt(const Parameters&, std::unique_ptr<Geometry>&);
    friend void setbasis(const Parameters&, std::unique_ptr<Basis>&, Geometry*);
    friend void setbasis(const Parameters&, std::unique_ptr<Basis>&, Geometry*, int, int, int, int);
    friend void setham(const Parameters&, std::unique_ptr<HamiltonianBase<dataType>>&, Geometry*, Basis*);
    friend void setmeasure(const Parameters&);
    friend void setpulse(const Parameters&, Pulse&);
    friend bool opt(const Parameters&, std::string key);
public:
    Parameters(){}
    Parameters(std::string configFile){config(configFile);}
    Parameters(std::string InputDir, std::vector<std::string> files);
    ~Parameters(){};
    void config(std::string configFile);
    void clear();
    bool read(const std::string& filename);
    void print(std::string filename) const;
    void print(std::ostream& os) const;

    template <typename T>
    std::optional<T> get(std::string key, REQUIRED opt = REQUIRED::YES) const;

    LATTICE getlatt() const;
    LATTICE_MODEL getmodel() const;

private:
    std::string rootDataPath,project;
    map<bool> mapb;
    map<int> mapi;
    map<double> mapd;
    map<Str> maps;
    map<VecStr> mapvecs;
    map<VecD> mapvecd;
    map<ArrD> maparrd;
};

template <typename T>
std::optional<T> Parameters::get(std::string key, REQUIRED opt) const {
    if constexpr (       std::is_same<T, bool>::value) {
        if (check(mapb, key)) return mapb.at(key);
    } else if constexpr (std::is_same<T, int>::value) {
        if (check(mapi, key)) return mapi.at(key);
    } else if constexpr (std::is_same<T, double>::value) {
        if (check(mapd, key)) return mapd.at(key);
    } else if constexpr (std::is_same<T, Str>::value) {
        if (check(maps, key)) return maps.at(key);
    } else if constexpr (std::is_same<T, VecStr>::value) {
        if (check(mapvecs, key)) return mapvecs.at(key);
    } else if constexpr (std::is_same<T, VecD>::value) {
        if (check(mapvecd, key)) return mapvecd.at(key);
    } else if constexpr (std::is_same<T, ArrD>::value) {
        if (check(maparrd, key)) return maparrd.at(key);
    }
    if (opt == REQUIRED::YES) {
        LOCATION(true);
        assert_msg(false, key + " not found!");
    }
    return std::nullopt;
}

void setpath(Parameters&);
void setlatt(const Parameters&, std::unique_ptr<Geometry>& latt);
void setbasis(const Parameters&, std::unique_ptr<Basis>&, Geometry*);
void setbasis(const Parameters&, std::unique_ptr<Basis>&, Geometry*, int nuf, int ndf, int kf, int pf);
void setham(const Parameters&, std::unique_ptr<HamiltonianBase<dataType>>& H, Geometry*, Basis*);
void setBasics(const Parameters&, std::unique_ptr<Geometry>& latt, std::unique_ptr<Basis>& B, std::unique_ptr<HamiltonianBase<dataType>>& H);

void setmeasure(const Parameters&);
void setpulse(const Parameters&, Pulse&);

void setPeierls(OperatorBase<dataType>& H, const std::vector<double>& pol);

bool opt(const Parameters &para, std::string key);

inline void removeComment(std::string& s, char delim = '#'){
    std::istringstream ins (s);
    std::string tmp;
    std::getline(ins,tmp,delim);
    s = std::move(tmp);
}

template<typename T>
std::vector<T> readVec(std::string& s, char sep=','){
    removeComment(s);
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
    for(int i = 0; i< (int)myvec.size()-1; ++i){
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

Orbital stringToOrb(std::string name, int id);

#endif // __PARAS_H__