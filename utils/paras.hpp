#pragma once

#include <iostream>
#include <iomanip>
#include <optional>
#include <map>
#include <type_traits>

#include "global/typeAlias.hpp"
#include "utils/io.hpp"
#include "utils/runtimeCheck.hpp"


template<typename T>
using map = std::map<std::string, T>;

template<typename T>
bool check(const map<T> &dic, std::string key) {
    return (dic.find(key) != dic.end());
}

enum class REQUIRED {YES, NO};

/**
 * @brief Parameter class. Handling the input
 * 
 */
class Parameters{
public:
    Parameters() = default;

    Parameters(std::string configFile){config(configFile);}

    Parameters(std::string InputDir, std::vector<std::string> files);

    void config(std::string configFile);

    void clear();

    bool read(const std::string& filename);

    void print(std::string filename) const;

    void print(std::ostream& os) const;

    template <typename T>
    std::optional<T> get(std::string key, REQUIRED opt = REQUIRED::YES) const;

private:
    std::string rootDataPath;

    std::string project;

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
std::ostream& operator<<(std::ostream& os, const map<T>& mymap){
    for(const auto& pair : mymap){
        os << pair.first << ":" << pair.second << std::endl;
    }
    return os;
}