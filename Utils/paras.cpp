#include "utils/paras.hpp"

Parameters::Parameters(std::string inputDir, std::vector<std::string> files) {
    for(auto const& file:files){
        if(!read(inputDir+"/"+file)) {
            std::cout<<"Parameters read err.\n";
            exit(1);
        }
    }
}

//TODO:Clear up
void Parameters::config(std::string configFile){
    clear();
    if(!read(configFile)) {
        std::cout << "Parameters read err.\n";
        exit(1);
    }
    assert(maps.find("rootDataPath")!=maps.end()); rootDataPath = maps["rootDataPath"];maps.erase("rootDataPath");
    assert(maps.find("project")!=maps.end()); project = maps["project"];maps.erase("project");
    assert(mapvecs.find("inputfiles")!=mapvecs.end());
    assert(maps.find("inputDir")!=maps.end()); std::string inputDir = maps["inputDir"];maps.erase("inputDir");
    assert(!mapvecs["inputfiles"].empty());
    for(auto const& file:mapvecs["inputfiles"]){
        if(!read(inputDir+"/"+file)){std::cout<<"Parameters read err.\n";exit(1);}
    }
    mapvecs.erase("inputfiles");
}


void Parameters::clear() {
    rootDataPath.clear(); 
    project.clear();
    
    mapi.clear(); 
    mapd.clear(); 
    maps.clear();

    mapvecd.clear(); 
    mapvecs.clear(); 
    maparrd.clear();
}

bool Parameters::read(const std::string& filename) {
    std::ifstream infile(filename);
    if(infile.is_open()){
        std::string line;
        while(std::getline(infile,line)){
            removeComment(line,'#');
            std::istringstream ins(line);
            std::string parsed;
            std::vector<std::string> vals;
            while(std::getline(ins,parsed,':')) {
                vals.push_back(parsed);
            }
            if(vals.size()==3) {
                std::string key = vals[0];
                std::string type = vals[1];
                std::istringstream ss(vals[2]);
                if (type == "b") {
                    bool val;
                    ss >> val;
                    mapb[key] = val;
                } else if (type == "i") {
                    int val;
                    ss >> val;
                    mapi[key] = val;
                } else if (type == "d") {
                    double val;
                    ss >> val;
                    mapd[key] = val;
                } else if (type == "s") {
                    std::string val;
                    ss >> val;
                    maps[key] = val;
                }
            } else if (vals.size() == 2) {
                std::string key = vals[0];
                std::string type = vals[1];
                std::string nextLine;
                if (type == "vecs") {
                    if (std::getline(infile,nextLine)) {
                        auto val = readVec<std::string>(nextLine);
                        mapvecs[key] = val;
                    }
                } else if (type == "vecd") {
                    if (std::getline(infile,nextLine)) {
                        auto val = readVec<double>(nextLine);
                        mapvecd[key] = val;
                    }
                } else if (type == "arrd") {
                    while (std::getline(infile,nextLine)) {
                        if(nextLine=="")break;
                        auto val = readVec<double>(nextLine);
                        maparrd[key].push_back(val);
                    }
                }
            }
        }
        infile.close();
        return true;
    } else {
        std::cout << filename << " failed to open!\n";
        return false;
    }
}

void Parameters::print(std::string filename) const {
    std::ofstream outfile(filename);
    print(outfile);
}
void Parameters::print(std::ostream& os) const {
    os << maps << mapi << mapd << mapvecs << mapvecd << maparrd;
}

bool opt(const Parameters &para, std::string key) {
    auto cond = para.get<bool>(key, REQUIRED::NO);
    if (cond) {
        return cond.value();
    } else {
        std::cout << key << " not found in opt(para,key) function!\n";
        return false;
    }
}