#include "paras.hpp"

Orbital stringToOrb(std::string name, int id) {
    ORBITAL orb;
    VecD coord;
    if (name == "single") {
        orb = ORBITAL::SINGLE;
        coord = VecD{0.0, 0.0, 0.0};
    } else if (name == "dx2y2") {
        orb = ORBITAL::Dx2y2;
        coord = VecD{0.0, 0.0, 0.0};
    } else if (name == "dz2") {
        orb = ORBITAL::Dz2;
        coord = VecD{0.0, 0.0, 0.0};
    } else if (name == "px") {
        orb = ORBITAL::Px;
        coord = VecD{0.5, 0.0, 0.0};
    } else if (name == "py") {
        orb = ORBITAL::Py;
        coord = VecD{0.0, 0.5, 0.0};
    } else if (name == "pzu") {
        orb = ORBITAL::Pzu;
        coord = VecD{0.0, 0.0, 0.5};
    } else if (name == "pzd") {
        orb = ORBITAL::Pzd;
        coord = VecD{0.0, 0.0, -0.5};
    } else {
        std::cout<<"Orbital: "<<name<<", is not defined!\n"
        exit(1);
    }
    return Orbital(orb, id, coord);
}

LATTICE Parameters::getlatt() const {
    std::string name = maps["lattice type"];
    if(name == "square") {
        return LATTICE::SQUARE;
    } else if(name == "triangular") {
        return LATTICE::TRIANGULAR;
    } else {
        std::cout<<"lattice type: "<<name<<" is not defined!";
        exit(1);
    }
}

LATTICE_MODEL Parameters::getmodel() const {
    std::string name = maps["model name"];
    if (name == "Hubbard") {
        return LATTICE_MODEL::HUBBARD;
    } else if (name == "Heisenberg") {
        return LATTICE_MODEL::HEISENBERG;
    } else if (name == "tJ") {
        return LATTICE_MODEL::tJ;
    } else {
        std::cout<<"Model name: "<<name<<", not defined!";
        exit(1);
    }
}

void Parameters::config(std::string configFile){
    clear();
    if(!read(configFile)){std::cout<<"Parameters read err.\n";exit(1);}
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


void Parameters::clear(){
    rootDataPath.clear(); project.clear();
    mapi.clear(); mapd.clear(); maps.clear();
    mapvecd.clear(); mapvecs.clear(); maparrd.clear();
}

bool Parameters::read(const std::string& filename){
    std::ifstream infile(filename);
    if(infile.is_open()){
        std::string line;
        while(std::getline(infile,line)){
            removeComment(line,'#');
            std::istringstream ins(line);
            std::string parsed;
            std::vector<std::string> vals;
            while(std::getline(ins,parsed,':')){
                vals.push_back(parsed);
            }
            if(vals.size()==3){
                std::string key = vals[0];
                std::string type = vals[1];
                std::istringstream ss(vals[2]);
                if(type=="i"){int val;ss>>val;mapi[key]=val;}
                else if(type=="d"){double val;ss>>val;mapd[key]=val;}
                else if(type=="s"){std::string val;ss>>val;maps[key]=val;}
            }else if(vals.size()==2){
                std::string key = vals[0];
                std::string type = vals[1];
                std::string nextLine;
                if(type=="vecs"){
                    if(std::getline(infile,nextLine)){
                        auto val = readVec<std::string>(nextLine);
                        mapvecs[key] = val;
                    }
                }else if(type=="vecd"){
                    if(std::getline(infile,nextLine)){
                        auto val = readVec<double>(nextLine);
                        mapvecd[key] = val;
                    }
                }else if(type=="arrd"){
                    while(std::getline(infile,nextLine)){
                        if(nextLine=="")break;
                        auto val = readVec<double>(nextLine);
                        maparrd[key].push_back(val);
                    }
                }
            }
        }
        infile.close();
        return true;
    }else{
        std::cout<<filename<<" failed to open!\n";
        return false;
    }
}

void Parameters::print(std::ostream& os){
    os<<maps<<mapi<<mapd<<mapvecs<<mapvecd<<maparrd;
}