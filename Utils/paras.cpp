#include "paras.hpp"

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
    os<<mapi<<mapd<<maps<<mapvecs<<mapvecd<<maparrd;
}