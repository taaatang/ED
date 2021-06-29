#ifndef __LINKS_H__
#define __LINKS_H__

#include <vector>
#include <iostream>

#include "Global/globalPara.hpp"
#include "Geometry/Vec3.hpp"
#include "Geometry/Geometry.hpp"

template <class T>
class Link{
public:
    Link( ) { }
    Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, bool is_Const_ = true, bool is_Ordered_ = false);
    Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, const std::vector<Vec3d>& vecs, bool is_Const_ = true, bool is_Ordered_ = false);
    ~Link( ) { };
    
    LINK_TYPE getLinkType( ) const { return link_type; }
    T getVal( ) const { return val; }
    int getlinkid( ) const { return linkid; }
    int getmatid( ) const { return matid; }
    int getLinkNum( ) const { return LinkList.size(); }
    int getvecid(int bondid) const { return LinkVecIdList.at(bondid); }
    Vec3d getvec(int vecid) const { return LinkVecs.at(vecid); }
    // rep orbital of the link
    ORBITAL getRepOrb( ) const { return orbList.at(0); }
    std::vector<ORBITAL> getOrbs( ) const { return orbList; }
    // number of orbs in one link
    int getLinkSize( ) const { return orbList.size() - 1; }
    // size of LinkVecs
    int getLinkVecNum( ) const { return LinkVecs.size(); }

    // is this a const link or time dependent link
    bool isConst( ) const { return is_Const; }
    bool isOrdered( ) const { return is_Ordered; }
    void setConst(bool cond) { is_Const = cond; }
    void setOrdered(bool cond) { is_Ordered = cond; }

    // begin/end iterator for LinkList
    std::vector<VecI>::iterator begin( ) { return LinkList.begin(); }
    std::vector<VecI>::iterator end( ) { return LinkList.end(); }
    const std::vector<VecI>& bond( ) const { return LinkList; }

    void setVal(T val_) { val = val_; }
    void setid(int linkid_, int matid_ = 0);
    void setmatid(int matid_ = 0) { matid = matid_;}

    // add one link (contains orbids in the same link)
    void push_back(std::vector<int> orbidList) { LinkList.push_back(orbidList); }
    // add a orb's relative coord to the rep orb
    Link<T>& addLinkVec(Vec3d vec);
    //generate all the links->linkList
    void genLinkMaps(Geometry* pt_lattice);
    void print(bool brief = true, std::ostream& os = std::cout) const;

private:
    LINK_TYPE link_type{LINK_TYPE::HOPPING_T};
    T val{0.0};
    // same group ids will be contained in the same sparse matrix
    int linkid{0}, matid{0};
    bool is_Const{true}, is_Ordered{false};
    std::vector<ORBITAL> orbList;
    std::vector<Vec3d> LinkVecs;
    std::vector<VecI> LinkList;
    // LinkVecIdList[i] is the id of LinkVec for i-th bond in LinkList
    VecI LinkVecIdList;
};

Link<dataType> HeisenbergLink(std::string name, const Geometry& latt);

std::vector<Link<dataType>> HubbardSingleBandLink(const Geometry& latt);

std::vector<Link<dataType>> HubbardMultiBandLink(const Geometry& latt);

std::vector<Link<dataType>> HubbardLink(const Geometry& latt);

std::vector<Link<dataType>> RamanChannel(std::string channel, double J1, double J2, const Geometry& latt);

/********
 * LINK *
 ********/
template <class T>
Link<T>::Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, bool is_Const_ , bool is_Ordered_) {
    link_type = link_type_;
    orbList = orbList_;
    val = val_;
    is_Const = is_Const_;
    is_Ordered = is_Ordered_;
}

template <class T>
Link<T>::Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, const std::vector<Vec3d>& vecs, bool is_Const_ , bool is_Ordered_ ) {
    link_type = link_type_;
    orbList = orbList_;
    val = val_;
    is_Const = is_Const_;
    is_Ordered = is_Ordered_;
    for (auto& vec : vecs) {
        addLinkVec(vec);
    }
}

template <class T>
void Link<T>::setid(int linkid_, int matid_ ) { 
    linkid = linkid_; 
    matid = matid_; 
}

template <class T>
Link<T>& Link<T>::addLinkVec(Vec3d vec){
    LinkVecs.push_back(vec); 
    return *this;
}

template <class T>
void Link<T>::genLinkMaps(Geometry* pt_lattice) {
    Vec3d coordi, coordf; 
    for (int orbid = 0; orbid < pt_lattice->getOrbNum(); ++orbid) {
        if(pt_lattice->is_Orbital(orbid,getRepOrb())) {
            coordi = pt_lattice->getOrbR(orbid);
            for (int j = 0; j < getLinkVecNum(); j += getLinkSize()) {
                VecI tmp;
                tmp.push_back(orbid);
                bool is_bond = true;
                for (int k = j; k < j + getLinkSize(); ++k) {
                    coordf = coordi + LinkVecs.at(k);
                    int orbidf;
                    if (pt_lattice->coordToOrbid(orbList.at(k-j+1), coordf, orbidf)){
                        assert(pt_lattice->getOrb(orbidf) == orbList.at(k-j+1));
                        tmp.push_back(orbidf);
                    }else{
                        is_bond = false;
                        break;
                    }
                }
                if(is_bond){
                    push_back(tmp);
                    LinkVecIdList.push_back(j);
                }
            }
        }   
    }
}

template <class T>
void Link<T>::print(bool brief, std::ostream& os) const {
    os<<"mat id:"<<matid<<", Link id:"<<linkid<<", value:"<<val<<"\n";
    for (auto orb:orbList) os<<orb<<" ";
    os<<"\n";
    os<<"link vec:\n";
    for (const auto& vec:LinkVecs) os<<vec<<"\n";
    os<<"bond num:"<<LinkList.size()<<"\n";
    if (!brief) {
        for (const auto& bond : LinkList) os<<bond<<"\n";
    }
}

#endif // __LINKS_H__