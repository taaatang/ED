//
//  Geometry.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Geometry_hpp
#define Geometry_hpp

#include "../Global/globalPara.hpp"
#include "../Utils/utils.hpp"
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>


/*
   ***************************
   *   Geometry Base Class   *
   ***************************
*/
struct Orbital{
    /*
        Describe an orbital inside a unit site.
    */
    int id; // 0 ~ Norb - 1
    int orbid; // 0 ~ unit cell size - 1
    int siteid; // 0 ~ Nsite -1
    ORBITAL orb;
    VecD coord;
    Orbital(ORBITAL orb_=ORBITAL::SINGLE, int orbid_=0, VecD coord_=VecD(3,0.0)):orbid(orbid_),orb(orb_){
        coord = coord_;
    }
};
struct Site{
    /*
     Describe a single unit site of a Lattice
     id holds the index of this site among all sites
     coord is the coordinate in the Lattice basis
     */
    int id; // 0 ~ Nsite-1
    VecD coord; // in unit of ax, ay; in x-y coordinates
    Site(int id_, std::vector<double> coord_):id(id_){
        coord = coord_;
    }
};

class Geometry{
protected:
    /*
        Instantiate an object describing the real/k space geometry of a Lattice
    */
    std::string name;
    int dim;
    int Nsite; // Number of sites
    int Norb; // Number of orbitals
    int Norb_enlg; // Number of orbitals in the enlarged Lattice (Used for mapping Site interacting with boundary Site back to Lattice).
    bool is_PBC; // Boundary condition
    bool is_TBC{false}; // twisted boundary condition
    VecD phase{0.0,0.0,0.0};
    VecD dPhase{0.0,0.0,0.0};
    std::vector<std::vector<int>> TransList; // Nsite*Norb. TransList[r][orbid] is the orbid' obtained by translating orbid by a vector r.
    std::vector<std::vector<cdouble>> TransPhaseList;
    PointGroup PG;
    std::vector<std::vector<int>> PGList; 
    std::vector<std::vector<cdouble>> CharacterList;
    
    // a1, a2 Lattice basis. ax, ay the corresponding values of Lattice basis in xy coordinates
    VecD a1, a2, a3, ax, ay, az;
    // unit vectors for the super Lattice
    VecD R1, R2, R3; 
    std::vector<VecD> TranVecs;
    // dual Lattice basis, ai*bj=2pi*delta(i,j)
    VecD b1, b2, b3, bx, by, bz; 
    VecD b10, b20, b30;
    double tol{1e-8};

    VecD xlist, ylist, zlist, kxlist, kylist, kzlist;

    // rotation and reflection center
    VecD center;

    std::vector<Site> Lattice; // A vector of Site in the Lattice
    std::vector<Site> KLattice; // Dual K space Lattice
    std::vector<Orbital> unitSite, boundary, orbs, enlgOrbs; // vectos of orbitals and orbitals in enlarged Lattice

public:
    Geometry(int dim_=3):dim(dim_), Nsite(0), Norb(0), Norb_enlg(0), is_PBC(true), PG(PointGroup::NONE){
        center = VecD{0.0,0.0,0.0};
        a1 = VecD{1.0,0.0,0.0}; a2 = VecD{0.0,1.0,0.0}; a3 = VecD{0.0,0.0,1.0};
        ax.resize(dim,0.0); ay.resize(dim,0.0); az.resize(dim,0.0);
        R1.resize(dim,0.0); R2.resize(dim,0.0); R3.resize(dim,0.0);
        b1 = VecD{1.0,0.0,0.0}; b2 = VecD{0.0,1.0,0.0}; b3 = VecD{0.0,0.0,1.0};
        b10 = b1; b20 = b2; b30 = b3;
        bx.resize(dim,0.0); by.resize(dim,0.0); bz.resize(dim,0.0);
    }
    ~Geometry(){};

    // add an Orbital to the unit cell
    Geometry& addOrb(Orbital orb) {unitSite.push_back(orb);return *this;}
    Geometry& addBoundary(Orbital orb) {assert(!is_PBC);boundary.push_back(orb);return *this;}
    // lattice dimension. default is 3
    int getDim() const {return dim;}
    // id of an orbital (0~unitCell size - 1)
    VecI getOrbID(ORBITAL orb) const {VecI ids; for(auto it = unitSite.begin(); it != unitSite.end(); it++) if (orb==(*it).orb) ids.push_back((*it).orbid); return ids;}
    // test if the orbital at position id is orb_test
    bool is_Orbital(int id, ORBITAL orb_test) const {return orbs.at(id).orb==orb_test;}
    // return the orbital at position id
    ORBITAL getOrb(int id) const {return orbs.at(id).orb;}
    // return the siteid unitcell coord
    void getSiteR(int siteid, double* r) const {for (int i = 0; i < dim; i++) r[i] = Lattice.at(siteid).coord[i];}
    // return the id orbital coord
    void getOrbR(int id, double* r) const {for (int i = 0; i < dim; i++) r[i] = orbs.at(id).coord[i];}
    // return the kid site coord
    void getK(int kid, double* k) const {for (int i = 0; i < dim; i++) k[i] = KLattice.at(kid).coord[i];}
    void getSiteRxy(int siteid, double* rxy) const {assert(dim==3); vecXAdd(Lattice.at(siteid).coord[0],ax.data(),Lattice.at(siteid).coord[1],ay.data(),Lattice.at(siteid).coord[2],az.data(),rxy,dim);}
    void getOrbRxy(int orbid, double* rxy) const {assert(dim==3); vecXAdd(orbs.at(orbid).coord[0],ax.data(),orbs.at(orbid).coord[1],ay.data(),orbs.at(orbid).coord[2],az.data(),rxy,dim);}
    VecD RtoRxy(VecD R) const {VecD Rxy(3); vecXAdd(R.at(0),ax.data(),R.at(1),ay.data(),R.at(2),az.data(),Rxy.data(),dim); return Rxy;}
    double RdotR(VecD R1, VecD R2)const {return vdotv(RtoRxy(R1), RtoRxy(R2));}
    void getKxy(int kid, double* kxy) const {assert(dim==3); vecXAdd(KLattice.at(kid).coord[0],bx.data(),KLattice.at(kid).coord[1],by.data(),KLattice.at(kid).coord[2],bz.data(),kxy,dim);}
    // return the id-th orbital's unit cell id
    int getSiteid(int id) const {return orbs.at(id).siteid;}
    // exp(i*2pi*k*r)
    cdouble expKR(int kid, int siteid) const {if (kid==-1) return 1.0; return std::exp(2*PI*CPLX_I*(KLattice.at(kid).coord[0]*Lattice.at(siteid).coord[0] + KLattice.at(kid).coord[1]*Lattice.at(siteid).coord[1]));}
    cdouble getChi(int PGRepInd, int OpInd) const {return CharacterList.at(PGRepInd).at(OpInd);}
    // Lattice name
    std::string getName() const {return name;}
    // number of unit cell
    int getSiteNum() const {return Nsite;}
    // number of orbitals in a unit cell
    int getUnitOrbNum() const {return unitSite.size();}
    // total number of orbitals
    int getOrbNum() const {return Norb;}
    // return the id' obtained by translate id by r
    int getOrbTran(int r, int id) const {return TransList.at(r).at(id);}
    cdouble getOrbTranPhase(int r, int id) const {return TransPhaseList.at(r).at(id);}
    // return number of ireps of the point group
    int getPGRepNum() const {return CharacterList.size();}
    int getPGOpNum(int PGRepInd) const {return CharacterList.at(PGRepInd).size();}
    // return the id' obtained by transform id by r-th Point Group Element
    int getOrbPG(int r, int id) const {return PGList.at(r).at(id);}
    // a1-a2 coord -> orbid
    bool coordToOrbid(double* coord, int &orbid) const;
    // orbital occupancy count
    void orbOCC(VecI& vec, VecI& occ) const {occ = VecI(getUnitOrbNum(),0); for(int i = 0; i < getOrbNum(); i++) if(vec[i]) occ.at(orbs.at(i).id) += 1;};
    void orbOCC(VecI& vecu, VecI& vecd, VecI& occ) const {occ = VecI(getUnitOrbNum(),0);for(int i = 0; i < getOrbNum(); i++){if(vecu[i])occ.at(orbs.at(i).orbid) += 1;if(vecd[i])occ.at(orbs.at(i).orbid) += 1;}}
    void orbOCC(VecI& vecu, VecI& vecd, VecI& occ, VecI& docc) const {
        occ = VecI(getUnitOrbNum(),0); docc = VecI(getUnitOrbNum(),0);
        for(int i = 0; i < getOrbNum(); i++) {
            if(vecu[i]) {
                occ.at(orbs.at(i).orbid) += 1;
                if(vecd[i])docc.at(orbs.at(i).orbid) += 1;
            }
            if(vecd[i])occ.at(orbs.at(i).orbid) += 1;
        }
    };
    // binary rep
    void orbOCC(ind_int repI, VecI& occ) const {occ = VecI(getUnitOrbNum(),0); for(int i = 0; i < getOrbNum(); i++) if(bitTest(repI,i)) occ.at(orbs.at(i).id) += 1;};
    void orbOCC(pairIndex pairRepI, VecI& occ) const {occ = VecI(getUnitOrbNum(),0);for(int i = 0; i < getOrbNum(); i++){if(bitTest(pairRepI.first,i))occ.at(orbs.at(i).orbid) += 1;if(bitTest(pairRepI.second,i))occ.at(orbs.at(i).orbid) += 1;}}
    void orbOCC(pairIndex pairRepI, VecI& occ, VecI& docc) const {
        occ = VecI(getUnitOrbNum(),0); docc = VecI(getUnitOrbNum(),0);
        for(int i = 0; i < getOrbNum(); i++) {
            if(bitTest(pairRepI.first,i)) {
                occ.at(orbs.at(i).orbid) += 1;
                if(bitTest(pairRepI.second,i))docc.at(orbs.at(i).orbid) += 1;
            }
            if(bitTest(pairRepI.second,i))occ.at(orbs.at(i).orbid) += 1;
        }
    };

    void printLattice() const;
    void printKLattice() const;
    void printOrbs() const;
    void printTrans() const;
    void printPG() const;
    void print() const;
    
    void setPBC(bool PBC){is_PBC = PBC;}
    void setName(std::string s){name = s;}
    // generate translation map for all r. TransList.at(r).at(id) = id->r
    bool rotate(int orbid, int& orbidf) const;
    VecD rotate(VecD vec) const;
    bool reflect(int orbid, int& orbidf) const;
    bool mirror(int orbid, int& orbidf) const;
    void genPGList();
    void genTransList();
    bool crossBoundx(const VecD& coord) const;
    bool crossBoundy(const VecD& coord) const;
    // construct Lattice, KLattice, orbs, enlg_orbs and TransList
    void construct();
};

class TriAngLattice: public Geometry{
public:
    // constructor for Lattice with D6 point group
    TriAngLattice(int N, bool PBC =  true);
    // constructor for N1*N2 Lattice
    TriAngLattice(int N1, int N2, bool PBC = true);
    ~TriAngLattice(){};
};

class SquareLattice: public Geometry{
public:
    SquareLattice(int N, bool PBC=true);
    SquareLattice(int N1, int N2, bool PBC=true, bool TBC=false, double phase_x=0.0, double phase_y=0.0);
    ~SquareLattice(){};
};

#endif
