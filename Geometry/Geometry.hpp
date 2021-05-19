//
//  Geometry.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Geometry_hpp
#define Geometry_hpp

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <complex>
#include <cassert>

#include "Global/globalPara.hpp"
#include "Geometry/Vec3.hpp"
#include "Algebra/algebra.hpp"
#include "Utils/bitop.hpp"


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
    Vec3d coord;
    Orbital(ORBITAL orb_=ORBITAL::SINGLE, int orbid_=0, Vec3d coord_=Vec3d{0.0, 0.0, 0.0}):orbid(orbid_),orb(orb_){
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
    Vec3d coord; // in unit of ax, ay; in x-y coordinates
    Site(int id_, Vec3d coord_):id(id_){
        coord = coord_;
    }
};

enum LATTICE {CHAIN, SQUARE, TRIANGULAR};

class Geometry {
public:
    Geometry( );
    virtual ~Geometry( ) { }

    bool check( ) const;

    // add an Orbital to the unit cell
    Geometry& addOrb(Orbital orb);
    Geometry& addBoundary(Orbital orb);

    LATTICE getType() const { return type; }
    // lattice dimension. default is 3
    int getDim( ) const { return dim; }

    const std::vector<Orbital>& getUnitCell( ) const { return unitSite; }
    
    // id of an orbital (0~unitCell size - 1)
    VecI getOrbID(ORBITAL orb) const;
    
    // test if the orbital at position id is orb_test
    bool is_Orbital(int id, ORBITAL orb_test) const { return orbs.at(id).orb == orb_test; }
    bool is_Orbital(int id, ORBITAL orb_test, int orbid) const;
    
    // return the orbital at position id
    ORBITAL getOrb(int id) const { return orbs.at(id).orb; }

    PointGroup getPG( ) const { return PG; }
    // return the siteid unitcell coord
    Vec3d getSiteR(int siteid) const { return Lattice.at(siteid).coord; }
    // return the id orbital coord
    Vec3d getOrbR(int id) const { return orbs.at(id).coord; }
    // return the kid site coord
    Vec3d getK(int kid) const { return KLattice.at(kid).coord; }

    Vec3d RtoRxy(const Vec3d &R) const;
    Vec3d KtoKxy(const Vec3d &K) const;
    Vec3d getSiteRxy(int siteid) const;
    Vec3d getOrbRxy(int orbid) const;
    Vec3d getKxy(int kid) const;
    
    double RdotR(Vec3d R1, Vec3d R2) const { return dot(RtoRxy(R1), RtoRxy(R2)); }
    
    // return the id-th orbital's unit cell id
    int getSiteid(int id) const {return orbs.at(id).siteid;}
    // exp(i*2pi*k*r)
    cdouble expKR(int kid, int siteid) const;
    cdouble twistPhase(int orbI, int orbJ) const; 
    cdouble getChi(int PGRepIdx, int OpIdx) const { return CharacterList.at(PGRepIdx).at(OpIdx); }
    
    // number of unit cell
    int getSiteNum( ) const { return Nsite; }
    // number of orbitals in a unit cell
    int getUnitOrbNum( ) const { return unitSite.size(); }
    // total number of orbitals
    int getOrbNum( ) const { return Norb; }
    // return the id' obtained by translate id by r
    int getOrbTran(int r, int orbIdx) const { return TransList.at(r).at(orbIdx); }
    cdouble getOrbTranPhase(int r, int id) const { return TransPhaseList.at(r).at(id); }
    // return number of ireps of the point group
    int getPGRepNum() const { return CharacterList.size(); }
    int getPGOpNum(int PGRepIdx) const { return CharacterList.at(PGRepIdx).size(); }
    // return the id' obtained by transform id by r-th Point Group Element
    int getOrbPG(int r, int id) const { return PGList.at(r).at(id); }
    // a1-a2 coord -> orbid
    bool coordToOrbid(ORBITAL orb, const Vec3d &coord, int &orbid) const;
    // orbital occupancy count
    void orbOCC(VecI &vec, VecI &occ) const; 
    void orbOCC(VecI &vecu, VecI &vecd, VecI &occ) const; 
    void orbOCC(VecI &vecu, VecI &vecd, VecI &occ, VecI &docc) const;
    // binary rep
    void orbOCC(idx_t repI, VecI &occ) const;
    void orbOCC(pairIdx_t pairRepI, VecI &occ) const;
    void orbOCC(pairIdx_t pairRepI, VecI &occ, VecI &docc) const;

    void printLattice() const;
    void printKLattice() const;
    void printOrbs() const;
    void printTrans() const;
    void printPG() const;
    void print() const;
    
    // generate translation map for all r. TransList.at(r).at(id) = id->r
    bool rotate(int orbid, int &orbidf) const;
    Vec3d rotate(const Vec3d &v) const;
    bool reflect(int orbid, int &orbidf) const;
    bool mirror(int orbid, int &orbidf) const;
    void genPGList();
    void genTransList();
    bool crossBoundx(const Vec3d &coord) const;
    bool crossBoundy(const Vec3d &coord) const;
    // construct Lattice, KLattice, orbs, enlg_orbs and TransList
    void construct();

protected:
    /*
        Instantiate an object describing the real/k space geometry of a Lattice
    */
    LATTICE type{LATTICE::SQUARE};
    std::string name;
    int dim{3};
    int Nsite; // Number of sites
    int Norb; // Number of orbitals
    // Number of orbitals in the enlarged Lattice (Used for mapping Site interacting with boundary Site back to Lattice).
    int Norb_enlg;
    bool is_PBC; // Boundary condition
    bool is_TBC{false}; // twisted boundary condition
    Vec3d phase{0.0, 0.0, 0.0};
    Vec3d dPhase{0.0, 0.0, 0.0};
    // Nsite*Norb. TransList[r][orbid] is the orbid' obtained by translating orbid by a vector r.
    std::vector<VecI> TransList;
    std::vector<VecZ> TransPhaseList;
    PointGroup PG;
    std::vector<VecI> PGList; 
    std::vector<VecZ> CharacterList;
    
    // a1, a2 Lattice basis. ax, ay the corresponding values of Lattice basis in xy coordinates
    Vec3d a1, a2, a3, ax, ay, az;
    // unit vectors for the super Lattice
    Vec3d R1, R2, R3; 
    std::vector<Vec3d> TranVecs;
    // dual Lattice basis, ai*bj=2pi*delta(i,j)
    Vec3d b1, b2, b3, bx, by, bz; 
    Vec3d b10, b20, b30;
    double tol{1e-8};

    VecD xlist, ylist, zlist, kxlist, kylist, kzlist;

    // rotation and reflection center
    Vec3d center;
    // A vector of Site in the Lattice
    std::vector<Site> Lattice;
    // Dual K space Lattice 
    std::vector<Site> KLattice;
    // vectos of orbitals and orbitals in enlarged Lattice
    std::vector<Orbital> unitSite, boundary, orbs, enlgOrbs;
};

class TriAngLattice: public Geometry{
public:
    // constructor for Lattice with D6 point group
    TriAngLattice(int N, bool PBC =  true);
    // constructor for N1*N2 Lattice
    TriAngLattice(int N1, int N2, bool PBC = true);
    ~TriAngLattice( ) { };
};

class SquareLattice: public Geometry{
public:
    SquareLattice(int N, bool PBC=true);
    SquareLattice(int N1, int N2, bool PBC=true, bool TBC=false, double phase_x=0.0, double phase_y=0.0);
    ~SquareLattice( ) { };
};

#endif
