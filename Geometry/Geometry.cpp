//
//  Geometry.cpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "Geometry.hpp"


/*
   ***************************
   *   Geometry Base Class   *
   ***************************
*/

Geometry::Geometry( ):Nsite(0), Norb(0), Norb_enlg(0), is_PBC(true), PG(PointGroup::NONE) {
    center = Vec3d{0.0,0.0,0.0};
    a1 = Vec3d{1.0,0.0,0.0}; a2 = Vec3d{0.0,1.0,0.0}; a3 = Vec3d{0.0,0.0,1.0};
    ax.fill(0.0); ay.fill(0.0); az.fill(0.0);
    R1.fill(0.0); R2.fill(0.0); R3.fill(0.0);
    b1 = Vec3d{1.0,0.0,0.0}; b2 = Vec3d{0.0,1.0,0.0}; b3 = Vec3d{0.0,0.0,1.0};
    b10 = b1; b20 = b2; b30 = b3;
    bx.fill(0.0); by.fill(0.0); bz.fill(0.0);
}

Geometry& Geometry::addOrb(Orbital orb) { 
    unitSite.push_back(orb); 
    return *this; 
}

Geometry& Geometry::addBoundary(Orbital orb) { 
    assert(!is_PBC); 
    boundary.push_back(orb); 
    return *this;
}

VecI Geometry::getOrbID(ORBITAL orb) const {
    VecI ids; 
    for (const auto& uorb:unitSite) {
        if (orb==uorb.orb) ids.push_back(uorb.orbid);
    }
    return ids;
}

bool Geometry::is_Orbital(int id, ORBITAL orb_test, int orbid) const {
    return orbs.at(id).orb == orb_test and orbs.at(id).orbid == orbid; 
}


Vec3d Geometry::RtoRxy(const Vec3d &R) const {
    return R[0] * ax + R[1] * ay + R[2] * az;
}

Vec3d Geometry::KtoKxy(const Vec3d &K) const {
    return K[0] * bx + K[1] * by + K[2] * bz;
}

Vec3d Geometry::getSiteRxy(int siteid) const {
    const Vec3d &v = Lattice.at(siteid).coord;
    return RtoRxy(v);
}

Vec3d Geometry::getOrbRxy(int orbid) const {
    const Vec3d &v = orbs.at(orbid).coord;
    return RtoRxy(v);
}

Vec3d Geometry::getKxy(int kid) const {
    const Vec3d &k = KLattice.at(kid).coord;
    return KtoKxy(k);
 }

cdouble Geometry::expKR(int kid, int siteid) const {
    if (kid == -1) {
        return 1.0; 
    }
    return std::exp(2 * PI * CPLX_I * dot(KLattice.at(kid).coord, Lattice.at(siteid).coord));
}

cdouble Geometry::twistPhase(int orbI, int orbJ) const {
    if (!is_TBC) {
        return 1.0;
    }
    int siteI = getSiteid(orbI);
    int siteJ = getSiteid(orbJ);
    const Vec3d &ri = Lattice.at(siteI).coord;
    const Vec3d &rj = Lattice.at(siteJ).coord;
    return std::exp(CPLX_I * dot(dPhase, rj - ri));
}

Generator<cdouble> Geometry::getGT(int kidx) const {
    Generator<cdouble> g;
    if (kidx == -1) {
        Transform<cdouble> t(1.0, getOrbNum());
        for (int i = 0; i < (int)t.size; ++i) {
            t[i] = i;
        }
        g.add(t);
    } else {
        for (int r = 0; r < getSiteNum(); ++r) {
            Transform<cdouble> t(expKR(kidx, r)/cdouble(getSiteNum()), getOrbNum());
            for (int i = 0; i < getOrbNum(); ++i) {
                t[i] = getOrbTran(r, i);
            }
            g.add(t);
        }
    }
    return g;
}

Generator<cdouble> Geometry::getGP(int pidx) const {
    Generator<cdouble> g;
    if (pidx == -1) {
        Transform<cdouble> t(1.0, getOrbNum());
        for (int i = 0; i < (int)t.size; ++i) {
            t[i] = i;
        }
        g.add(t);
    } else {
        for (int p = 0; p < getPGOpNum(pidx); ++p) {
            Transform<cdouble> t(getChi(pidx, p)/cdouble(getPGOpNum(pidx)), getOrbNum());
            for (int i = 0; i < getOrbNum(); ++i) {
                t[i] = getOrbPG(p, i);
            }
            g.add(t);
        }
    }
    return g;
}

bool Geometry::coordToOrbid(ORBITAL orb, const Vec3d &coord, int &orbid) const {
    bool cond;
    for (const auto& Orb:enlgOrbs) {
        if (Orb.orb != orb) continue;
        cond = true;
        for (int i = 0; i < 3; ++i) {
            if(std::abs(coord[i]-Orb.coord[i])>tol) {
                cond = false; 
                break;
            }
        }
        if (cond) {
            orbid = Orb.id; 
            return cond;
        } else {
            continue;
        }
    }
    return false;
}

void Geometry::orbOCC(VecI& vec, VecI& occ) const {
    occ = VecI(getUnitOrbNum(),0); 
    for(int i = 0; i < getOrbNum(); ++i) {
        if (vec[i]) occ.at(orbs.at(i).orbid) += 1;
    }
}

void Geometry::orbOCC(VecI& vecu, VecI& vecd, VecI& occ) const {
    occ = VecI(getUnitOrbNum(),0);
    for(int i = 0; i < getOrbNum(); ++i) {
        if (vecu[i]) occ.at(orbs.at(i).orbid) += 1;
        if (vecd[i]) occ.at(orbs.at(i).orbid) += 1;
    }
}

void Geometry::orbOCC(VecI& vecu, VecI& vecd, VecI& occ, VecI& docc) const {
    occ = VecI(getUnitOrbNum(),0); 
    docc = VecI(getUnitOrbNum(),0);
    for(int i = 0; i < getOrbNum(); ++i) {
        if (vecu[i]) {
            occ.at(orbs.at(i).orbid) += 1;
            if (vecd[i]) docc.at(orbs.at(i).orbid) += 1;
        }
        if (vecd[i]) {
            occ.at(orbs.at(i).orbid) += 1;
        }
    }
}

void Geometry::orbOCC(idx_t repI, VecI& occ) const {
    occ = VecI(getUnitOrbNum(),0); 
    for(int i = 0; i < getOrbNum(); ++i) {
        if(bitTest(repI,i)) {
            occ.at(orbs.at(i).id) += 1;
        }
    }
}
void Geometry::orbOCC(pairIdx_t pairRepI, VecI& occ) const {
    occ = VecI(getUnitOrbNum(),0);
    for(int i = 0; i < getOrbNum(); ++i) {
        if (bitTest(pairRepI.first,i)) {
            occ.at(orbs.at(i).orbid) += 1;
        }
        if (bitTest(pairRepI.second,i)) {
            occ.at(orbs.at(i).orbid) += 1;
        }
    }
}
void Geometry::orbOCC(pairIdx_t pairRepI, VecI& occ, VecI& docc) const {
    occ = VecI(getUnitOrbNum(),0); docc = VecI(getUnitOrbNum(),0);
    for(int i = 0; i < getOrbNum(); ++i) {
        if(bitTest(pairRepI.first,i)) {
            occ.at(orbs.at(i).orbid) += 1;
            if(bitTest(pairRepI.second,i))docc.at(orbs.at(i).orbid) += 1;
        }
        if(bitTest(pairRepI.second,i))occ.at(orbs.at(i).orbid) += 1;
    }
}

/**
 * @brief check there is no same ORBITAL with same coordinates. 
 * This will ensure the search (ORBITAL, coord) --> id 
 * 
 * @return true passed check
 * @return false 
 */
bool Geometry::check( ) const {
    for (size_t i = 0; i < enlgOrbs.size(); ++i) {
        for (size_t j = i + 1; j < enlgOrbs.size(); ++j) {
            if (enlgOrbs[i].orb == enlgOrbs[j].orb) {
                bool cond = true;
                for (int r = 0; r < dim; ++r){
                    if(std::abs(enlgOrbs[i].coord[r]-enlgOrbs[j].coord[r])>tol) {cond = false; break;}
                } 
                if (cond) return false;
            }
        }
    }
    return true;
}

bool Geometry::crossBoundx(const Vec3d &coord) const {
    for(auto const& orb:orbs) {
        if(std::abs(coord.at(0)-orb.coord.at(0)) < tol) return false;
    }
    return true;
}

bool Geometry::crossBoundy(const Vec3d &coord) const {
    for(auto const& orb:orbs) {
        if(std::abs(coord.at(1)-orb.coord.at(1)) < tol) return false;
    }
    return true;    
}

void Geometry::genTransList( ) {
    Vec3d coordi, coordr, coordf;
    TransList.clear();
    std::vector<int> tmp;
    std::vector<cdouble> phasetmp;
    for (int r = 0; r < getSiteNum(); ++r) {
        tmp.clear();
        phasetmp.clear();
        coordr = getSiteR(r);
        for(int orbidi = 0; orbidi < getOrbNum(); ++orbidi) {
            coordi = getOrbR(orbidi);
            coordf = coordi + coordr;
            int orbidf;
            if (coordToOrbid(orbs.at(orbidi).orb, coordf, orbidf)) {
                cdouble tranPhase = 1.0;
                if(crossBoundx(coordf)) tranPhase *= std::exp(-CPLX_I*phase.at(0));
                if(crossBoundy(coordf)) tranPhase *= std::exp(-CPLX_I*phase.at(1));
                tmp.push_back(orbidf);
                phasetmp.push_back(tranPhase);
            } else {
                std::cout<<"translation position not found for orbid = "<<orbidi<<", transVecid = "<<r<<"\n";
                exit(1);
            }
        }
        TransList.push_back(tmp);
        TransPhaseList.push_back(phasetmp);
    }
}

bool Geometry::rotate(int orbid, int& orbidf) const {
    Vec3d coordi, coordr, coordrp, coordf;
    coordi = getOrbR(orbid);
    coordr = coordi - center;
    switch(PG){
        case PointGroup::D3: case PointGroup::C3:
            /*
                a1->a2,a2->-a1-a2
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1-a2) = -x2*a1 + (x1-x2)*a2
            */
            coordrp[0] = -coordr[1]; 
            coordrp[1] = coordr[0]-coordr[1]; 
            coordrp[2] = coordr[2];
            break;
        case PointGroup::D4: case PointGroup::D4m: case PointGroup::D4m5: case PointGroup::C4:
            /*
                a1->a2,a2->-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1) = -x2*a1 + x1*a2
            */
            coordrp[0] = -coordr[1]; 
            coordrp[1] = coordr[0]; 
            coordrp[2] = coordr[2];
            break;
        case PointGroup::D6: case PointGroup::C6:
            /*
                a1->a2,a2->a2-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(a2-a1) = -x2*a1 + (x1+x2)*a2
            */
            coordrp[0] = -coordr[1]; 
            coordrp[1] = coordr[0] + coordr[1]; 
            coordrp[2] = coordr[2];  
            break;
        case PointGroup::NONE:
            return false;
            break;
        default:
            return false;
            break;   
    }
    coordf = center + coordrp;
    return coordToOrbid(orbs.at(orbid).orb, coordf, orbidf);
}

Vec3d Geometry::rotate(const Vec3d &coordr) const {
    Vec3d coordrp;
    switch(PG){
        case PointGroup::D3: case PointGroup::C3:
            /*
                a1->a2,a2->-a1-a2
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1-a2) = -x2*a1 + (x1-x2)*a2
            */
            coordrp[0] = -coordr[1]; 
            coordrp[1] = coordr[0]-coordr[1]; 
            coordrp[2] = coordr[2];
            break;
        case PointGroup::D4: case PointGroup::D4m: case PointGroup::D4m5: case PointGroup::C4:
            /*
                a1->a2,a2->-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1) = -x2*a1 + x1*a2
            */
            coordrp[0] = -coordr[1]; 
            coordrp[1] = coordr[0]; 
            coordrp[2] = coordr[2];
            break;
        case PointGroup::D6: case PointGroup::C6:
            /*
                a1->a2,a2->a2-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(a2-a1) = -x2*a1 + (x1+x2)*a2
            */
            coordrp[0] = -coordr[1]; 
            coordrp[1] = coordr[0] + coordr[1]; 
            coordrp[2] = coordr[2];  
            break;
        case PointGroup::NONE:
            return coordr;
            break;
        default:
            return coordr;
            break;   
    }
    return coordrp;
}
bool Geometry::reflect(int orbid, int& orbidf) const {
    Vec3d coordi, coordr, coordrp, coordf;
    coordi = getOrbR(orbid);
    coordr = coordi - center;
    switch(PG){
        case PointGroup::D3:
            /*
                a1->a1,a2->-a1-a2
                x1*a1 + x2*a2 -> x1*a1 + x2*(-a1-a2) = (x1-x2)*a1 + (-x2)*a2
            */
            coordrp[0] = coordr[0]-coordr[1]; 
            coordrp[1] = -coordr[1]; 
            coordrp[2] = coordr[2];
            break;
        case PointGroup::D4: case PointGroup::D4m: case PointGroup::D4m5:
            /*
                a1->a1,a2->-a2
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1) = -x2*a1 + x1*a2
            */
            coordrp[0] = coordr[0]; 
            coordrp[1] = -coordr[1]; 
            coordrp[2] = coordr[2];
            break;
        case PointGroup::D6:
            /*
                a1->a1,a2->a1-a2
                x1*a1 + x2*a2 -> x1*a1 + x2*(a1-a2) = (x1+x2)*a1 + (-x2)*a2
            */
            coordrp[0] = coordr[0]+coordr[1]; 
            coordrp[1] = -coordr[1]; 
            coordrp[2] = coordr[2];  
            break;
        case PointGroup::NONE: case PointGroup::C3: case PointGroup::C4: case PointGroup::C6:
            return false;
            break;
        default:
            return false;
            break;   
    }
    coordf = center + coordrp;
    return coordToOrbid(orbs.at(orbid).orb, coordf, orbidf);
}
bool Geometry::mirror(int orbid, int& orbidf) const {
    Vec3d coordi, coordr, coordrp, coordf;
    coordi = getOrbR(orbid);
    coordr = coordi - center;
    switch(PG){
        case PointGroup::D4m5:
            /*
                a1->a1, a2->a2, a3->-a3
            */
            coordrp[0] = coordr[0]; 
            coordrp[1] = coordr[1]; 
            coordrp[2] = -coordr[2];
            break;
        default:
            return false;
            break;   
    }
    coordf = center + coordrp;
    return coordToOrbid(orbs.at(orbid).orb, coordf, orbidf);
}
void Geometry::genPGList(){
    int PGdeg;
    cdouble e,ez; // conjugate pair
    switch(PG){
        case PointGroup::D3: 
            PGdeg = 6; 
            e = std::exp(CPLX_I*2.0*PI/3.0);
            ez =  std::exp(-CPLX_I*2.0*PI/3.0);
            //                E,  R1, R2, z, zR1, zR2
            CharacterList = {{1,  1,  1,  1,  1,  1},//A1:R=1,z=1
                             {1,  1,  1, -1, -1, -1},//A2:R=1,z=-1
                             {1,  e, ez},//Ea:R=e
                             {1, ez, e} //Eb:R=ez
                            };
            break;
        case PointGroup::D4: 
            PGdeg = 8; 
            e = std::exp(CPLX_I*2.0*PI/4.0);
            ez =  std::exp(-CPLX_I*2.0*PI/4.0);
            //                E,  R1, R2, R3, z, zR1, zR2, zR3
            CharacterList = {{1,  1,  1,  1,  1,  1,  1,   1},//A1:R=1,z=1
                             {1,  1,  1,  1, -1, -1, -1,  -1},//A2:R=1,z=-1
                             {1, -1,  1, -1,  1, -1,  1,  -1},//B1:R=-1,z=1
                             {1, -1,  1, -1, -1,  1, -1,   1},//B2:R=-1,z=-1
                             {1,  e, -1, ez},//Ea:R=e,inv=-1
                             {1, ez, -1,  e} //Eb:R=ez,inv=-1
                            };
            break;
        case PointGroup::D4m: 
            PGdeg = 2;
        //                    E, Z*R
            CharacterList = {{1, 1},    //A1
                             {1, -1}};  //A2
            break;
        case PointGroup::D4m5: 
            PGdeg = 2;
        //                    E, M*Z*R
            CharacterList = {{1, 1},    //A1
                             {1, -1}};  //A2
            break;
        case PointGroup::D6: 
            PGdeg = 12;
            e = std::exp(CPLX_I*2.0*PI/6.0);
            ez =  std::exp(-CPLX_I*2.0*PI/6.0);
            //                E,  R1, R2, R3, R4, R5, z, zR1, zR2, zR3, zR4, zR5
            CharacterList = {{1,  1,  1,  1,  1,  1,  1,  1,  1,   1,   1,   1},//A1:R=1,z=1
                             {1,  1,  1,  1,  1,  1, -1, -1, -1,  -1,  -1,  -1},//A2:R=1,z=-1
                             {1, -1,  1, -1,  1, -1,  1, -1,  1,  -1,   1,  -1},//B1:R=-1,z=1
                             {1, -1,  1, -1,  1, -1, -1,  1, -1,   1,  -1,   1},//B2:R=-1,z=-1
                             {1,  e,-ez, -1, -e, ez},//E1a:R=e,inv=-1
                             {1, ez, -e, -1,-ez,  e},//E1b:R=e^5,inv=-1
                             {1, -ez,-e,  1,-ez, -e},//E2a:R=e^2,inv=1
                             {1, -e,-ez,  1, -e, -ez} //E2b:R=e^4,inv=1
                            };
            break;
        default: PGdeg = 0; break;
    }
    PGList.resize(PGdeg);
    int orbid, orbidf;
    switch (PG){
        case PointGroup::D3: case PointGroup::D4: case PointGroup::D6:
            for(orbid = 0; orbid < getOrbNum(); orbid++) PGList.at(0).push_back(orbid);    
            // rotation
            for(int i = 1; i < PGdeg/2; ++i){
                for(int j = 0; j < getOrbNum(); j++){
                    orbid = PGList.at(i-1).at(j);
                    if(rotate(orbid,orbidf)) PGList.at(i).push_back(orbidf);
                    else{std::cout<<"rotation of orbital:"<<orbid<<" not found!"<<std::endl; exit(1);}
                }
            }
            // reflection
            for(int i = PGdeg/2; i < PGdeg; ++i){
                for(int j = 0; j < getOrbNum(); j++){
                    orbid = PGList.at(i-PGdeg/2).at(j);
                    if(reflect(orbid,orbidf)) PGList.at(i).push_back(orbidf);
                    else{std::cout<<"reflection of orbital:"<<orbid<<" not found!"<<std::endl; exit(1);}
                }
            }
            break;
        case PointGroup::D4m:
            for(orbid = 0; orbid < getOrbNum(); orbid++) PGList.at(0).push_back(orbid);
            // reflection * rotation
            for(orbid = 0; orbid < getOrbNum(); orbid++){
                int orbtmp;
                if(rotate(orbid,orbtmp)) if(reflect(orbtmp,orbidf)) PGList.at(1).push_back(orbidf);
            }
            break;
        case PointGroup::D4m5:
            for(orbid = 0; orbid < getOrbNum(); orbid++) PGList.at(0).push_back(orbid);
            // reflection * rotation
            for(orbid = 0; orbid < getOrbNum(); orbid++){
                int orbtmp1, orbtmp2;
                if(rotate(orbid,orbtmp1)) if(reflect(orbtmp1,orbtmp2)) if(mirror(orbtmp2,orbidf)) PGList.at(1).push_back(orbidf);
            }
            break;
        default: break;
    }
}

void Geometry::construct(){
    if (PG==PointGroup::D4){
        if(getUnitOrbNum()==1) PG=PointGroup::D4;
        else if(getUnitOrbNum()==3) PG=PointGroup::D4m;
        else if(getUnitOrbNum()==5) PG=PointGroup::D4m5;
        else PG=PointGroup::NONE;
    }
    if (boundary.size() > 0) assert(!is_PBC);
    Norb = Nsite * unitSite.size() + boundary.size();
    Norb_enlg = is_PBC ? Norb * TranVecs.size() : Norb;
    assert(Nsite>0 and Norb>0 and Norb_enlg>=Norb);
    assert((int)xlist.size() == getSiteNum() and (int)ylist.size() == getSiteNum() and (int)zlist.size( )== getSiteNum());
    Vec3d vsite;
    Vec3d vorb;
    int id = 0;
    // construc Lattice and orbs
    for (int siteid = 0; siteid < getSiteNum(); ++siteid) {
        vsite = xlist.at(siteid) * a1 + ylist.at(siteid) * a2 + zlist.at(siteid) * a3; 
        Lattice.push_back(Site{siteid, vsite});
        for (int i = 0; i < (int)unitSite.size(); ++i){
            vorb = vsite + unitSite.at(i).coord;
            orbs.push_back(Orbital{unitSite.at(i).orb, unitSite.at(i).orbid, vorb});
            orbs.at(id).siteid = siteid;
            orbs.at(id).id = id;
            ++id;
        }
    }
    // ! siteid for boundary orb
    for(auto orb:boundary){
        orb.id = id;
        ++id;
        orbs.push_back(orb);
    }
    // construc enlarged orbs (periodic boundary condition)
    int count = 0;
    for (int R = 0; R < (int)TranVecs.size(); ++R) {
        for (int r = 0; r < (int)orbs.size(); ++r) {
            vorb = orbs.at(r).coord + TranVecs.at(R);
            enlgOrbs.push_back(Orbital{orbs.at(r).orb, orbs.at(r).orbid, vorb});
            enlgOrbs.at(count).id = orbs.at(r).id;
            enlgOrbs.at(count).siteid = orbs.at(r).siteid;
            ++count;
        }
    }
    assert((int)enlgOrbs.size() == Norb_enlg);
    // construct corresponding k-space lattice
    if (is_PBC) {
        assert((int)kxlist.size() == getSiteNum() and (int)kylist.size() == getSiteNum() and (int)kzlist.size() == getSiteNum());
        Vec3d vsite_phase;
        for (int siteid = 0; siteid < getSiteNum(); ++siteid){
            vsite = kxlist.at(siteid) * b10 + kylist.at(siteid) * b20 + kzlist.at(siteid) * b30;
            vsite_phase = vsite + dPhase;
            KLattice.push_back(Site{siteid, vsite_phase});
        }
        // generate transformation matrix for translation operation
        genTransList();
    } 
    if (PG != PointGroup::NONE) genPGList();
    // return resources no longer needed
    xlist.clear(); ylist.clear(); zlist.clear();
    kxlist.clear(); kylist.clear(); kzlist.clear();
}

void Geometry::printLattice() const {
    // print Lattice
    std::cout<<"Lattice Unit Cell:\n";
    for (int i = 0; i < Nsite; i ++){
        std::cout<<"siteID:"<<Lattice[i].id<<", coord:"<<Lattice[i].coord<<"\n";
    }
    std::cout<<"\n";
}
void Geometry::printOrbs() const {
    std::cout<<"Lattice Orbitals:\n";
    for (int i = 0; i < Norb; i ++){
        std::cout<<"orbID:"<<orbs[i].id<<", from site"<<orbs[i].siteid<<", coord:"<<orbs[i].coord<<"\n";
    }
    std::cout<<"\n";
}
void Geometry::printKLattice() const {
    // print kLattice
    std::cout<<"k-Lattice:\n";
    for (int i = 0; i < Nsite; i ++){
        std::cout<<"siteID:"<<KLattice[i].id<<", coord:"<<KLattice[i].coord<<"\n";
    }
    std::cout<<"\n";
}
void Geometry::printTrans() const {
    std::cout<<"TransList:\n";
    for (int r = 0; r < Nsite; r++){
        std::cout<<"r("<<r<<"):\n";
        for (int i = 0; i < Norb; ++i){
            std::cout<<getOrbTran(r,i)<<" "<<"phase:"<<getOrbTranPhase(r,i)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}
void Geometry::printPG() const {
    std::cout<<"Point Group Transformation List:\n";
    for (int r = 0; r < (int)PGList.size(); r++){
        std::cout<<"PG("<<r<<"):\n";
        for (int i = 0; i < (int)PGList[r].size(); ++i){
            std::cout<<getOrbPG(r,i)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}
void Geometry::print() const {
    printLattice();
    printOrbs();
    printTrans();
    printKLattice();
    printPG();
}

/*
    **********************
    * Triangular Lattice *
    **********************
*/

//constructor with T x D6 Symm
TriAngLattice::TriAngLattice(int numSites, bool PBC){
    type = LATTICE::TRIANGULAR;
    PG = PointGroup::D6;
    is_PBC = PBC;
    Nsite = numSites;
    name = "TriAng_D6_N"+std::to_string(Nsite);
    ax = Vec3d {1.0, 0.0, 0.0};
    ay = Vec3d {0.5, std::sqrt(3.0)/2.0, 0.0};
    az = Vec3d {0.0, 0.0, 1.0};
    bx = Vec3d {2*PI, -2*PI/std::sqrt(3.0), 0.0};
    by = Vec3d {0.0, 4*PI/std::sqrt(3.0), 0.0};
    bz = Vec3d {0.0, 0.0, 2*PI};
    // unitSite.push_back(Orbital{0,0,VecD{0.0,0.0,0.0},ORBITAL::SINGLE});
    switch(numSites){
        case 9:{
            R1 = 3.0 * a1 - 3.0 * a2;
            R2 = 3.0 * a1 + 0.0 * a2;
            b10 = 1.0/3.0 * b1 + 1.0/3.0 * b2;
            b20 = 0.0 * b1 + 1.0/3.0 * b2;
            xlist = VecD {+0, +1, -1, 0, 1, -1, 0, 1, -1};
            ylist = VecD {-1, -1, +0, 0, 0, +1, 1, 1, +2};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, +0.0, +1.0, +2.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, -1.0, -1.0, -1.0};
            break;
        }
        case 12:{
            R1 = 4.0 * a1 - 2.0 * a2;
            R2 = 2.0 * a1 + 2.0 * a2;
            b10 = 1.0/6.0 * b1 - 1.0/6.0 * b2;
            b20 = 1.0/6.0 * b1 + 2.0/6.0 * b2;
            xlist = VecD {+0, +1, +2, -1, 0, 1, 2, -1, 0, 1, -1, 0};
            ylist = VecD {-1, -1, -1, +0, 0, 0, 0, +1, 1, 1, +2, 2};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, -1.0, -1.0};
            break;
        }
        //TODO:fix 21-sites point group
        case 21:{
            R1 = 5.0 * a1 - 4.0 * a2;
            R2 = 4.0 * a1 + 1.0 * a2;
            b10 = 1.0/21.0 * b1 - 4.0/21.0 * b2;
            b20 = 4.0/21.0 * b1 + 5.0/21.0 * b2;
            xlist = VecD {+0, +1, +2, -1, +0, +1, +2, -2, -1, 0, 1, 2, -2, -1, 0, 1, -2, -1, 0, 1, -2};
            ylist = VecD {-2, -2, -2, -1, -1, -1, -1, +0, +0, 0, 0, 0, +1, +1, 1, 1, +2, +2, 2, 2, +3};
            kxlist = VecD {0.0, 1.0, 2.0, -2.0, -1.0, -1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0, 2.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, 0.0, +0.0, +0.0, +1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 3.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0};
            break;
        }
        case 27:{
            R1 = 6.0 * a1 - 3.0 * a2;
            R2 = 3.0 * a1 + 3.0 * a2;
            b10 = 1.0/9.0 * b1 - 1.0/9.0 * b2;
            b20 = 1.0/9.0 * b1 + 2.0/9.0 * b2;
            xlist = VecD {+0, +1, +2, +3, -1, +0, +1, +2, +3, -2, -1, 0, 1, 2, 3, -2, -1, 0, 1, 2, -2, -1, 0, 1, -2, -1, 0};
            ylist = VecD {-2, -2, -2, -2, -1, -1, -1, -1, -1, +0, +0, 0, 0, 0, 0, +1, +1, 1, 1, 1, +2, +2, 2, 2, +3, +3, 3};
            kxlist = VecD {0.0, 1.0, 2.0, -2.0, -1.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, -1.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, 0.0, +0.0, +0.0, +1.0, +1.0, 1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0};
            break;
        }
        case 36:{
            R1 = 6.0 * a1 - 0.0 * a2;
            R2 = 0.0 * a1 + 6.0 * a2;
            b10 = 1.0/6.0 * b1 + 0.0 * b2;
            b20 = 0.0 * b1 + 1.0/6.0 * b2;
            xlist = VecD {+0, +1, +2, -2, -1, +0, +1, +2, +3, -2, -1, +0, +1, +2, +3, -3, -2, -1, +0, +1, +2, -3, -2, -1, +0, +1, +2, -4, -3, -2, -1, 0, 1, -3, -2, -1};
            ylist = VecD {-3, -3, -3, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, +0, +0 ,+0 ,+0, +0, +0, +1, +1, +1, +1, +1, +1, +2, +2, +2, +2, 2, 2, +3, +3, +3};
            kxlist = VecD {0.0, 1.0, 2.0, 3.0, -2.0, -1.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 2.0, -3.0, -2.0, -1.0, +0.0, +1.0, +2.0, -3.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0};
            kylist = VecD {0.0, 0.0, 0.0, 0.0, +0.0, +0.0, +1.0, +1.0, 1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -3.0, -3.0};
            break;
        }
        default: 
            std::cout<<"N="<<Nsite<<" must be one of the values: 9, 12, 21, 27, 36!"<<std::endl;
            exit(1);
            break;
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    Vec3d vtmp;
    vtmp.fill(0.0);
    TranVecs.push_back(vtmp);
    if(is_PBC) {
        VecD c1s {1.0, 0.0, -1.0, -1.0,  0.0,  1.0};
        VecD c2s {0.0, 1.0,  1.0,  0.0, -1.0, -1.0};
        for (int i = 0; i < (int)c1s.size(); ++i) {
            vtmp = c1s.at(i) * R1 + c2s.at(i) * R2;
            TranVecs.push_back(vtmp);
        }
    }
    
    // construct();
};

//constructor with T without D6 (possible C6?)
TriAngLattice::TriAngLattice(int N1, int N2, bool PBC) {
    type = LATTICE::TRIANGULAR;
    // default PBC condition
    is_PBC = PBC;
    Nsite = N1 * N2;
    name = "TriAng"+std::to_string(N1)+"x"+std::to_string(N2);

    ax = Vec3d {1.0, 0.0, 0.0};
    ay = Vec3d {0.5, std::sqrt(3.0)/2.0, 0.0};
    az = Vec3d {0.0, 0.0, 1.0};
    bx = Vec3d {2*PI, -2*PI/std::sqrt(3.0), 0.0};
    by = Vec3d {0.0, 4*PI/std::sqrt(3.0), 0.0};
    bz = Vec3d {0.0, 0.0, 2*PI};
    
    // unitSite.push_back(Orbital{0,0,VecD{0.0,0.0,0.0},ORBITAL::SINGLE});

    R1 = (double)N1 * a1;
    R2 = (double)N2 * a2;
    b10 = 1.0/(double)N1 * b1;
    b20 = 1.0/(double)N2 * b2;
    for (int y = 0; y < N2; ++y) {
        for (int x = 0; x < N1; ++x) {
            xlist.push_back(x);
            ylist.push_back(y);
            kxlist.push_back(x);
            kylist.push_back(y);
        }
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    Vec3d vtmp;
    vtmp.fill(0.0);
    TranVecs.push_back(vtmp);
    if(is_PBC){
        VecD c1s, c2s;
        if (N1==1 and N2>1){
            c1s = VecD{0.0,  0.0};
            c2s = VecD{1.0, -1.0};
        }
        else if (N1>1 and N2==1){
            c1s = VecD{1.0, -1.0};
            c2s = VecD{0.0,  0.0};
        }
        else if (N1>1 and N2>1){
            c1s = VecD{1.0, 1.0, 0.0, -1.0, -1.0, -1.0,  0.0,  1.0};
            c2s = VecD{0.0, 1.0, 1.0,  1.0,  0.0, -1.0, -1.0, -1.0};
        }
        for (int i = 0; i < (int)c1s.size(); ++i) {
            vtmp = c1s.at(i) * R1 + c2s.at(i) * R2;
            TranVecs.push_back(vtmp);
        }
    }
    
    // construct();
};


/*
    ******************
    * Square Lattice *
    ******************
*/
SquareLattice::SquareLattice(int N1, int N2, bool PBC, bool TBC, double phase_x, double phase_y){
    type = LATTICE::SQUARE;
    if(N1==N2) {
        PG = PointGroup::D4;
        center = {double(N1-1)/2.0, double(N2-1)/2.0, 0.0};
    }
    is_PBC = PBC;
    is_TBC = TBC;
    phase = {phase_x, phase_y, 0.0};
    dPhase = {phase_x/N1, phase_y/N2, 0.0};
    Nsite = N1 * N2;
    name = "Square"+std::to_string(N1)+"x"+std::to_string(N2);
    ax = {1.0, 0.0, 0.0};
    ay = {0.0, 1.0, 0.0};
    az = {0.0, 0.0, 2.429/1.89};
    bx = {2*PI, 0.0, 0.0};
    by = {0.0, 2*PI, 0.0};
    bz = {0.0, 0.0, 2*PI};

    R1 = (double)N1 * a1;
    R2 = (double)N2 * a2;
    b10 = 1.0/(double)N1 * b1;
    b20 = 1.0/(double)N2 * b2;
    for (int y = 0; y < N2; y++){
        for (int x = 0; x < N1; x++){
            xlist.push_back(x);
            ylist.push_back(y);
            kxlist.push_back(x);
            kylist.push_back(y);
        }
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    Vec3d vtmp;
    vtmp.fill(0.0);
    TranVecs.push_back(vtmp);
    if(is_PBC){
        VecD c1s, c2s;
        if (N1==1 and N2>1){
            c1s = VecD{0.0,  0.0};
            c2s = VecD{1.0, -1.0};
        }
        else if (N1>1 and N2==1){
            c1s = VecD{1.0, -1.0};
            c2s = VecD{0.0,  0.0};
        }
        else if (N1>1 and N2>1){
            c1s = VecD{1.0, 1.0, 0.0, -1.0, -1.0, -1.0,  0.0,  1.0};
            c2s = VecD{0.0, 1.0, 1.0,  1.0,  0.0, -1.0, -1.0, -1.0};
        }    
        for (int i = 0; i < (int)c1s.size(); ++i) {
            vtmp = c1s.at(i) * R1 + c2s.at(i) * R2;
            TranVecs.push_back(vtmp);
        }
    }
}

SquareLattice::SquareLattice(int N, bool PBC){
    type = LATTICE::SQUARE;
    PG = PointGroup::D4;
    is_PBC = PBC;
    Nsite = N;
    name = "Square"+std::to_string(N);
    ax = {1.0, 0.0, 0.0};
    ay = {0.0, 1.0, 0.0};
    az = {0.0, 0.0, 2.429/1.89};
    bx = {2*PI, 0.0, 0.0};
    by = {0.0, 2*PI, 0.0};
    bz = {0.0, 0.0, 2*PI};

    switch(Nsite){
        case 8:{
            R1 = 2.0 * a1 - 2.0 * a2;
            R2 = 2.0 * a1 + 2.0 * a2;
            b10 = 1.0/4.0 * b1 - 1.0/4.0 * b2;
            b20 = 1.0/4.0 * b2 + 1.0/4.0 * b2;
            xlist = VecD {+0, +1, +2, -1, +0, +1, +0, +1};
            ylist = VecD {+0, +0, +0, +0, +1, +1, -1, -1};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, +0.0, +0.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, +2.0, -1.0};
            break;
        }
        default: 
            std::cout<<"N="<<Nsite<<" must be one of the values: 8!\n";
            exit(1);
            break;
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    Vec3d vtmp;
    vtmp.fill(0.0);
    TranVecs.push_back(vtmp);
    if(is_PBC){
        auto c1s = VecD{1.0, 1.0, 0.0, -1.0, -1.0, -1.0,  0.0,  1.0};
        auto c2s = VecD{0.0, 1.0, 1.0,  1.0,  0.0, -1.0, -1.0, -1.0};
        for (int i = 0; i < (int)c1s.size(); ++i) {
            vtmp = c1s.at(i) * R1 + c2s.at(i) * R2;
            TranVecs.push_back(vtmp);
        } 
    }
}