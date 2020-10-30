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

bool Geometry::coordToOrbid(double* coord, int &orbid, double tol) const {
    bool cond;
    for (auto it = enlgOrbs.begin(); it != enlgOrbs.end(); it++){
        cond = true;
        for (int i = 0; i < dim; i++){
            if(std::abs(coord[i]-(*it).coord[i])>tol) {cond = false; break;}
        }
        if (cond) {orbid = (*it).id; return cond;}
        else{continue;}
    }
    return false;
}

void Geometry::genTransList(){
    VecD coordi(3), coordr(3), coordf(3);
    TransList.clear();
    std::vector<int> tmp;
    for (int r = 0; r < getSiteNum(); r++){
        tmp.clear();
        getSiteR(r,coordr.data());
        for(int orbidi = 0; orbidi < getOrbNum(); orbidi++){
            getOrbR(orbidi,coordi.data());
            vecXAdd(1.0, coordi.data(), 1.0, coordr.data(), coordf.data(), dim);
            int orbidf;
            if (coordToOrbid(coordf.data(), orbidf)){
                tmp.push_back(orbidf);
            }else{
                std::cout<<"translation position not found for orbid = "<<orbidi<<", transVecid = "<<r<<std::endl;
                exit(1);
            }
        }
        TransList.push_back(tmp);
    }
}

bool Geometry::rotate(int orbid, int& orbidf) const {
    VecD coordi(3), coordr(3), coordrp(3), coordf(3);
    getOrbR(orbid,coordi.data());
    vecXAdd(1.0, coordi.data(), -1.0, center.data(), coordr.data(), dim);
    switch(PG){
        case PointGroup::D3: case PointGroup::C3:
            /*
                a1->a2,a2->-a1-a2
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1-a2) = -x2*a1 + (x1-x2)*a2
            */
            coordrp[0] = -coordr[1]; coordrp[1] = coordr[0]-coordr[1]; coordrp[2] = coordr[2];
            break;
        case PointGroup::D4: case PointGroup::D4m: case PointGroup::D4m5: case PointGroup::C4:
            /*
                a1->a2,a2->-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1) = -x2*a1 + x1*a2
            */
            coordrp[0] = -coordr[1]; coordrp[1] = coordr[0]; coordrp[2] = coordr[2];
            break;
        case PointGroup::D6: case PointGroup::C6:
            /*
                a1->a2,a2->a2-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(a2-a1) = -x2*a1 + (x1+x2)*a2
            */
            coordrp[0] = -coordr[1]; coordrp[1] = coordr[0] + coordr[1]; coordrp[2] = coordr[2];  
            break;
        case PointGroup::NONE:
            return false;
            break;
        default:
            return false;
            break;   
    }
    vecXAdd(1.0, center.data(), 1.0, coordrp.data(), coordf.data(), dim);
    return coordToOrbid(coordf.data(), orbidf);
}

VecD Geometry::rotate(VecD coordr) const {
    VecD coordrp(3);
    switch(PG){
        case PointGroup::D3: case PointGroup::C3:
            /*
                a1->a2,a2->-a1-a2
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1-a2) = -x2*a1 + (x1-x2)*a2
            */
            coordrp[0] = -coordr[1]; coordrp[1] = coordr[0]-coordr[1]; coordrp[2] = coordr[2];
            break;
        case PointGroup::D4: case PointGroup::D4m: case PointGroup::D4m5: case PointGroup::C4:
            /*
                a1->a2,a2->-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1) = -x2*a1 + x1*a2
            */
            coordrp[0] = -coordr[1]; coordrp[1] = coordr[0]; coordrp[2] = coordr[2];
            break;
        case PointGroup::D6: case PointGroup::C6:
            /*
                a1->a2,a2->a2-a1
                x1*a1 + x2*a2 -> x1*a2 + x2*(a2-a1) = -x2*a1 + (x1+x2)*a2
            */
            coordrp[0] = -coordr[1]; coordrp[1] = coordr[0] + coordr[1]; coordrp[2] = coordr[2];  
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
    VecD coordi(3), coordr(3), coordrp(3), coordf(3);
    getOrbR(orbid,coordi.data());
    vecXAdd(1.0, coordi.data(), -1.0, center.data(), coordr.data(), dim);
    switch(PG){
        case PointGroup::D3:
            /*
                a1->a1,a2->-a1-a2
                x1*a1 + x2*a2 -> x1*a1 + x2*(-a1-a2) = (x1-x2)*a1 + (-x2)*a2
            */
            coordrp[0] = coordr[0]-coordr[1]; coordrp[1] = -coordr[1]; coordrp[2] = coordr[2];
            break;
        case PointGroup::D4: case PointGroup::D4m: case PointGroup::D4m5:
            /*
                a1->a1,a2->-a2
                x1*a1 + x2*a2 -> x1*a2 + x2*(-a1) = -x2*a1 + x1*a2
            */
            coordrp[0] = coordr[0]; coordrp[1] = -coordr[1]; coordrp[2] = coordr[2];
            break;
        case PointGroup::D6:
            /*
                a1->a1,a2->a1-a2
                x1*a1 + x2*a2 -> x1*a1 + x2*(a1-a2) = (x1+x2)*a1 + (-x2)*a2
            */
            coordrp[0] = coordr[0]+coordr[1]; coordrp[1] = -coordr[1]; coordrp[2] = coordr[2];  
            break;
        case PointGroup::NONE: case PointGroup::C3: case PointGroup::C4: case PointGroup::C6:
            return false;
            break;
        default:
            return false;
            break;   
    }
    vecXAdd(1.0, center.data(), 1.0, coordrp.data(), coordf.data(), dim);
    return coordToOrbid(coordf.data(), orbidf);
}
bool Geometry::mirror(int orbid, int& orbidf) const {
    VecD coordi(3), coordr(3), coordrp(3), coordf(3);
    getOrbR(orbid,coordi.data());
    vecXAdd(1.0, coordi.data(), -1.0, center.data(), coordr.data(), dim);
    switch(PG){
        case PointGroup::D4m5:
            /*
                a1->a1, a2->a2, a3->-a3
            */
            coordrp[0] = coordr[0]; coordrp[1] = coordr[1]; coordrp[2] = -coordr[2];
            break;
        default:
            return false;
            break;   
    }
    vecXAdd(1.0, center.data(), 1.0, coordrp.data(), coordf.data(), dim);
    return coordToOrbid(coordf.data(), orbidf);
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
            for(int i = 1; i < PGdeg/2; i++){
                for(int j = 0; j < getOrbNum(); j++){
                    orbid = PGList.at(i-1).at(j);
                    if(rotate(orbid,orbidf)) PGList.at(i).push_back(orbidf);
                    else{std::cout<<"rotation of orbital:"<<orbid<<" not found!"<<std::endl; exit(1);}
                }
            }
            // reflection
            for(int i = PGdeg/2; i < PGdeg; i++){
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
    if(boundary.size()>0)assert(!is_PBC);
    Norb = Nsite * unitSite.size()+boundary.size();
    Norb_enlg = is_PBC?Norb*TranVecs.size():Norb;
    assert(Nsite>0 and Norb>0 and Norb_enlg>=Norb);
    assert(xlist.size()==getSiteNum() and ylist.size()==getSiteNum() and zlist.size()==getSiteNum());
    VecD vsite(3,0.0);
    VecD vorb(3,0.0);
    int id = 0;
    // construc Lattice and orbs
    for (int siteid = 0; siteid < getSiteNum(); siteid++){
        vecXAdd(xlist.at(siteid), a1.data(), ylist.at(siteid), a2.data(), zlist.at(siteid), a3.data(), vsite.data(), getDim());
        Lattice.push_back(Site{siteid, vsite});
        for (int i = 0; i < unitSite.size(); i++){
            vecXAdd(1.0, vsite.data(), 1.0, unitSite.at(i).coord.data(), vorb.data(), getDim());
            orbs.push_back(Orbital{unitSite.at(i).orb,unitSite.at(i).orbid,vorb});
            orbs.at(id).siteid = siteid;
            orbs.at(id).id = id;
            id++;
        }
    }
    for(auto orb:boundary){
        orb.id = id;
        id++;
        orbs.push_back(orb);
    }
    // construc enlarged orbs (periodic boundary condition)
    int count = 0;
    for (int R = 0; R < TranVecs.size(); R++){
        for (int r = 0; r < orbs.size(); r++){
            vecXAdd(1.0, TranVecs.at(R).data(), 1.0, orbs.at(r).coord.data(), vorb.data(), getDim());
            enlgOrbs.push_back(Orbital{orbs.at(r).orb, orbs.at(r).orbid, vorb});
            enlgOrbs.at(count).id = orbs.at(r).id;
            enlgOrbs.at(count).siteid = orbs.at(r).siteid;
            count++;
        }
    }
    assert(enlgOrbs.size()==Norb_enlg);
    // construct corresponding k-space lattice
    if(is_PBC){
        assert(kxlist.size()==getSiteNum() and kylist.size()==getSiteNum() and kzlist.size()==getSiteNum());
        for (int siteid = 0; siteid < getSiteNum(); siteid++){
            vecXAdd(kxlist.at(siteid), b10.data(), kylist.at(siteid), b20.data(), kzlist.at(siteid), b30.data(), vsite.data(), getDim());
            KLattice.push_back(Site{siteid, vsite});
        }
        // generate transformation matrix for translation operation
        genTransList();
    } 
    if(PG != PointGroup::NONE) genPGList();
    // return resources no longer needed
    xlist.clear(); ylist.clear(); zlist.clear();
    kxlist.clear(); kylist.clear(); kzlist.clear();
}

void Geometry::printLattice() const {
    // print Lattice
    std::cout<<"Lattice Unit Cell:"<<std::endl;
    for (int i = 0; i < Nsite; i ++){
        std::cout<<"siteID:"<<Lattice[i].id<<", coord:["<<Lattice[i].coord[0]<<", "<<Lattice[i].coord[1]<<", "<<Lattice[i].coord[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::printOrbs() const {
    std::cout<<"Lattice Orbitals:"<<std::endl;
    for (int i = 0; i < Norb; i ++){
        std::cout<<"orbID:"<<orbs[i].id<<", from site"<<orbs[i].siteid<<", coord:["<<orbs[i].coord[0]<<", "<<orbs[i].coord[1]<<", "<<orbs[i].coord[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::printKLattice() const {
    // print kLattice
    std::cout<<"k-Lattice:"<<std::endl;
    for (int i = 0; i < Nsite; i ++){
        std::cout<<"siteID:"<<KLattice[i].id<<", coord:["<<KLattice[i].coord[0]<<", "<<KLattice[i].coord[1]<<", "<<KLattice[i].coord[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::printTrans() const {
    std::cout<<"TransList:"<<std::endl;
    for (int r = 0; r < Nsite; r++){
        std::cout<<"r("<<r<<"):"<<std::endl;
        for (int i = 0; i < Norb; i++){
            std::cout<<getOrbTran(r,i)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}
void Geometry::printPG() const {
    std::cout<<"Point Group Transformation List:"<<std::endl;
    for (int r = 0; r < PGList.size(); r++){
        std::cout<<"PG("<<r<<"):"<<std::endl;
        for (int i = 0; i < PGList[r].size(); i++){
            std::cout<<getOrbPG(r,i)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
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
    PG = PointGroup::D6;
    is_PBC = PBC;
    Nsite = numSites;
    name = "TriAng_D6_N"+std::to_string(Nsite);
    ax = VecD {1.0, 0.0, 0.0};
    ay = VecD {0.5, std::sqrt(3.0)/2.0, 0.0};
    az = VecD {0.0, 0.0, 1.0};
    bx = VecD {2*PI, -2*PI/std::sqrt(3.0), 0.0};
    by = VecD {0.0, 4*PI/std::sqrt(3.0), 0.0};
    bz = VecD {0.0, 0.0, 2*PI};
    // unitSite.push_back(Orbital{0,0,VecD{0.0,0.0,0.0},ORBITAL::SINGLE});
    switch(numSites){
        case 9:{
            vecXAdd(3.0, a1.data(), -3.0, a2.data(), R1.data(), getDim());
            vecXAdd(3.0, a1.data(), 0.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/3.0, b1.data(), 1.0/3.0, b2.data(), b10.data(), getDim());
            vecXAdd(0.0, b1.data(), 1.0/3.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, -1, 0, 1, -1, 0, 1, -1};
            ylist = VecD {-1, -1, +0, 0, 0, +1, 1, 1, +2};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, +0.0, +1.0, +2.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, -1.0, -1.0, -1.0};
            break;
        }
        case 12:{
            vecXAdd(4.0, a1.data(), -2.0, a2.data(), R1.data(), getDim());
            vecXAdd(2.0, a1.data(), 2.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/6.0, b1.data(), -1.0/6.0, b2.data(), b10.data(), getDim());
            vecXAdd(1.0/6.0, b1.data(), 2.0/6.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, -1, 0, 1, 2, -1, 0, 1, -1, 0};
            ylist = VecD {-1, -1, -1, +0, 0, 0, 0, +1, 1, 1, +2, 2};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, -1.0, -1.0};
            break;
        }
        case 21:{
            vecXAdd(5.0, a1.data(), -4.0, a2.data(), R1.data(), getDim());
            vecXAdd(4.0, a1.data(), 1.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/21.0, b1.data(), -4.0/21.0, b2.data(), b10.data(), getDim());
            vecXAdd(4.0/21.0, b1.data(), 5.0/21.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, -1, +0, +1, +2, -2, -1, 0, 1, 2, -2, -1, 0, 1, -2, -1, 0, 1, -2};
            ylist = VecD {-2, -2, -2, -1, -1, -1, -1, +0, +0, 0, 0, 0, +1, +1, 1, 1, +2, +2, 2, 2, +3};
            kxlist = VecD {0.0, 1.0, 2.0, -2.0, -1.0, -1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0, 2.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, 0.0, +0.0, +0.0, +1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 3.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0};
            break;
        }
        case 27:{
            vecXAdd(6.0, a1.data(), -3.0, a2.data(), R1.data(), getDim());
            vecXAdd(3.0, a1.data(), 3.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/9.0, b1.data(), -1.0/9.0, b2.data(), b10.data(), getDim());
            vecXAdd(1.0/9.0, b1.data(), 2.0/9.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, +3, -1, +0, +1, +2, +3, -2, -1, 0, 1, 2, 3, -2, -1, 0, 1, 2, -2, -1, 0, 1, -2, -1, 0};
            ylist = VecD {-2, -2, -2, -2, -1, -1, -1, -1, -1, +0, +0, 0, 0, 0, 0, +1, +1, 1, 1, 1, +2, +2, 2, 2, +3, +3, 3};
            kxlist = VecD {0.0, 1.0, 2.0, -2.0, -1.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, -1.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, -2.0, -1.0, +0.0, +1.0, -2.0, -1.0, +0.0};
            kylist = VecD {0.0, 0.0, 0.0, +0.0, +0.0, +1.0, +1.0, 1.0, 1.0, 1.0, 1.0, +2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0};
            break;
        }
        case 36:{
            vecXAdd(6.0, a1.data(), 0.0, a2.data(), R1.data(), getDim());
            vecXAdd(0.0, a1.data(), 6.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/6.0, b1.data(), 0.0, b2.data(), b10.data(), getDim());
            vecXAdd(0.0, b1.data(), 1.0/6.0, b2.data(), b20.data(), getDim());
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
    VecD vtmp(3);
    vecXAdd(0.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    if(is_PBC){
        vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), 1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }
    
    // construct();
};

//constructor with T without D6 (possible C6?)
TriAngLattice::TriAngLattice(int N1, int N2, bool PBC){
    // default PBC condition
    is_PBC = PBC;
    Nsite = N1 * N2;
    name = "TriAng"+std::to_string(N1)+"x"+std::to_string(N2);

    ax = VecD {1.0, 0.0, 0.0};
    ay = VecD {0.5, std::sqrt(3.0)/2.0, 0.0};
    az = VecD {0.0, 0.0, 1.0};
    bx = VecD {2*PI, -2*PI/std::sqrt(3.0), 0.0};
    by = VecD {0.0, 4*PI/std::sqrt(3.0), 0.0};
    bz = VecD {0.0, 0.0, 2*PI};
    
    // unitSite.push_back(Orbital{0,0,VecD{0.0,0.0,0.0},ORBITAL::SINGLE});

    vecXAdd((double)N1, a1.data(), 0.0, a2.data(), R1.data(), getDim());
    vecXAdd(0.0, a1.data(), (double)N2, a2.data(), R2.data(), getDim());
    vecXAdd(1.0/(double)N1, b1.data(), 0.0, b2.data(), b10.data(), getDim());
    vecXAdd(0.0, b1.data(), 1.0/(double)N2, b2.data(), b20.data(), getDim());
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
    VecD vtmp(3);
    vecXAdd(0.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    if(is_PBC){
        if (N1==1 and N2>1){
            vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(0.0, R1.data(), -1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        }
        else if (N1>1 and N2==1){
            vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        }
        else if (N1>1 and N2>1){
            vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(1.0, R1.data(), 1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(0.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        }
    }
    
    // construct();
};


/*
    ******************
    * Square Lattice *
    ******************
*/
SquareLattice::SquareLattice(int N1, int N2, bool PBC){
    if(N1==N2) {
        PG = PointGroup::D4;
        center = {double(N1-1)/2.0, double(N2-1)/2.0, 0.0};
    }
    is_PBC = PBC;
    Nsite = N1 * N2;
    name = "Square"+std::to_string(N1)+"x"+std::to_string(N2);
    ax = VecD {1.0, 0.0, 0.0};
    ay = VecD {0.0, 1.0, 0.0};
    az = VecD {0.0, 0.0, 1.0};
    bx = VecD {2*PI, 0.0, 0.0};
    by = VecD {0.0, 2*PI, 0.0};
    bz = VecD {0.0, 0.0, 2*PI};

    vecXAdd((double)N1, a1.data(), 0.0, a2.data(), R1.data(), getDim());
    vecXAdd(0.0, a1.data(), (double)N2, a2.data(), R2.data(), getDim());
    vecXAdd(1.0/(double)N1, b1.data(), 0.0, b2.data(), b10.data(), getDim());
    vecXAdd(0.0, b1.data(), 1.0/(double)N2, b2.data(), b20.data(), getDim());
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
    VecD vtmp(3);
    vecXAdd(0.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    if(is_PBC){
        if (N1==1 and N2>1){
            vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(0.0, R1.data(), -1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        }
        else if (N1>1 and N2==1){
            vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        }
        else if (N1>1 and N2>1){
            vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(1.0, R1.data(), 1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(-1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(0.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
            vecXAdd(1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        }    
    }
}

SquareLattice::SquareLattice(int N, bool PBC){
    PG = PointGroup::D4;
    is_PBC = PBC;
    Nsite = N;
    name = "Square"+std::to_string(N);
    ax = VecD {1.0, 0.0, 0.0};
    ay = VecD {0.0, 1.0, 0.0};
    az = VecD {0.0, 0.0, 1.0};
    bx = VecD {2*PI, 0.0, 0.0};
    by = VecD {0.0, 2*PI, 0.0};
    bz = VecD {0.0, 0.0, 2*PI};

    switch(Nsite){
        case 8:{
            vecXAdd(2.0, a1.data(), -2.0, a2.data(), R1.data(), getDim());
            vecXAdd(2.0, a1.data(), 2.0, a2.data(), R2.data(), getDim());
            vecXAdd(1.0/4.0, b1.data(),-1.0/4.0, b2.data(), b10.data(), getDim());
            vecXAdd(1.0/4.0, b1.data(), 1.0/4.0, b2.data(), b20.data(), getDim());
            xlist = VecD {+0, +1, +2, -1, +0, +1, +0, +1};
            ylist = VecD {+0, +0, +0, +0, +1, +1, -1, -1};
            kxlist = VecD {0.0, 1.0, -1.0, -1.0, 0.0, 1.0, +0.0, +0.0};
            kylist = VecD {0.0, 0.0, +0.0, +1.0, 1.0, 1.0, +2.0, -1.0};
            break;
        }
        default: 
            std::cout<<"N="<<Nsite<<" must be one of the values: 8!"<<std::endl;
            exit(1);
            break;
    }
    zlist.resize(Nsite,0.0);
    kzlist.resize(Nsite,0.0);
    VecD vtmp(3);
    vecXAdd(0.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    if(is_PBC){
        vecXAdd(1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(1.0, R1.data(), 1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 1.0 , R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), 0.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(-1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(0.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
        vecXAdd(1.0, R1.data(), -1.0, R2.data(), vtmp.data(), getDim()); TranVecs.push_back(vtmp);
    }
}