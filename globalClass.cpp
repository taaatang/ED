//
//  globalClass.cpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "globalClass.hpp"


/*
   ***************************
   *   Geometry Base Class   *
   ***************************
*/
Geometry::Geometry(){};

Geometry::~Geometry(){};

std::string Geometry::getName(){
    return name;
}

int Geometry::getSiteNum(){
    return N;
}

int* Geometry::indToCoord(int i){
    /*
     return the coord of site i
     */
    if (i<N and i>=0) return lattice_[i].coord;
    std::cerr<<"Index is out side the range [0,N)"<<std::endl;
    exit(EXIT_FAILURE);
};

double* Geometry::indToCoordxy(int i){
    /*
     return the real space coord of site i
     */
    if (i<N and i>=0) return lattice_[i].coordxy;
    std::cerr<<"Index is out side the range [0,N)"<<std::endl;
    exit(EXIT_FAILURE);
};

bool Geometry::coordToInd(int* coord, int &ind){
    /*
     Find the index of site with coord
     */
    for (int i = 0; i < N_enlg; i++){
        if (coord[0]==enlgLattice_[i].coord[0] and coord[1] == enlgLattice_[i].coord[1]){
            ind = enlgLattice_[i].id;
            return true;
        }
    };
    return false;
//    std::cerr<<"Coordinates not found in the enlargerd lattice!"<<std::endl;
//    exit(EXIT_FAILURE);
    };

void Geometry::genBonds(int linkNum, BasisLat *links, std::string bondName){
    if (bondMapNum > MAX_BOND_NUM){
        std::cerr<<"number of bondMap exceeds MAX_BOND_NUM!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    bondNames[bondMapNum] = bondName;
    bondMaps_[bondMapNum] = new BondMap(N);
    int ind;
    BasisLat coord;
    for (int j = 0; j < linkNum; j++){
        for (int i = 0; i < N; i++){
            coord[0] = lattice_[i].coord[0] + links[j][0];
            coord[1] = lattice_[i].coord[1] + links[j][1];
            if (coordToInd(coord, ind)){
                (bondMaps_[bondMapNum]->at(i)).push_back(ind);
            }
        }
    }
    bondMapNum++;
}

void Geometry::genTransList(){
    BasisLat coord;
    TransList.clear();
    std::vector<int> tmp;
    for (int r = 0; r < N; r++){
        tmp.clear();
        for(int i = 0; i < N; i++){
            int j;
            coord[0] = lattice_[i].coord[0] + lattice_[r].coord[0];
            coord[1] = lattice_[i].coord[1] + lattice_[r].coord[1];
            if (coordToInd(coord, j)){
                tmp.push_back(j);
            }else{
                std::cout<<"translation position not found for i = "<<i<<", r = "<<r<<std::endl;
                exit(1);
            }
        }
        TransList.push_back(tmp);
    }
}

void Geometry::print(){
    // print Lattice_
    std::cout<<"Lattice:"<<std::endl;
    for (int i = 0; i < N; i ++){
        std::cout<<"siteID:"<<lattice_[i].id<<", coord:["<<lattice_[i].coord[0]<<", "<<lattice_[i].coord[1]<<"], "<<", coordxy:["<<lattice_[i].coordxy[0]<<", "<<lattice_[i].coordxy[1]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
    // print kLattice_
    std::cout<<"k-Lattice:"<<std::endl;
    for (int i = 0; i < N; i ++){
        std::cout<<"siteID:"<<KLattice_[i].id<<", coord:["<<KLattice_[i].coord[0]<<", "<<KLattice_[i].coord[1]<<"], "<<", coordxy:["<<KLattice_[i].coordxy[0]<<", "<<KLattice_[i].coordxy[1]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
    // print J1 bond
    int bondid = 0;
    std::cout<<"J1 Bonds num:"<<bondMaps_[bondid]->size()<<std::endl;
    for (int i = 0; i < bondMaps_[bondid]->size(); i++){
        std::cout<<"siteID:"<<lattice_[i].id<<", links to:";
        for (int j = 0; j < bondMaps_[bondid]->at(i).size(); j++){
            std::cout<<bondMaps_[bondid]->at(i).at(j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    // print J2 bond
    bondid = 1;
    std::cout<<"J2 Bonds num:"<<bondMaps_[bondid]->size()<<std::endl;
    for (int i = 0; i < bondMaps_[bondid]->size(); i++){
        std::cout<<"siteID:"<<lattice_[i].id<<", links to:";
        for (int j = 0; j < bondMaps_[bondid]->at(i).size(); j++){
            std::cout<<bondMaps_[bondid]->at(i).at(j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    // print TransList
    int rID = 1;
    for (int i = 0; i < N; i++){
        std::cout<<i<<"->"<<TransList.at(rID).at(i)<<std::endl;
    }
}

/*
    ********************************************
    * Triangular Lattice derived from Geometry *
    ********************************************
*/

//constructor with T36 and C6
TriAngLattice::TriAngLattice(int numSites){
    is_PBC = true;
    N = numSites;
    if (is_PBC) N_enlg = N * 7;
    name = "TriAng_C6_N"+std::to_string(N);
    bondMapNum = 0;
    a1[0] = 1; a1[1] = 0;
    a2[0] = 0; a2[1] = 1;
    ax[0] = 1.0; ax[1] = 0.0;
    ay[0] = 0.5; ay[1] = std::sqrt(3.0)/2.0;
    b1[0] = 1.0; b1[1] = 0.0;
    b2[0] = 0.0; b2[1] = 1.0;
    int vecSize = 2;
    
    switch(numSites){
        case 12:
            break;
        case 36:
            kx[0] = 2*PI/6;
            kx[1] = -2*PI/std::sqrt(3.0)/6;
            ky[0] = 0.0;
            ky[1] = 4*PI/std::sqrt(3.0)/6;
            
            lattice_ = new(std::nothrow) Site[N]; assert(lattice_!=NULL);
            enlgLattice_ = new(std::nothrow) Site[N_enlg]; assert(enlgLattice_!=NULL);
            KLattice_ = new(std::nothrow) kSite[N]; assert(KLattice_!=NULL);
            
          
            BasisLat TranVecs[7];
            vecXAdd(0, a1, 0, a2, TranVecs[0], vecSize);
            
            vecXAdd(6, a1, 0, a2, TranVecs[1], vecSize);
            vecXAdd(0, a1, 6 , a2, TranVecs[2], vecSize);
            vecXAdd(-6, a1, 6 , a2, TranVecs[3], vecSize);

            vecXAdd(-6, a1, 0, a2, TranVecs[4], vecSize);
            vecXAdd(0, a1, -6, a2, TranVecs[5], vecSize);
            vecXAdd(6, a1, -6, a2, TranVecs[6], vecSize);
            
            BasisLat rowInitCoord[7];
            int rowCount[7] = {3, 6 ,6 ,6, 6, 6, 3};
            vecXAdd(0, a1, -3, a2, rowInitCoord[0], vecSize);
            vecXAdd(-2, a1, -2, a2, rowInitCoord[1], vecSize);
            vecXAdd(-2, a1, -1, a2, rowInitCoord[2], vecSize);
            vecXAdd(-3, a1, 0, a2, rowInitCoord[3], vecSize);
            vecXAdd(-3, a1, 1, a2, rowInitCoord[4], vecSize);
            vecXAdd(-4, a1, 2, a2, rowInitCoord[5], vecSize);
            vecXAdd(-3, a1, 3, a2, rowInitCoord[6], vecSize);
            // construct lattice
            int counter = 0;
            for (int row = 0; row < 7; row++){
                BasisXY initR;
                vecXAdd(rowInitCoord[row][0], ax, rowInitCoord[row][1], ay, initR, vecSize);
                for (int col = 0; col < rowCount[row]; col++){
                    lattice_[counter].id = counter;
                    lattice_[counter].orbital = 0;
                    vecXAdd(col, a1, 1, rowInitCoord[row], lattice_[counter].coord, vecSize);
                    vecXAdd((double)col, ax, 1.0, initR, lattice_[counter].coordxy, vecSize);
                    counter++;
                }
            }
            // cherck counter==numSites
            assert(counter==N);
            // construct enlarged lattice
            counter = 0;
            for (int k = 0; k < 7; k++){
                for (int i = 0; i < N; i++){
                    enlgLattice_[counter].id = counter % N;
                    vecXAdd(1, lattice_[i].coord, 1, TranVecs[k], enlgLattice_[counter].coord, vecSize);
                    counter++;
                }
            }
            assert(counter==7*N);
            
            // construct k-lattice
            BasisXY bInitCoord[8], kInitCoord[8];
            int krowCount[8] = {2, 5, 6, 6, 6, 6, 4, 1};
            // row 0
            vecXAdd(-2.0/6, b1, -3.0/6, b2, bInitCoord[0], vecSize);
            vecXAdd(-2.0, kx, -3.0, ky, kInitCoord[0], vecSize);
            // row 1
            vecXAdd(-3.0/6, b1, -2.0/6, b2, bInitCoord[1], vecSize);
            vecXAdd(-3.0, kx, -2.0, ky, kInitCoord[1], vecSize);
            // row 2
            vecXAdd(-3.0/6, b1, -1.0/6, b2, bInitCoord[2], vecSize);
            vecXAdd(-3.0, kx, -1.0, ky, kInitCoord[2], vecSize);
            // row 3
            vecXAdd(-2.0/6, b1, 0.0/6, b2, bInitCoord[3], vecSize);
            vecXAdd(-2.0, kx, 0.0, ky, kInitCoord[3], vecSize);
            // row 4
            vecXAdd(-2.0/6, b1, 1.0/6, b2, bInitCoord[4], vecSize);
            vecXAdd(-2.0, kx, 1.0, ky, kInitCoord[4], vecSize);
            // row 5
            vecXAdd(-1.0/6, b1, 2.0/6, b2, bInitCoord[5], vecSize);
            vecXAdd(-1.0, kx, 2.0, ky, kInitCoord[5], vecSize);
            // row 6
            vecXAdd(0.0/6, b1, 3.0/6, b2, bInitCoord[6], vecSize);
            vecXAdd(0.0, kx, 3.0, ky, kInitCoord[6], vecSize);
            // row 7
            vecXAdd(2.0/6, b1, 4.0/6, b2, bInitCoord[7], vecSize);
            vecXAdd(2.0, kx, 4.0, ky, kInitCoord[7], vecSize);
            
            counter = 0;
            for (int row = 0; row < 8; row++){
                for (int col = 0; col < krowCount[row]; col++){
                    KLattice_[counter].id = counter;
                    KLattice_[counter].orbital = 0;
                    vecXAdd((double)col/6, b1, 1.0, bInitCoord[row], KLattice_[counter].coord, vecSize);
                    vecXAdd((double)col, kx, 1.0, kInitCoord[row], KLattice_[counter].coordxy, vecSize);
                    counter++;
                }
            }
            break;
        default: break;
    }
    // construct J1*Si*Sj and J2*Si*Sj bonds
    int  J1Num = 3;
    BasisLat J1Links[J1Num];
    vecXAdd(1, a1, 0, a2, J1Links[0], vecSize);
    vecXAdd(0, a1, 1, a2, J1Links[1], vecSize);
    vecXAdd(1, a1, -1, a2, J1Links[2], vecSize);

    int  J2Num = 3;
    BasisLat J2Links[J1Num];
    vecXAdd(1, a1, 1, a2, J2Links[0], vecSize);
    vecXAdd(-1, a1, 2, a2, J2Links[1], vecSize);
    vecXAdd(2, a1, -1, a2, J2Links[2], vecSize);

    genBonds(J1Num, J1Links, "J1");
    genBonds(J2Num, J2Links, "J2");
    
    // construct translation map list
    genTransList();
};

//constructor with T without C6
TriAngLattice::TriAngLattice(int N1, int N2){
    is_PBC = true;
    Nx = N1;
    Ny = N2;
    N = N1 * N2;
    if (is_PBC) N_enlg = N * 9;
    name = "TriAngNxNy";
    bondMapNum = 0;
    a1[0] = 1;
    a1[1] = 0;
    a2[0] = 0;
    a2[1] = 1;
    ax[0] = 1.0;
    ax[1] = 0.0;
    ay[0] = 0.5;
    ay[1] = std::sqrt(3.0)/2.0;
    b1[0] = 1.0;
    b1[1] = 0.0;
    b2[0] = 0.0;
    b2[1] = 1.0;
    kx[0] = 2*PI/N1;
    kx[1] = -2*PI/std::sqrt(3.0)/N1;
    ky[0] = 0.0;
    ky[1] = 4*PI/std::sqrt(3.0)/N2;
    lattice_ = new(std::nothrow) Site[N]; assert(lattice_!=NULL);
    enlgLattice_ = new(std::nothrow) Site[N_enlg]; assert(enlgLattice_!=NULL);
    KLattice_ = new(std::nothrow) kSite[N]; assert(KLattice_!=NULL);
    
    int vecSize = 2;
    BasisLat TranVecs[9];
    vecXAdd(0, a1, 0, a2, TranVecs[0], vecSize);
    
    vecXAdd(N1, a1, 0, a2, TranVecs[1], vecSize);
    vecXAdd(N1, a1, N2, a2, TranVecs[2], vecSize);
    vecXAdd(0, a1, N2 , a2, TranVecs[3], vecSize);
    vecXAdd(-N1, a1, N2 , a2, TranVecs[4], vecSize);

    vecXAdd(-N1, a1, 0, a2, TranVecs[5], vecSize);
    vecXAdd(-N1, a1, -N2, a2, TranVecs[6], vecSize);
    vecXAdd(0, a1, -N2, a2, TranVecs[7], vecSize);
    vecXAdd(N1, a1, -N2, a2, TranVecs[8], vecSize);
    
    // construct lattice
    int counter = 0;
    for (int j = 0; j < N2; j++){
        for (int i = 0; i < N1; i++ ){
            lattice_[counter].id = counter;
            lattice_[counter].orbital = 0;
            vecXAdd(i, a1, j, a2, lattice_[counter].coord, vecSize);
            vecXAdd((double)i, ax, (double)j, ay, lattice_[counter].coordxy, vecSize);
            counter++;
        }
    }
    
    // construct enlarged lattice
    counter = 0;
    for(int k = 0; k < 9; k++){
        for (int j = 0; j < N2; j++){
            for (int i = 0; i < N1; i++ ){
                enlgLattice_[counter].id = counter % N;
                vecXAdd(i, a1, j, a2, enlgLattice_[counter].coord, vecSize);
                vecXAdd(1, TranVecs[k], 1, enlgLattice_[counter].coord, enlgLattice_[counter].coord, vecSize);
                counter++;
            }
        }
    }
    
    // construc kspace-lattice dual to lattice
    counter = 0;
    for (int j = 0; j < N2; j++){
        for (int i = 0; i < N1; i++ ){
            KLattice_[counter].id = counter;
            KLattice_[counter].orbital = 0;
            vecXAdd((double)i/N1, b1, (double)j/N2, b2, KLattice_[counter].coord, vecSize);
            vecXAdd((double)i, kx, (double)j, ky, KLattice_[counter].coordxy, vecSize);
            counter++;
        }
    }
    
    // construct J1*Si*Sj and J2*Si*Sj bonds
    int  J1Num = 3;
    BasisLat J1Links[J1Num];
    vecXAdd(1, a1, 0, a2, J1Links[0], vecSize);
    vecXAdd(0, a1, 1, a2, J1Links[1], vecSize);
    vecXAdd(1, a1, -1, a2, J1Links[2], vecSize);

    int  J2Num = 3;
    BasisLat J2Links[J1Num];
    vecXAdd(1, a1, 1, a2, J2Links[0], vecSize);
    vecXAdd(-1, a1, 2, a2, J2Links[1], vecSize);
    vecXAdd(2, a1, -1, a2, J2Links[2], vecSize);

    genBonds(J1Num, J1Links, "J1");
    genBonds(J2Num, J2Links, "J2");
    
    // construct translation map list
    genTransList();
};

TriAngLattice::~TriAngLattice(){
    delete [] lattice_;
    delete [] enlgLattice_;
    for (int i = 0; i < bondMapNum; i++){
        delete (bondMaps_[i]);
    }
}

/*
    ***************
    * Basis Class *
    ***************
*/
// initialize
Basis::Basis(int sitedim, int* dimList, Geometry *pt_lat){
    /*
     site is the Hilbert space dimension of a single site.
     dimList is of length sitedim, containing the number of sites in the state of corresponding site Hilbert Space dimension.
    */
    kIndex = -1;
    pt_lattice = pt_lat;
    N = 0;
    siteDim = sitedim;
    kIndex = -1;
    // initialize vec
    for (ind_int i = 0; i < siteDim; i++){
        N += dimList[i];
        for (int j = 0; j < dimList[i]; j++) vec.push_back(i);
    }
    initVec_ = new int[N];
    finalVec_ = new int[N];
    finalInd = new ind_int[N];
    // initialize mul
    for (int i = 0; i < N; i++) mul.push_back(pow((ind_int) siteDim, (ind_int) (N - i - 1)));
    // calculate totDim
    totDim = 1;
    ind_int n = N;
    for (int i = 0; i < siteDim; i++){
        totDim *= combination<ind_int>(n, dimList[i]);
        n -= dimList[i];
    }
}
// generate Basis for the full Hilbert Space
void Basis::gen(){
    kIndex = -1;
    indexList.clear();
    indexList.reserve(totDim);
    ind_int counter = 0;
    do{
        indexList.push_back(vecToInd(vec.data()));
        counter++;
    }while (std::next_permutation(vec.begin(), vec.end()));
    assert(counter == totDim);
}

// generate Basis for the subspace labeled by kInd
void Basis::gen(int kInd){
    if (kInd==-1) gen();
    else{
        kIndex = kInd;
        subDim = 0;
        indexList.clear();
        double norm;
        normList.clear();
        do{
            ind_int initInd = vecToInd(vec.data());
            if (isKRep(kIndex, initInd, pt_lattice, norm)){
                indexList.push_back(initInd);
                subDim++;
            #ifdef KEEP_BASIS_NORM
                normList.push_back(norm);
            #endif
            }
        }while (std::next_permutation(vec.begin(), vec.end()));
    }
}

// construct subspace basis from reps loaded from file
void Basis::gen(int kInd, std::string basisfile){
    kIndex = kInd;
    subDim = 0;
    indexList.clear();
    read<ind_int>(&indexList, basisfile);
    subDim = indexList.size();
}

// construct subspace basis and norm from reps loaded from file
void Basis::gen(int kInd, std::string basisfile, std::string normfile){
    kIndex = kInd;
    subDim = 0;
    indexList.clear();
    read<ind_int>(&indexList, basisfile);
    subDim = indexList.size();
    
    normList.clear();
    read<double>(&normList, normfile);
    assert(normList.size()==indexList.size());
}

// generate Reps for translation symm
void Basis::genRep(){
    kIndex = -1;
    subDim = 0;
    repIndList.clear();
    do{
        ind_int initInd = vecToInd(vec.data());
        indToVec(initInd, initVec_);
        // smallest index in the cycle?
        bool is_min = true;
        for (int r = 0; r < pt_lattice->N; r++){
            for(int i = 0; i < pt_lattice->N; i++){
                finalVec_[pt_lattice->TransList.at(r).at(i)] = initVec_[i];
            }
            if (initInd > vecToInd(finalVec_)) {
                is_min = false;
                break;
            }
        }
        if (is_min) repIndList.push_back(initInd);
    }while (std::next_permutation(vec.begin(), vec.end()));
}

void Basis::saveBasis(std::string basisfile){
    std::ofstream outfile;
    ind_int *d_pt = indexList.data();
    ind_int size = indexList.size();
    save<ind_int>(d_pt, size, &outfile, basisfile);
}

void Basis::saveBasis(std::string basisfile, std::string normfile){
    std::ofstream outfile;
    save<ind_int>(indexList.data(), (ind_int)indexList.size(), &outfile, basisfile);
    save<double>(normList.data(), (ind_int)normList.size(), &outfile, normfile);
}


Basis::~Basis(){
    delete [] initVec_;
    delete [] finalVec_;
    delete [] finalInd;
}

int Basis::getSiteDim(){return siteDim;}

ind_int Basis::getTotDim(){return totDim;}

ind_int Basis::vecToInd(int *v){
    ind_int result = 0;
    for (int i = 0; i < N; i++){if (v[i] != 0) result += v[i] * mul[i];}
    return result;
}

void Basis::indToVec(ind_int index, int* output){
    for (int i = N - 1; i >= 0; i--){
        if (index > 0){
            output[i] = index % siteDim;
            index /= siteDim;
        }
        else{
            output[i] = 0;
        }
    }
}
bool Basis::search(ind_int index, ind_int &ind){
    // Binary search
    ind_int begin = 0, middle, end = subDim;
    while (end >= begin){
        middle = (end + begin) / 2;
        // in case index>indexList[subdim-1], make sure middle<=subdim-1
        if (middle>=subDim) break;
        if (index > indexList.at(middle)){
            begin = middle + 1;
            continue;
        }
        else if (index < indexList.at(middle)){
            // in case we are using unsigned integer and index<indexList[0], make sure middle-1 is not negative
            if (middle>=1){
                end = middle - 1;
                continue;
            }else{
                break;
            }
        }
        else if (index == indexList.at(middle)){
            ind = middle;
            return true;
        }
    }
    if (kIndex==-1){
        std::cerr<<"index not found in indexList (full hilbert space kIndex=-1)!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    return false;
}

ind_int Basis::search(ind_int index){
    // Binary search
    ind_int begin = 0, middle, end = totDim;
    while (end >= begin){
        middle = (end + begin) / 2;
        if (middle>=subDim) break;
        if (index > indexList.at(middle)){
            begin = middle + 1;
            continue;
        }
        else if (index < indexList.at(middle)){
            if (middle>=1){
                end = middle -1;
                continue;
            }else{
                break;
            }
        }
        else if (index == indexList.at(middle)){
            return middle;
        }
    }
    std::cerr<<"index not found in indexList_!"<<std::endl;
    exit(EXIT_FAILURE);
}


/*
    *******************************************
    * Implementation for translation symmetry *
    *******************************************
*/
bool Basis::isKRep(int kInd, ind_int initInd, Geometry* pt_lattice, double& norm){
    // full hilbert space
    if (kInd==-1) return true;
    indToVec(initInd, initVec_);
    BasisLat coord;
    // smallest index in the cycle?
    for (int r = 0; r < pt_lattice->N; r++){
        for(int i = 0; i < pt_lattice->N; i++){
            finalVec_[pt_lattice->TransList.at(r).at(i)] = initVec_[i];
        }
        if (initInd > vecToInd(finalVec_)) return false;
    }
    // norm > 0?
    norm = kNorm(kInd, initInd, pt_lattice);
    // sqrt(infinitismal)>>infinitesimal
    if (std::real(norm)*std::real(norm)>INFINITESIMAL){
        return true;
    }
    return false;
}
double Basis::kNorm(int kInd, ind_int initInd, Geometry* pt_lattice){
    double kx = pt_lattice->KLattice_[kIndex].coord[0];
    double ky = pt_lattice->KLattice_[kIndex].coord[1];
    indToVec(initInd, initVec_);
    BasisLat coord;
    dataType norm = 0.0;
    for (int r = 0; r < pt_lattice->N; r++){
        int rx = pt_lattice->lattice_[r].coord[0];
        int ry = pt_lattice->lattice_[r].coord[1];
        for(int i = 0; i < pt_lattice->N; i++){
            finalVec_[pt_lattice->TransList.at(r).at(i)] = initVec_[i];
        }
        if (vecToInd(finalVec_)==initInd) norm += std::exp(2*PI*CPLX_I*(kx*rx+ky*ry));
    }
    norm /= pt_lattice->N;
    assert(std::abs(std::imag(norm))<INFINITESIMAL);
    return std::sqrt(std::real(norm));
}
void Basis::genTranslation(Geometry* pt_lattice, ind_int ind, ind_int* finalInd){
    indToVec(ind, initVec_);
    BasisLat coord;
    for (int r = 0; r < pt_lattice->N; r++){
        for(int i = 0; i < pt_lattice->N; i++){
            finalVec_[pt_lattice->TransList.at(r).at(i)] = initVec_[i];
        }
        finalInd[r] = vecToInd(finalVec_);
    }
}

/*
    ************************
    * Spin Operators Class *
    ************************
*/
SpinOperator::SpinOperator(){
    spinDim = SPIN_DIM;
    double s = (double)(spinDim - 1)/2.0;
    double m = s;
    for (int i = 0; i < spinDim; i++){
        szMat.push_back(m);
        spMat.push_back(std::sqrt(s*(s+1.0)-m*(m+1.0)));
        smMat.push_back(std::sqrt(s*(s+1.0)-m*(m-1.0)));
        m -= 1.0;
    }
}

SpinOperator::SpinOperator(int dim){
    spinDim = dim;
    double s = (double)(spinDim - 1)/2.0;
    double m = s;
    for (int i = 0; i < spinDim; i++){
        szMat.push_back(m);
        spMat.push_back(std::sqrt(s*(s+1.0)-m*(m+1.0)));
        smMat.push_back(std::sqrt(s*(s+1.0)-m*(m-1.0)));
        m -= 1.0;
    }
}

SpinOperator::~SpinOperator(){};

/*
    *********************
    * Hamiltonian Class *
    *********************
*/

// For spin-1/2 Heisenburg
HeisenbergOneHalf::HeisenbergOneHalf(ind_int totDim, dataType J1, dataType J2){
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    parameters.push_back(J1); // J1
    parameters.push_back(J2); // J2
    dim = totDim;
    nlocmax = (dim + workerNum - 1)/workerNum;
    ntot = nlocmax * workerNum;
    startRow = workerID * nlocmax;
    endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
    nloc = endRow - startRow;
#ifdef SAXPY
    matrix_num = 2;
    assert(matrix_num<=MAX_MATRIX_NUM);
    pt_valList[0] = &valListJ1;
    pt_valList[1] = &valListJ2;
    pt_colList[0] = &colListJ1;
    pt_colList[1] = &colListJ2;
    pt_rowColInitList[0] = &rowColInitListJ1;
    pt_rowColInitList[1] = &rowColInitListJ2;
#endif
}

HeisenbergOneHalf::~HeisenbergOneHalf(){};

void HeisenbergOneHalf::setJ2(dataType J2){
    parameters.at(1) = J2;
}

#ifdef SAXPY
// generate Hamiltonian in the subspacd labeled by kIndex
 void HeisenbergOneHalf::genSubMatMap(int kIndex, Geometry *pt_lattice, Basis *pt_Basis, int couplingNum, int polarNum){
     assert(couplingNum==matrix_num);
     if (kIndex==-1){
         genMat(pt_lattice, pt_Basis, couplingNum, polarNum);
         return;
     }
     clear();
     for (int matID = 0; matID < matrix_num; matID++){
         pt_rowColInitList[matID]->reserve(nloc+1);
         pt_valList[matID]->reserve(nloc*polarNum*pt_lattice->N/2+nloc);
         pt_colList[matID]->reserve(nloc*polarNum*pt_lattice->N/2+nloc);
     }
     MAP rowMap;
     MAPIT it;

     double kx = pt_lattice->KLattice_[kIndex].coord[0];
     double ky = pt_lattice->KLattice_[kIndex].coord[1];
     ind_int *finalIndList_ = new(std::nothrow) ind_int[pt_lattice->N]; assert(finalIndList_!=NULL);
     int *initVec_ = new(std::nothrow) int[pt_lattice->N]; assert(initVec_!=NULL);
     ind_int counter;
     double initNorm, finalNorm;
     
     // calculate <R1k|H*Pk|R2k>/norm1/norm2
     for (int couplingID = 0; couplingID < couplingNum; couplingID++){
         counter = 0;
         pt_rowColInitList[couplingID]->push_back(counter);
         for (ind_int rowID = startRow; rowID < endRow; rowID++){
             // clean map
             rowMap.clear();
             pt_Basis->genTranslation(pt_lattice, pt_Basis->indexList.at(rowID), finalIndList_);
             initNorm = pt_Basis->getNorm(rowID);
             for (int i = 0; i < pt_lattice->N; i++){
                 if (DEBUG) std::cout<<"trans vec:"<<i<<" begins..."<<std::endl;
                 pt_Basis->indToVec(finalIndList_[i], initVec_);
                 int rx = pt_lattice->lattice_[i].coord[0];
                 int ry = pt_lattice->lattice_[i].coord[1];
                 cdouble factor = std::exp(2*PI*CPLX_I*(kx*rx+ky*ry))/pt_lattice->N/initNorm;
                 for (int siteID = 0; siteID < pt_lattice->N; siteID++){
                     for (int polarID = 0; polarID < polarNum; polarID++){
                         int siteIDp = pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID);
                         // sz.siteID * sz.siteIDp
                         szsz(siteID, siteIDp, parameters.at(couplingID)*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
                         // 1/2 * sm.siteID * sp.siteIDp
                         spsm(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
                         // 1/2 * sp.siteID * sm.siteIDp
                         smsp(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
                     }
                 }
             }
             newRow(couplingID, &rowMap, counter);
         }
     }
     delete [] initVec_;
     delete [] finalIndList_;
 }

void HeisenbergOneHalf::genMat(Geometry *pt_lattice, Basis *pt_Basis, int couplingNum, int polarNum){
    clear();
    for (int matID = 0; matID < matrix_num; matID++){
         pt_rowColInitList[matID]->reserve(nloc+1);
         pt_valList[matID]->reserve(nloc*polarNum*pt_lattice->N/2+nloc);
         pt_colList[matID]->reserve(nloc*polarNum*pt_lattice->N/2+nloc);
     }
    MAP rowMap;
    MAPIT it;

    ind_int counter;
    int *initVec_ = new int[pt_lattice->N];
    for (int couplingID = 0; couplingID < couplingNum; couplingID++){
        counter = 0;
        pt_rowColInitList[couplingID]->push_back(counter);
        for (ind_int rowID = startRow; rowID < endRow; rowID++){
            rowMap.clear();
            pt_Basis->indToVec(pt_Basis->indexList.at(rowID), initVec_);
            for (int siteID = 0; siteID < pt_lattice->N; siteID++){
                for (int polarID = 0; polarID < polarNum; polarID++){
                    int siteIDp = pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID);
                    // sz.siteID * sz.siteIDp
                    szsz(siteID, siteIDp, parameters.at(couplingID), rowID, initVec_, pt_Basis, &rowMap);
                    // 1/2 * sm.siteID * sp.siteIDp
                    spsm(siteID, siteIDp, parameters.at(couplingID)/2.0, rowID, initVec_, pt_Basis, &rowMap);
                    // 1/2 * sp.siteID * sm.siteIDp
                    smsp(siteID, siteIDp, parameters.at(couplingID)/2.0, rowID, initVec_, pt_Basis, &rowMap);
                }
            }
            newRow(couplingID, &rowMap, counter);
        }
    }
    delete [] initVec_;
}
#else
// generate Hamiltonian in the subspacd labeled by kIndex
    void HeisenbergOneHalf::genSubMatMap(int kIndex, Geometry *pt_lattice, Basis *pt_Basis, int couplingNum, int polarNum){
    if (kIndex==-1){
        genMat(pt_lattice, pt_Basis, couplingNum, polarNum);
        return;
    }
        clear();
        rowColInitList.reserve(nloc+1);
        colList.reserve(nloc*couplingNum*polarNum*pt_lattice->N/2+nloc);
        valList.reserve(nloc*couplingNum*polarNum*pt_lattice->N/2+nloc);
        MAP rowMap;
        MAPIT it;
        
        double kx = pt_lattice->KLattice_[kIndex].coord[0];
        double ky = pt_lattice->KLattice_[kIndex].coord[1];
        ind_int *finalIndList_ = new(std::nothrow) ind_int[pt_lattice->N]; assert(finalIndList_!=NULL);
        int *initVec_ = new(std::nothrow) int[pt_lattice->N]; assert(initVec_!=NULL);
        ind_int counter = 0;
        double initNorm, finalNorm;
        rowColInitList.push_back(counter);
        
        // calculate <R1k|H*Pk|R2k>/norm1/norm2
        for (ind_int rowID = startRow; rowID < endRow; rowID++){
            rowMap.clear();
            if (DEBUG) std::cout<<"row:"<<rowID<<" begins..."<<std::endl;
            pt_Basis->genTranslation(pt_lattice, pt_Basis->indexList.at(rowID), finalIndList_);
            initNorm = pt_Basis->getNorm(rowID);
            for (int i = 0; i < pt_lattice->N; i++){
                pt_Basis->indToVec(finalIndList_[i], initVec_);
                int rx = pt_lattice->lattice_[i].coord[0];
                int ry = pt_lattice->lattice_[i].coord[1];
                cdouble factor = std::exp(2*PI*CPLX_I*(kx*rx+ky*ry))/pt_lattice->N/initNorm;
                for (int siteID = 0; siteID < pt_lattice->N; siteID++){
                    for (int couplingID = 0; couplingID < couplingNum; couplingID++){
                        for (int polarID = 0; polarID < polarNum; polarID++){
                            int siteIDp = pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID);
                            // sz.siteID * sz.siteIDp
                            szsz(siteID, siteIDp, parameters.at(couplingID)*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
                            // 1/2 * sm.siteID * sp.siteIDp
                            spsm(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
                            // 1/2 * sp.siteID * sm.siteIDp
                            smsp(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
                        }
                    }
                }
            }
            newRow(&rowMap, counter);
        }
        delete [] initVec_;
        delete [] finalIndList_;
    }

    void HeisenbergOneHalf::genMat(Geometry *pt_lattice, Basis *pt_Basis, int couplingNum, int polarNum){
        clear();
        rowColInitList.reserve(nloc+1);
        colList.reserve(nloc*couplingNum*polarNum*pt_lattice->N/2+nloc);
        valList.reserve(nloc*couplingNum*polarNum*pt_lattice->N/2+nloc);
        MAP rowMap;
        MAPIT it;

        ind_int counter = 0;
        dataType dval;
        rowColInitList.push_back(counter);
        int *initVec_ = new int[pt_lattice->N];
        for (ind_int rowID = startRow; rowID < endRow; rowID++){
            rowMap.clear();
            pt_Basis->indToVec(pt_Basis->indexList.at(rowID), initVec_);
            for (int siteID = 0; siteID < pt_lattice->N; siteID++){
                for (int couplingID = 0; couplingID < couplingNum; couplingID++){
                    for (int polarID = 0; polarID < polarNum; polarID++){
                        int siteIDp = pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID);
                        // sz.siteID * sz.siteIDp
                        szsz(siteID, siteIDp, parameters.at(couplingID), rowID, initVec_, pt_Basis, &rowMap);
                        // 1/2 * sm.siteID * sp.siteIDp
                        spsm(siteID, siteIDp, parameters.at(couplingID)/2.0, rowID, initVec_, pt_Basis, &rowMap);
                        // 1/2 * sp.siteID * sm.siteIDp
                        smsp(siteID, siteIDp, parameters.at(couplingID)/2.0, rowID, initVec_, pt_Basis, &rowMap);
                    }
                }
            }
            newRow(&rowMap, counter);
        }
        delete [] initVec_;
    }
#endif

//// Off-diagonal part
//void HeisenbergOneHalf::genMatOffDiag(Geometry *pt_lattice, Basis *pt_Basis, int couplingID, int polarNum){
//    ind_int counter = 0;
//    rowColInitList.push_back(counter);
//    int *initVec_ = new int[pt_lattice->N];
//    for (ind_int rowID = startRow; rowID < endRow; rowID++){
//        pt_Basis->indToVec(pt_Basis->indexList_[rowID], initVec_);
//        for (int siteID = 0; siteID < pt_lattice->N; siteID++){
//            for (int polarID = 0; polarID < polarNum; polarID++){
//                if (initVec_[siteID] == 1 && initVec_[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)] == 0){
//                    ind_int finalIndex = pt_Basis->indexList_[rowID] - pt_Basis->mul[siteID] + pt_Basis->mul[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)];
//                    ind_int colID = pt_Basis->search(finalIndex);
//                    valList.push_back(parameters.at(couplingID)/2);
//                    colList.push_back(colID);
//                    counter++;
//
//                }
//                if (initVec_[siteID] == 0 && initVec_[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)] == 1){
//                    ind_int finalIndex = pt_Basis->indexList_[rowID] + pt_Basis->mul[siteID] - pt_Basis->mul[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)];
//                    ind_int colID = pt_Basis->search(finalIndex);
//                    valList.push_back(parameters.at(couplingID)/2);
//                    colList.push_back(colID);
//                    counter++;
//                }
//            }
//        }
//    rowColInitList.push_back(counter);
//    }
//    delete [] initVec_;
//}

//// Off-diagonal part with specific polarization (rj - ri). necessary for time evolve and Raman. conjugate parts are also needed to be stored seperately.
//void HeisenbergOneHalf::genMatPol(Geometry* pt_lattice, Basis *pt_Basis, int couplingID, int polarID){
//    ind_int counter = 0;
//    ind_int counterConj = 0;
//    rowColInitList.push_back(counter);
//    rowColInitListConj.push_back(counterConj);
//    int *initVec_ = new int[pt_lattice->N];
//    for (ind_int rowID = startRow; rowID < endRow; rowID++){
//        pt_Basis->indToVec(pt_Basis->indexList_[rowID], initVec_);
//        for (int siteID = 0; siteID < pt_lattice->N; siteID++){
//            if (initVec_[siteID] == 1 && initVec_[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)] == 0){
//                ind_int finalIndex = pt_Basis->indexList_[rowID] - pt_Basis->mul[siteID] + pt_Basis->mul[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)];
//                ind_int colID = pt_Basis->search(finalIndex);
//                valList.push_back(parameters.at(couplingID)/2);
//                colList.push_back(colID);
//                counter++;
//
//            }
//            if (initVec_[siteID] == 0 && initVec_[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)] == 1){
//                ind_int finalIndex = pt_Basis->indexList_[rowID] + pt_Basis->mul[siteID] - pt_Basis->mul[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)];
//                ind_int colID = pt_Basis->search(finalIndex);
//                valListConj.push_back(parameters.at(couplingID)/2);
//                colListConj.push_back(colID);
//                counterConj++;
//            }
//        }
//        rowColInitList.push_back(counter);
//        rowColInitListConj.push_back(counterConj);
//    }
//    delete[] initVec_;
//}

//// Diagonal Part
//void HeisenbergOneHalf::genDiagMat(Geometry *pt_lattice, Basis *pt_Basis, int couplingID, int polarNum){
//    dataType ddata;
//    int *initVec_ = new int[pt_lattice->N];
//    for (ind_int rowID = startRow; rowID < endRow; rowID++){
//        ddata = 0;
//        pt_Basis->indToVec(pt_Basis->indexList_[rowID], initVec_);
//        for (int siteID = 0; siteID < pt_lattice->N; siteID++){
//            for (int polarID = 0; polarID < polarNum; polarID++){
//                ddata += parameters.at(couplingID) * sz[initVec_[siteID]] * sz[initVec_[pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID)]];
//            }
//        }
//        diagValList.push_back(ddata);
//    }
//    delete [] initVec_;
//}

//template<class T>
//void HeisenbergOneHalf::save(std::vector<T> *pt_vec, std::ofstream *pt_file, std::string filename){
//    pt_file->open(filename,std::ios::binary);
//    if (pt_file->is_open()){
//        pt_file->write(reinterpret_cast<char*>(&(pt_vec->at(0))),pt_vec->size() * sizeof(T));
//        pt_file->close();
//    }else{
//        std::cout<<filename<<" failed to open!"<<std::endl;
//    }
//}
//void HeisenbergOneHalf::saveAll(std::string dataDir){
//    std::ofstream outFile;
//    save<dataType>(&diagValList, &outFile, dataDir + "Hdiag" + std::to_string(workerID)+".dat");
//    save<dataType>(&valList, &outFile, dataDir + "H" + std::to_string(workerID)+".dat");
//    save<ind_int>(&colList, &outFile, dataDir + "colList" + std::to_string(workerID)+".dat");
//    save<ind_int>(&rowColInitList, &outFile, dataDir + "rowColInitList" + std::to_string(workerID)+".dat");
//}

/*
    *******************
    * SzSz Correlator *
    *******************
*/

SzqOneHalf::SzqOneHalf(ind_int totDim){
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    dim = totDim;
    nlocmax = (dim + workerNum - 1)/workerNum;
    ntot = nlocmax * workerNum;
    startRow = workerID * nlocmax;
    endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
    nloc = endRow - startRow;
}

SzqOneHalf::~SzqOneHalf(){}

void SzqOneHalf::genMat(Geometry* pt_lattice, Basis* pt_Basis, BasisXY q){
    clear();
    diagValList.reserve(nloc);
    cdouble dval;
    int *initVec_ = new(std::nothrow) int[pt_lattice->N]; assert(initVec_!=NULL);
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        dval = 0.0;
        pt_Basis->indToVec(pt_Basis->indexList.at(rowID), initVec_);
        for (int siteID = 0; siteID < pt_lattice->N; siteID++){
            dval += szMat.at(initVec_[siteID]) * std::exp(-CPLX_I * (q[0]*pt_lattice->lattice_[siteID].coordxy[0] + q[1]*pt_lattice->lattice_[siteID].coordxy[1]));
        }
        diagValList.push_back(dval);
    }
    delete [] initVec_;
}

// generate
void SzqOneHalf::genMat(Geometry* pt_lattice, Basis* pt_B1, Basis* pt_B2, BasisXY q){
    clear();
    rowColInitList.reserve(nloc+1);
    colList.reserve(nloc);
    valList.reserve(nloc);
    cdouble dval;
    int *initVec_ = new(std::nothrow) int[pt_lattice->N]; assert(initVec_!=NULL);
    switch(PARTITION){
        case ROW_PARTITION:{
            ind_int counter = 0;
            rowColInitList.push_back(counter);
            ind_int colID;
            for (ind_int rowID = startRow; rowID < endRow; rowID++){
                if (pt_B1->search(pt_B2->indexList.at(rowID),colID)){
                    dval = 0.0;
                    pt_B1->indToVec(pt_B1->indexList.at(colID), initVec_);
                    for (int siteID = 0; siteID < pt_lattice->N; siteID++){
                        dval += szMat.at(initVec_[siteID]) * std::exp(-CPLX_I * (q[0]*pt_lattice->lattice_[siteID].coordxy[0] + q[1]*pt_lattice->lattice_[siteID].coordxy[1]));
                    }
                    dval *= pt_B2->getNorm(rowID)/pt_B1->getNorm(colID);
                    colList.push_back(colID);
                    valList.push_back(dval);
                    counter++;
                }
                rowColInitList.push_back(counter);
            }
            break;
        }
        case COL_PARTITION:{
            ind_int counter = 0;
            rowColInitList.push_back(counter);
            ind_int colID;
            for (ind_int rowID = startRow; rowID < endRow; rowID++){
                if (pt_B1->search(pt_B2->indexList.at(rowID),colID)){
                    dval = 0.0;
                    pt_B2->indToVec(pt_B2->indexList.at(rowID), initVec_);
                    for (int siteID = 0; siteID < pt_lattice->N; siteID++){
                        dval += szMat[initVec_[siteID]] * std::exp(-CPLX_I * (q[0]*pt_lattice->lattice_[siteID].coordxy[0] + q[1]*pt_lattice->lattice_[siteID].coordxy[1]));
                    }
                    dval *= pt_B1->getNorm(rowID)/pt_B2->getNorm(colID);
                    colList.push_back(colID);
                    valList.push_back(dval);
                    counter++;
                }
                rowColInitList.push_back(counter);
            }
            break;
        }
    }
    delete [] initVec_;
}

/*
    *****************
    * SS Correlator *
    *****************
*/
SSOneHalf::SSOneHalf(ind_int totDim){
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    dim = totDim;
    nlocmax = (dim + workerNum - 1)/workerNum;
    ntot = nlocmax * workerNum;
    startRow = workerID * nlocmax;
    endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
    nloc = endRow - startRow;
}

SSOneHalf::~SSOneHalf(){}

// generate Si*Sj in full hilbert space
void SSOneHalf::genPairMat(Geometry* pt_lattice, Basis* pt_Basis, int siteID, int initSiteID){
    clear();
    rowColInitList.reserve(nloc+1);
    colList.reserve(nloc*3);
    valList.reserve(nloc*3);
    MAP rowMap;
    MAPIT it;

    ind_int counter = 0;
    rowColInitList.push_back(counter);
    int *initVec_ = new int[pt_lattice->N];
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        rowMap.clear();
        pt_Basis->indToVec(pt_Basis->indexList.at(rowID), initVec_);
        // sz.siteID * sz.siteIDp
        szsz(siteID, initSiteID, 1.0, rowID, initVec_, pt_Basis, &rowMap);
        // 1/2 * sm.siteID * sp.siteIDp
        spsm(siteID, initSiteID, 0.5, rowID, initVec_, pt_Basis, &rowMap);
        // 1/2 * sp.siteID * sm.siteIDp
        smsp(siteID, initSiteID, 0.5, rowID, initVec_, pt_Basis, &rowMap);

        newRow(&rowMap, counter);
    }
    delete [] initVec_;
}

// generate matrix in subsapce labeled by kIndex for sum.r:Sr*Sr+dr, dr is labeled by rIndex
void SSOneHalf::genPairMat(int kIndex, Geometry* pt_lattice, Basis* pt_Basis, int rIndex){
    clear();
    rowColInitList.reserve(nloc+1);
    colList.reserve(nloc*pt_lattice->N/2+nloc);
    valList.reserve(nloc*pt_lattice->N/2+nloc);
    MAP rowMap;
    MAPIT it;

    double kx = pt_lattice->KLattice_[kIndex].coord[0];
    double ky = pt_lattice->KLattice_[kIndex].coord[1];
    int drx = pt_lattice->lattice_[rIndex].coord[0];
    int dry = pt_lattice->lattice_[rIndex].coord[1];
    ind_int *finalIndList_ = new(std::nothrow) ind_int[pt_lattice->N]; assert(finalIndList_!=NULL);
    int *initVec_ = new(std::nothrow) int[pt_lattice->N]; assert(initVec_!=NULL);

    ind_int counter = 0;
    cdouble dval;
    double initNorm, finalNorm;

    rowColInitList.push_back(counter);
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        rowMap.clear();
        initNorm = pt_Basis->getNorm(rowID);
        pt_Basis->genTranslation(pt_lattice, pt_Basis->indexList.at(rowID), finalIndList_);
        for (int i = 0; i < pt_lattice->N; i++){
            pt_Basis->indToVec(finalIndList_[i], initVec_);
            dval = 0.0;
            int rx = pt_lattice->lattice_[i].coord[0];
            int ry = pt_lattice->lattice_[i].coord[1];
            cdouble factor = std::exp(2*PI*CPLX_I*(kx*rx+ky*ry))/pt_lattice->N/initNorm;
            for (int siteID = 0; siteID < pt_lattice->N; siteID++){
                int siteIDp;
                BasisLat coord;
                coord[0] = pt_lattice->lattice_[siteID].coord[0] + drx;
                coord[1] = pt_lattice->lattice_[siteID].coord[1] + dry;
                assert(pt_lattice->coordToInd(coord, siteIDp));
                // sz.siteID * sz.siteIDp
                szsz(siteID, siteIDp, factor, finalIndList_[i], initVec_, pt_Basis, &rowMap);
                // 1/2 * sm.siteID * sp.siteIDp
                spsm(siteID, siteIDp, factor/2.0, finalIndList_[i], initVec_, pt_Basis, &rowMap);
                // 1/2 * sp.siteID * sm.siteIDp
                smsp(siteID, siteIDp, factor/2.0, finalIndList_[i], initVec_, pt_Basis, &rowMap);
            }
        }
        newRow(&rowMap, counter);
    }
    delete [] initVec_;
    delete [] finalIndList_;
}
