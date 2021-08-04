#ifndef __BITOP_H__
#define __BITOP_H__

#include "Global/globalType.hpp"

/*
    ******************
    * Bit Operations *
    * ****************
*/
inline idx_t bitMask(int pos){
    return idx_t{1}<<pos;
}
inline bool bitTest(idx_t n, int pos){
    return (n>>pos) & idx_t{1};
}
inline void bitSet(idx_t& n, int pos){
    n |= bitMask(pos);
}
inline void bitFlip(idx_t& n, int pos){
    n ^= bitMask(pos);
}
inline idx_t bitCount(idx_t& n, const VecI& idxs){
    idx_t sum = 0;
    for(auto it=idxs.begin(); it!=idxs.end(); it++) sum += n>>(*it) & idx_t{1};
    return sum; 
}
inline int bitCount(idx_t n, int len) {
    int count = 0;
    for (int i = 0; i < len; ++i) {
        if (bitTest(n, i)) {
            ++count;
        }
    }
    return count;
}
inline void bitPrint(idx_t n, int range){
    for(int pos=range-1;pos>=0;pos--)std::cout<<((n>>pos) & idx_t{1});
    std::cout<<std::endl;
}
// calculate (-1)^# of permutations -> fermion sign
inline double seqSign(std::vector<int>& seq){
    int counter = 0;
    for(auto it = seq.begin(); it != seq.end(); it++){
        for (auto itp = it + 1; itp != seq.end(); itp++){
            if (*it > *itp) counter++;
        }
    }
    if (counter%2==0) return 1.0;
    return -1.0;
}
// From Yao's code
template<class UnsignedType>
inline UnsignedType nextLexicographicalNumber(UnsignedType x) {
    if(x==0)return x+1;
    UnsignedType t = (x | (x - 1)) + 1; //find pivot position and set it to 1
    return t | ((((t & -t) / (x & -x)) >> 1) - 1); //reverse bits after the pivot
}

#endif // __BITOP_H__