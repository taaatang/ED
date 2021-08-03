#pragma once

#include <vector>
#include "Geometry/Transform.hpp"

using idx_t = uint64_t;
using sidx_t = int64_t;

class BasisStateInterface {
public:
	virtual bool next() = 0;

	virtual double transform(const Transform<cdouble> &u) = 0;

	virtual idx_t getTotDim() = 0;

public:
    static void configureBase(int nSite_, int nu_, int nd_, int maxPhPerSite_);

    static int getNSite() { return nSite; }

    static bool isConfigured() { return configured; }

    static int getNu() { return nu; }

    static int getNd() { return nd; }

    static int getMaxPhPerSite() { return maxPhPerSite; }

protected:
    static bool configured;

	static int nSite;

	static int nu;

	static int nd;

	static int maxPhPerSite;
};

struct BinaryState {
public:
    BinaryState() = default;

	BinaryState(int nSite, int nParticle);

	void next();

	int count(int i) const { return bitTest(state, i); }

	idx_t& operator()(){ return state;}

public:
	static BinaryState min(int nSite, int nParticle);

	static BinaryState max(int nSite, int nParticle);

	static idx_t getDim(int nSite, int nParticle);

	void print(std::ostream& os, int n) const;

public:
    idx_t state{0};
};

bool operator<(const BinaryState& lhs, const BinaryState& rhs);

bool operator==(const BinaryState& lhs, const BinaryState& rhs);

std::ostream& operator<<(std::ostream& os, BinaryState b);

struct PhononState {
public:
    PhononState() = default;

    PhononState(int nSite, uint8_t maxPhononPerSite);

    void next(int maxPhPerSite);

    int count(int i) const { return state.at(i); };

    std::vector<uint8_t>& operator()(){ return state; };

public:
    static PhononState min(int nSite);

    static PhononState max(int nSite, int nMaxPhPerSite);

    static idx_t getDim(int nSite, int nMaxPhPerSite);

public:
    std::vector<uint8_t> state{};
};

bool operator<(const PhononState& lhs, const PhononState& rhs);

bool operator==(const PhononState& lhs, const PhononState& rhs);

std::ostream& operator<<(std::ostream& os, const PhononState& p);

class SpinBasis : public BasisStateInterface{
public:
    SpinBasis();

    bool next() override;

    double transform(const Transform<cdouble> &u) override;

    static void configure(int nSite_, int nu_, int nd_, int maxPhPerSite_);

    idx_t getTotDim() override;

public:
    BinaryState state;

private:
    static bool isMinMaxSet;

    static BinaryState min;

    static BinaryState max;
};

bool operator<(const SpinBasis& lhs, const SpinBasis& rhs);

bool operator==(const SpinBasis& lhs, const SpinBasis& rhs);

std::ostream& operator<<(std::ostream& os, const SpinBasis& s);

class ElectronBasis : public BasisStateInterface {
public:
    ElectronBasis();

    ElectronBasis(int n_, int nu_, int n_d);

    ElectronBasis(BinaryState up_, BinaryState dn_) : up(up_), dn(dn_) { }

    static void configure(int nSite_, int nu_, int nd_, int maxPhPerSite_);

    bool next() override;

    double transform(const Transform<cdouble> &u) override;

    idx_t getTotDim() override;

    bool isDoubleOcc() const { return  up.state & dn.state; }

public:
	BinaryState up;

	BinaryState	dn;

private:
    static bool isMinMaxSet;

	static BinaryState upMin;

    static BinaryState upMax;

	static BinaryState dnMin;

    static BinaryState dnMax;

    static ElectronBasis min;

    static ElectronBasis max;

    static bool allowDoubleOcc;

public:
    static const BinaryState &getUpMin() {
        return upMin;
    }

    static const BinaryState &getUpMax() {
        return upMax;
    }

    static const BinaryState &getDnMin() {
        return dnMin;
    }

    static const BinaryState &getDnMax() {
        return dnMax;
    }

    static const ElectronBasis &getMin() {
        return min;
    }

    static const ElectronBasis &getMax() {
        return max;
    }

    static bool isAllowDoubleOcc() {
        return allowDoubleOcc;
    }

    static void setAllowDoubleOcc(bool allowDoubleOcc_) {
        ElectronBasis::allowDoubleOcc = allowDoubleOcc_;
    }
};

bool operator<(const ElectronBasis& lhs, const ElectronBasis& rhs);

bool operator==(const ElectronBasis& lhs, const ElectronBasis& rhs);

std::ostream& operator<<(std::ostream& os, const ElectronBasis& e);

class ElectronPhononBasis : public BasisStateInterface {
public:
    ElectronPhononBasis();

    ElectronPhononBasis(int n_, int nu_, int nd_, int nPho_);

    ElectronPhononBasis(ElectronBasis el_, PhononState ph_) : el(el_), ph(std::move(ph_)) { }

    static void configure(int nSite_, int nu_, int nd_, int maxPhPerSite_);

    bool next() override;

    double transform(const Transform<cdouble> &u) override;

    idx_t getTotDim() override;

public:
    ElectronBasis el;

    PhononState ph;

private:
    static bool isMinMaxSet;

    static PhononState phMin;

    static PhononState phMax;

    static ElectronPhononBasis min;

    static ElectronPhononBasis max;
};

bool operator<(const ElectronPhononBasis& lhs, const ElectronPhononBasis& rhs);

bool operator==(const ElectronPhononBasis& lhs, const ElectronPhononBasis& rhs);

std::ostream& operator<<(std::ostream& os, const ElectronPhononBasis& ep);
