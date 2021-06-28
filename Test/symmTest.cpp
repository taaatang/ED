#include "Geometry/Geometry.hpp"
using namespace std;
int main() {
	TriAngLattice latt(36);
	latt.addOrb({});
	latt.construct();
	
	cout << "Test Gt * Gt = Gt\n";
	printLine();
	for (int kidx = 0; kidx < latt.getSiteNum(); kidx++) {
		auto T1 = latt.getGT(kidx);
		auto T2 = T1 * T1;
		if (T1 == T2) {
			std::cout << "kidx = " << kidx << ", Gt * Gt = Gt" << std::endl;
		} else {
			std::cout << "kidx = " << kidx << ", Gt * Gt != Gt" << std::endl;	
		}
	}

	cout << "Test Pt * Pt = Pt\n";
	printLine();
	for (int kidx = 0; kidx < latt.getPGRepNum(); kidx++) {
		auto T1 = latt.getGP(kidx);
		auto T2 = T1 * T1;
		if (T1 == T2) {
			std::cout << "pidx = " << kidx << ", Gt * Gt = Gt" << std::endl;
		} else {
			std::cout << "pidx = " << kidx << ", Gt * Gt != Gt" << std::endl;	
		}
		cout << T1 << endl;
	}	
	
	cout << "Test Gt * Gp = Gp * Gt\n";
	printLine();
	for (int kidx = 0; kidx < latt.getSiteNum(); kidx++) {
		printLine();
		auto T = latt.getGT(kidx);
		for (int pidx = 0; pidx < latt.getPGRepNum(); pidx++) {
			auto P = latt.getGP(pidx);
			// auto TP = T * P;
			// auto PT = P * T;
			if (commute(T, P)) {
				std::cout << "kidx = " << kidx << ", pidx = " << pidx << ": commute" << std::endl;
			} else {
				std::cout << "kidx = " << kidx << ", pidx = " << pidx << ": not commute" << std::endl;	
			}
		}	
	}	
}