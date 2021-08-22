/**
 * @file progressBar.hpp
 * @author Ta Tang (tatang.physics@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-08-22
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once

#include <iostream>

class ProgressBar {
public:
	ProgressBar(uint32_t tot, uint32_t segCount, bool allow = true) : total(tot), segmentCount(segCount), segmentLength(tot / segCount), allowPrint(allow) { 
		init();
	}

	void progress() {
		if (allowPrint) {
			count++;
			if (count < total) {
				if (count % segmentLength == 0) {
					std::cout << ">" << std::flush;
				}
			} else {
				if (!finished) {
					std::cout << "\nFinished!\n" << std::flush;
					finished = true;
				}
			}
		}
	}

private:
	void init() {
		if (allowPrint) {
			std::cout << "Total Work: " << std::string(segmentCount, '<') << std::endl;
			std::cout << "Progress  : "<<std::flush;
		}
	}

private:
	uint32_t total;
	uint32_t segmentCount;
	uint32_t segmentLength;
	bool allowPrint;
	uint32_t count{0};
	bool finished{false};
};