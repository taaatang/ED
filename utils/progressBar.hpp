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
	ProgressBar(const std::string& taskName, uint32_t tot, uint32_t segCount, bool allow = true) : task(taskName), total(tot), segmentCount(segCount), segmentLength(tot / segCount), allowPrint(allow) { 
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
					std::cout << '\n' << task << " finished!\n" << std::flush;
					finished = true;
				}
			}
		}
	}

private:
	void init() {
		if (allowPrint) {
			std::cout << "Task Name : " << task << std::endl;
			std::cout << "Total Work: " << std::string(segmentCount, '<') << std::endl;
			std::cout << "Progress  : "<<std::flush;
		}
	}

private:
	std::string task;
	uint32_t total;
	uint32_t segmentCount;
	uint32_t segmentLength;
	bool allowPrint;
	uint32_t count{0};
	bool finished{false};
};