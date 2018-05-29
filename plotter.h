#pragma once

#include <opencv2/opencv.hpp>
#include <string>
#include <vector>
#include <algorithm>

struct Line {
	std::vector<std::pair<int, int>> points;
	cv::Scalar color;
};

const int step_he = 40;
const int step_wi = 60;
const int y_line = 50;
const int x_line = 50;
const int text_e = 25;
const int EPS = 5;
const int THIN2 = 2;

void plot_graph(int minY, int maxY, int stepY, int stepsX, int stepX, std::vector<Line> & plots, std::string fn) {
	int h = (maxY - minY + 4 * stepY) / stepY * step_he;
	int w = stepsX * step_wi + y_line + x_line;

	cv::Mat mat(h, w, CV_8UC3, cv::Scalar(50, 50, 50));

	cv::arrowedLine(mat, cv::Point(15, h), cv::Point(15, 10), cvScalar(255, 255, 0), THIN2, 8, 0, 0.02);
	int cnt = 0;
	int zro_cnt = 0;
	for (int i = minY; i <= maxY; i += stepY) {
		std::string text = std::to_string(i);
		if (i == 0) {
			cv::arrowedLine(mat, cv::Point(15, h - cnt * step_he - 2 * step_he), cv::Point(w - 15, h - cnt * step_he - 2 * step_he), cvScalar(255, 255, 0), THIN2, 8, 0, 0.02);
			zro_cnt = cnt;

			cv::putText(mat, text, cvPoint(text_e, h - cnt * step_he - 2 * step_he - EPS), cv::FONT_HERSHEY_TRIPLEX, 0.4, cvScalar(255, 255, 255), 1, CV_AA);
		} 
		else {
			cv::putText(mat, text, cvPoint(text_e, h - cnt * step_he - 2 * step_he), cv::FONT_HERSHEY_TRIPLEX, 0.4, cvScalar(255, 255, 255), 1, CV_AA);
		}

		cnt++;
	}

	for (int i = 1; i <= stepsX; i++) {
		std::string text = std::to_string(i * stepX);
		cv::putText(mat, text, cvPoint(text_e + step_wi * i, h - zro_cnt * step_he - 2 * step_he - EPS), cv::FONT_HERSHEY_TRIPLEX, 0.4, cvScalar(255, 255, 255), 1, CV_AA);

		cnt++;
	}

	for (Line &l : plots) {
		bool first = true;
		std::pair<int, int> last;
		for (std::pair<int, int> pnt : l.points) {
			pnt.second += (abs(minY) / stepY);
			if (first) {
				last = pnt;
				first = false;
			}
			else {
				cv::line(mat, cv::Point(text_e + last.first * step_wi, h - last.second * step_he - 2 * step_he), cv::Point(text_e + pnt.first * step_wi, h - pnt.second * step_he - 2 * step_he), l.color, THIN2);
				last = pnt;
			}
		}
	}

	imwrite(fn + ".png", mat);
}
