#pragma once

#include <opencv2/opencv.hpp>
#include <string>
#include <vector>
#include <algorithm>

const int column_width = 16;
const int column_step_inner = 4;
const int column_step_outer = 10;
const int h_step = 20;
const int h_eps = 100;
const int x_eps = 20;
const int start_eps = column_width / 2;
const int Y_axis = 20;
const int text_eps = 3;
const int text_eps_Y = 4;
const int text_eps_x = 20;
int total_width(int size) {
	return size * (2 * column_width + column_step_inner + column_step_outer) - column_step_outer;
}

int total_height(int num) {
	return h_step * num;
}

int max_el(std::vector<int> &hcd_vals, std::vector<int> &cid_vals) {
	int res = 0;

	for (int x : hcd_vals)
		res = std::max(x, res);
	for (int x : cid_vals)
		res = std::max(x, res);

	return res;
}

void make_gist(int size, std::vector<int> &hcd_vals, std::vector<int> &cid_vals, std::string fn) {
	int max_h = max_el(hcd_vals, cid_vals);

	int baseline = 0;
	int h = total_height(max_h) + h_eps;
	int w = total_width(size) + 2 * x_eps + Y_axis;

	cv::Mat mat(h, w, CV_8UC3, cv::Scalar(50, 50, 50));

	
	for (int i = 0; i <= max_h / 5 + 1; i++) {
		std::string text = std::to_string(i * 5);
		cv::putText(mat, text, cvPoint(text_eps, h - i * 5 * h_step - text_eps_Y- text_eps_x), cv::FONT_HERSHEY_TRIPLEX, 0.4, cvScalar(255, 255, 255), 1, CV_AA);
	}

	for (int i = 0; i < size; i++) {
		std::string text = std::to_string(i + 1);
		int x = Y_axis + x_eps + start_eps + i * (2 * column_width + column_step_inner + column_step_outer) + column_step_inner;
		cv::putText(mat, text, cvPoint(x, h - text_eps_x / 4), cv::FONT_HERSHEY_TRIPLEX, 0.4, cvScalar(255, 255, 255), 1, CV_AA);
	}

	for (int i = 0; i < size; i++) {
		int x = Y_axis + x_eps + start_eps + i * (2 * column_width + column_step_inner + column_step_outer);
		int hh = h_step * hcd_vals[i];
		
		cv::rectangle(mat, cv::Point(x - column_width / 2, h - text_eps_x), cv::Point(x + column_width / 2, h - hh - text_eps_x), cvScalar(0, 255, 204), -1);
	}

	for (int i = 0; i < size; i++) {
		int x = Y_axis + x_eps + start_eps + i * (2 * column_width + column_step_inner + column_step_outer) + column_width + column_step_inner;
		int hh = h_step * cid_vals[i];
		cv::rectangle(mat, cv::Point(x - column_width / 2, h - text_eps_x), cv::Point(x + column_width / 2, h - hh - text_eps_x), cvScalar(255, 51, 204), -1);
	}

	imwrite(fn + ".png", mat);

}
