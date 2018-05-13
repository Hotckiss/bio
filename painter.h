#pragma once

#include <opencv2/opencv.hpp>
#include <string>
#include <vector>

const int HEIGHT_SPACE = 20;
const int HEIGHT_HALF = HEIGHT_SPACE / 2;
const int H_EPS = 10;
const int THIN = 2;
const int X_EPS = 6;
const int LEN = 6;
const int R = 3;
struct DetectionEntry {
	int pos;
	int type;      // false - prefix; true - suffix
	int lose_type; // 0 - noth, 1 -- nh3, 2-- h2o, 3 -- nh3+h2o
};

std::string remove_non_amino(std::string pept) {
	while (true) {
		size_t pos = pept.find("+57.021");
		if (pos != std::string::npos) {
			pept = pept.substr(0, pos) + pept.substr(pos + 7, pept.length());
		}
		else {
			break;
		}
	}

	return pept;
}
std::string convert_pep(std::string pep) {
	std::string ans;
	for (char c : pep) {
		ans += c;
		ans += " ";
	}

	if (ans.length() > 0)
		ans = ans.substr(0, ans.length() - 1);
	return ans;
}

cv::Size position_of_base(int pos, std::string clear_pep) {
	int baseline = 0;

	return cv::getTextSize(convert_pep(clear_pep.substr(0, pos + 1)), cv::FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);
}

void make_picture_with_ions(std::string pep, std::vector<DetectionEntry> &cid_cuts, std::vector<DetectionEntry> &hcd_cuts, int cid_scan_num, int hcd_scan_num) {
	std::string clear_pep = remove_non_amino(pep);
	std::string text_pep = convert_pep(clear_pep);

	int baseline = 0;
	cv::Size s = cv::getTextSize(text_pep, cv::FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);
	int h = s.height + HEIGHT_SPACE;
	int w = s.width;

	cv::Mat mat(h, w, CV_8UC3, cv::Scalar(50, 50, 50));

	cv::putText(mat, text_pep, cvPoint(0, h - HEIGHT_HALF), cv::FONT_HERSHEY_TRIPLEX, 1.5, cvScalar(255, 255, 255), 1, CV_AA); // text - while

	for (DetectionEntry pos : hcd_cuts) {
		cv::Size s1 = position_of_base(pos.pos, clear_pep);
		cv::Size s2 = cv::getTextSize(" ", cv::FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);

		int x = s1.width + s2.width / 2;
		line(mat, cv::Point(x - X_EPS, 0 + H_EPS), cv::Point(x - X_EPS, h - H_EPS), cvScalar(0, 255, 204), THIN); // hcd - yellow
		if (pos.type == false) {
			line(mat, cv::Point(x - X_EPS + LEN, h - H_EPS), cv::Point(x - X_EPS, h - H_EPS), cvScalar(0, 255, 204), THIN);
			if (pos.lose_type == 1) {
				cv::circle(mat, cv::Point(x - X_EPS + LEN / 2, h - H_EPS / 3), R, cv::Scalar(255, 110, 0), CV_FILLED);
			} else if (pos.lose_type == 2) {
				cv::circle(mat, cv::Point(x - X_EPS + LEN / 2, h - H_EPS / 3), R, cv::Scalar(0, 255, 0), CV_FILLED);
			} else if (pos.lose_type == 3) {
				cv::circle(mat, cv::Point(x - X_EPS + LEN / 2, h - H_EPS / 3), R, cv::Scalar(0, 0, 255), CV_FILLED);
			}
		}
		else {
			line(mat, cv::Point(x - X_EPS, 0 + H_EPS), cv::Point(x - X_EPS - LEN, 0 + H_EPS), cvScalar(0, 255, 204), THIN);
			if (pos.lose_type == 1) {
				cv::circle(mat, cv::Point(x - X_EPS - LEN / 2, 0 +  H_EPS / 3), R, cv::Scalar(255, 110, 0), CV_FILLED);
			} else if (pos.lose_type == 2) {
				cv::circle(mat, cv::Point(x - X_EPS - LEN / 2, 0 + H_EPS / 3), R, cv::Scalar(0, 255, 0), CV_FILLED);
			} else if (pos.lose_type == 3) {
				cv::circle(mat, cv::Point(x - X_EPS - LEN / 2, 0 + H_EPS / 3), R, cv::Scalar(0, 0, 255), CV_FILLED);
			}
		}
	}

	for (DetectionEntry pos : cid_cuts) {
		cv::Size s1 = position_of_base(pos.pos, clear_pep);
		cv::Size s2 = cv::getTextSize(" ", cv::FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);

		int x = s1.width + s2.width / 2;
		line(mat, cv::Point(x + X_EPS, 0 + H_EPS), cv::Point(x + X_EPS, h - H_EPS), cvScalar(255, 51, 204), THIN); // cid - pink
		if (pos.type == false) {
			line(mat, cv::Point(x + X_EPS + LEN, h - H_EPS), cv::Point(x + X_EPS, h - H_EPS), cvScalar(255, 51, 204), THIN);
			// green - lose nh3, red - lose both blue - lose water
			if (pos.lose_type == 1) {
				cv::circle(mat, cv::Point(x + X_EPS + LEN / 2, h - H_EPS / 3), R, cv::Scalar(255, 110, 0), CV_FILLED);
			} else if (pos.lose_type == 2) {
				cv::circle(mat, cv::Point(x + X_EPS + LEN / 2, h - H_EPS / 3), R, cv::Scalar(0, 255, 0), CV_FILLED);
			} else if (pos.lose_type == 3) {
				cv::circle(mat, cv::Point(x + X_EPS + LEN / 2, h - H_EPS / 3), R, cv::Scalar(0, 0, 255), CV_FILLED);
			}
		}
		else {
			line(mat, cv::Point(x + X_EPS, 0 + H_EPS), cv::Point(x + X_EPS - LEN, 0 + H_EPS), cvScalar(255, 51, 204), THIN);
			if (pos.lose_type == 1) {
				cv::circle(mat, cv::Point(x + X_EPS - LEN / 2, 0 + H_EPS / 3), R, cv::Scalar(255, 110, 0), CV_FILLED);
			} else if (pos.lose_type == 2) {
				cv::circle(mat, cv::Point(x + X_EPS - LEN / 2, 0 + H_EPS / 3), R, cv::Scalar(0, 255, 0), CV_FILLED);
			} else if (pos.lose_type == 3) {
				cv::circle(mat, cv::Point(x + X_EPS - LEN / 2, 0 + H_EPS / 3), R, cv::Scalar(0, 0, 255), CV_FILLED);
			}
		}
	}

	std::string pic_name = std::to_string(hcd_scan_num) + "_" + std::to_string(cid_scan_num) + std::string(".png");

	imwrite(pic_name, mat);
}

void make_picture(std::string pep, std::vector<int> &cid_cuts, std::vector<int> &hcd_cuts, int cid_scan_num, int hcd_scan_num) {
	std::string clear_pep = remove_non_amino(pep);
	std::string text_pep = convert_pep(clear_pep);

	int baseline = 0;
	cv::Size s = cv::getTextSize(text_pep, cv::FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);
	int h = s.height + HEIGHT_SPACE;
	int w = s.width;

	cv::Mat mat(h, w, CV_8UC3, cv::Scalar(50, 50, 50));

	cv::putText(mat, text_pep, cvPoint(0, h - HEIGHT_HALF), cv::FONT_HERSHEY_TRIPLEX, 1.5, cvScalar(255, 255, 255), 1, CV_AA); // text - while

	for (int pos : hcd_cuts) {
		cv::Size s1 = position_of_base(pos, clear_pep);
		cv::Size s2 = cv::getTextSize(" ", cv::FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);

		int x = s1.width + s2.width / 2;
		line(mat, cv::Point(x - X_EPS, 0 + H_EPS), cv::Point(x - X_EPS, h - H_EPS), cvScalar(0, 255, 204), THIN); // hcd - yellow
	}

	for (int pos : cid_cuts) {
		cv::Size s1 = position_of_base(pos, clear_pep);
		cv::Size s2 = cv::getTextSize(" ", cv::FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);

		int x = s1.width + s2.width / 2;
		line(mat, cv::Point(x + X_EPS, 0 + H_EPS), cv::Point(x + X_EPS, h - H_EPS), cvScalar(255, 51, 204), THIN); // cid - pink
	}

	std::string pic_name = std::to_string(hcd_scan_num) + "_" + std::to_string(cid_scan_num) + std::string(".png");

	imwrite(pic_name, mat);
}
