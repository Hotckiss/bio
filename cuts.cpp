#include <opencv2/opencv.hpp>
#include <string>
#include <vector>
using namespace std;
using namespace cv;

const int HEIGHT_SPACE = 20;
const int HEIGHT_HALF = HEIGHT_SPACE / 2;
const int H_EPS = 5;
const int THIN = 4;
const int X_EPS = 4;
string remove_non_amino(string pept) {
	while (true) {
		size_t pos = pept.find("+57.021");
		if (pos != string::npos) {
			pept = pept.substr(0, pos) + pept.substr(pos + 7, pept.length());
		}
		else {
			break;
		}
	}

	return pept;
}
string convert_pep(string pep) {
	string ans;
	for (char c : pep) {
		ans += c;
		ans += " ";
	}

	if (ans.length() > 0)
		ans = ans.substr(0, ans.length() - 1);
	return ans;
}

Size position_of_base(int pos, string clear_pep) {
	int baseline = 0;

	return getTextSize(convert_pep(clear_pep.substr(0, pos + 1)), FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);
}

void make_picture(string pep, vector<int> &cid_cuts, vector<int> &hcd_cuts, int cid_scan_num, int hcd_scan_num) {
	string clear_pep = remove_non_amino(pep);
	string text_pep = convert_pep(clear_pep);

	int baseline = 0;
	Size s = getTextSize(text_pep, FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);
	int h = s.height + HEIGHT_SPACE;
	int w = s.width;

	Mat mat(h, w, CV_8UC3, Scalar(50, 50, 50));

	putText(mat, text_pep, cvPoint(0, h - HEIGHT_HALF), FONT_HERSHEY_TRIPLEX, 1.5, cvScalar(255, 255, 255), 1, CV_AA); // text - while

	for (int pos : hcd_cuts) {
		Size s1 = position_of_base(pos, clear_pep);
		Size s2 = getTextSize(" ", FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);

		int x = s1.width + s2.width / 2;
		line(mat, Point(x - X_EPS, 0 + H_EPS), Point(x - X_EPS, h - H_EPS) , cvScalar(0, 255, 204), THIN); // hcd - yellow
	}

	for (int pos : cid_cuts) {
		Size s1 = position_of_base(pos, clear_pep);
		Size s2 = getTextSize(" ", FONT_HERSHEY_TRIPLEX, 1.5, 1, &baseline);

		int x = s1.width + s2.width / 2;
		line(mat, Point(x + X_EPS, 0 + H_EPS), Point(x + X_EPS, h - H_EPS), cvScalar(255, 51, 204), THIN); // cid - pink
	}

	string pic_name = to_string(hcd_scan_num) + "_" + to_string(cid_scan_num) + string(".png");

	imwrite(pic_name, mat);
}
