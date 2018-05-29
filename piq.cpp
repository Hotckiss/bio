#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include "painter.h"
#include "gister.h"

using namespace std;

class peptide {
public:
	int scan_number;
	string frag_method;
	long double precursor_mass, precursor_charge;
	vector<long double> actual_masses;
};

const long double eps_cid = 2e-1;
const long double eps_hcd = 2e-2;
const long double val = 57.021;
const long double h2o = 18.010157;
const long double nh3 = 17.030215;

vector<peptide> all_peptides;
map<int, string> pept_by_scan;
map<int, long double> evalue_by_scan;
map<char, long double> amino_mass;
vector<pair<int, string>> listed_scans;

vector<long double> get_prefix_sums(string pept) {
	size_t pos = pept.find("+57.021");
	if (pos != string::npos) {
		pept = pept.substr(0, pos) + pept.substr(pos + 7, pept.length());
	}

	vector<long double> ans;
	pos--;
	if (pept.empty())
		return ans;
	ans.push_back(amino_mass[pept[0]]);
	if (pos == 0) {
		ans[ans.size() - 1] += val;
	}
	for (int i = 1; i < pept.length(); i++) {
		ans.push_back(ans[ans.size() - 1] + amino_mass[pept[i]]);
		if (pos == i) {
			ans[ans.size() - 1] += val;
		}
	}

	return ans; //18.010157
}

vector<long double> get_suffix_sums(string pept) {
	size_t pos = pept.find("+57.021");
	if (pos != string::npos) {
		pept = pept.substr(0, pos) + pept.substr(pos + 7, pept.length());
	}

	vector<long double> ans;
	pos--;
	if (pept.empty())
		return ans;
	ans.push_back(amino_mass[pept[pept.length() - 1]]);
	if (pos == pept.length() - 1) {
		ans[ans.size() - 1] += val;
	}
	for (int i = pept.length() - 2; i >= 0; i--) {
		ans.push_back(ans[ans.size() - 1] + amino_mass[pept[i]]);
		if (pos == i) {
			ans[ans.size() - 1] += val;
		}
	}

	return ans;
}
void init() {
	amino_mass['A'] = 71.03711;
	amino_mass['R'] = 156.10111;
	amino_mass['N'] = 114.04293;
	amino_mass['D'] = 115.02694;
	amino_mass['C'] = 103.00919;
	amino_mass['E'] = 129.04259;
	amino_mass['Q'] = 128.05858;
	amino_mass['G'] = 57.02146;
	amino_mass['H'] = 137.05891;
	amino_mass['I'] = 113.08406;
	amino_mass['L'] = 113.08406;
	amino_mass['K'] = 128.09496;
	amino_mass['M'] = 131.04049;
	amino_mass['F'] = 147.06841;
	amino_mass['P'] = 97.05276;
	amino_mass['S'] = 87.03203;
	amino_mass['T'] = 101.04768;
	amino_mass['W'] = 186.07931;
	amino_mass['Y'] = 163.06333;
	amino_mass['V'] = 99.06841;
}

void parse() {
	while (true) {
		string str;
		getline(cin, str, '\n');

		if (str == "BEGIN IONS") {
			getline(cin, str, '\n');
			getline(cin, str, '\n');
			str = str.substr(6, str.length());

			peptide p;
			p.scan_number = atoi(str.c_str());
			getline(cin, str, '\n');
			str = str.substr(11, str.length());
			p.frag_method = str;

			getline(cin, str, '\n');
			getline(cin, str, '\n');
			str = str.substr(17, str.length());
			p.precursor_charge = atoi(str.c_str());

			getline(cin, str, '\n');
			str = str.substr(15, str.length());
			p.precursor_mass = stod(str.c_str());

			bool is_end = false;
			while (!is_end) {
				getline(cin, str, '\n');
				if (str == "END IONS") {
					is_end = true;
					continue;
				}

				char * pch = strtok(const_cast<char*>(str.c_str()), "\t");

				p.actual_masses.push_back(stod(pch));
			}
			getline(cin, str, '\n');

			all_peptides.push_back(p);
		}
		else {
			break;
		}
	}
}

void fill_pept_by_scan() {
	ifstream input("table2.txt", ifstream::in);

	int n;
	input >> n;

	for (int i = 0; i < n; i++) {
		int sn;
		string pept;
		long double evalue;
		input >> sn >> pept >> evalue;

		listed_scans.push_back({sn, pept});
		pept_by_scan[sn] = pept;
		evalue_by_scan[sn] = evalue;
	}

	input.close();
}

inline bool has_small_evalue(int scan) {
	return evalue_by_scan.find(scan) != evalue_by_scan.end();
}

void resume_annotated() {
	cout << "SELECTED DETECTION ERROR = " << eps_hcd << " for HCD" << "; and " << eps_cid << " for CID\n\n";
	for (peptide &p : all_peptides) {
		string pept(pept_by_scan[p.scan_number]);
		if (pept.empty() || !has_small_evalue(p.scan_number))
			continue;
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);
		vector<int> detected_pref;
		vector<int> detected_suf;
		vector<double> detected_pref1;
		vector<double> detected_suf1;
		double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
		for (long double mass : p.actual_masses) {
			for (int i = 0; i < pref_sums.size(); i++)
				if (abs(mass - pref_sums[i]) < eps) {
					detected_pref.push_back(i);
					detected_pref1.push_back(mass);
				}

			for (int i = 0; i < suf_sums.size(); i++) // 17.021 - N; 18.010157 - H2O
				if (abs(mass - suf_sums[i] - 18.010157) < eps) {
					detected_suf.push_back(i);
					detected_suf1.push_back(mass);
				}
		}
		if (detected_pref.size() + detected_suf.size() > 0) {
			cout << "Detected for peptide " << pept << ":\n" << "SCAN = " << p.scan_number << "; FRAG METHOD = " << p.frag_method << "; PRECURSOR MASS = " << p.precursor_mass << "; CHARGE = " << p.precursor_charge << "\n\n";

			size_t pos = pept.find("+57.021");
			if (pos != string::npos) {
				pept = pept.substr(0, pos) + pept.substr(pos + 7, pept.length());
			}

			for (int i = 0; i < detected_pref.size(); i++) {
				cout << "Prefix " << pept.substr(0, detected_pref[i] + 1) << "; MASS = " << pref_sums[detected_pref[i]] << "; ACTUAL = " << detected_pref1[i] << "; DELTA = " << pref_sums[detected_pref[i]] - detected_pref1[i] << "\n";
			}

			for (int i = 0; i < detected_suf.size(); i++) {
				cout << "Suffix " << pept.substr(pept.length() - 1 - detected_suf[i], pept.length()) << "; MASS = " << suf_sums[detected_suf[i]] << "; ACTUAL = " << detected_suf1[i] << "; DELTA = " << suf_sums[detected_suf[i]] + 18.010157 - detected_suf1[i] << "\n";
			}

			cout << "------------------------\n\n\n";
		}
	}
}

void build_gists_equal_length() {
	vector<vector<int>> hcds, cids;
	hcds.resize(100);
	cids.resize(100);
	for (int i = 0; i < 100; i++) {
		hcds[i].resize(i + 1);
		cids[i].resize(i + 1);
	}

	for (int i = 0; i < all_peptides.size(); i++) {
		peptide p1 = all_peptides[i];

		if (!has_small_evalue(p1.scan_number))
			continue;

		//amino sequence
		string pept(pept_by_scan[p1.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);

		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		// checking that jth position annotated
		// pref[i] = mass of prefix of length i + 1;
		// suff[i] = mass of suffix of length i + 1;
		vector<DetectionEntry> v1; // annotated positions in p1;
		vector<int> vv1;
		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;
			
			double eps = (p1.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p1.actual_masses) {
				if (abs(mass - pref_m) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
				else if (abs(mass - pref_m_h2o) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
				else if (abs(mass - pref_m_nh3) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
				else if (abs(mass - pref_m_h2o_nh3) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
				else if (abs(mass - suf_m) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
				else if (abs(mass - suf_m_h2o) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
				else if (abs(mass - suf_m_nh3) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
				else if (abs(mass - suf_m_h2o_nh3) < eps) {
					if (p1.frag_method == "HCD") hcds[clear_pept.length()][j]++;
					else cids[clear_pept.length()][j]++;
				}
			}
		}
	}

	for (int i = 1; i < 100; i++) {
		int s = 0;
		for (int j = 0; j < hcds[i].size(); j++) {
			s += hcds[i][j];
			s += cids[i][j];
		}

		if (s == 0) continue;

		make_gist(i - 1, hcds[i], cids[i], "gist_" + to_string(i));
	}
}

void build_gists_close_length(int lb, int ub) {
	vector<int> vals_hcd, vals_cid;
	vals_cid.resize(10);
	vals_hcd.resize(10);
	for (int i = 0; i < all_peptides.size(); i++) {
		peptide p1 = all_peptides[i];

		if (!has_small_evalue(p1.scan_number))
			continue;

		//amino sequence
		string pept(pept_by_scan[p1.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);
		if (!(clear_pept.length() >= lb && clear_pept.length() <= ub)) continue;
		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		vector<int> vv1;
		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;

			double eps = (p1.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p1.actual_masses) {
				if (abs(mass - pref_m) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
				else if (abs(mass - pref_m_h2o) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
				else if (abs(mass - pref_m_nh3) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
				else if (abs(mass - pref_m_h2o_nh3) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
				else if (abs(mass - suf_m) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
				else if (abs(mass - suf_m_h2o) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
				else if (abs(mass - suf_m_nh3) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
				else if (abs(mass - suf_m_h2o_nh3) < eps) {
					if (p1.frag_method == "HCD") vals_hcd[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
					else vals_cid[(int)(10.0 * j / ((clear_pept.length() - 1)))]++;
				}
			}
		}
	}

	make_gist(10, vals_hcd, vals_cid, "close_gist_" + to_string(lb) + "_" + to_string(ub));
}

void build_annotated_positions() {
	ofstream ofs("percentage.txt", ofstream::out);

	for (int i = 0; i < all_peptides.size() - 1; i++) {
		peptide p1 = all_peptides[i];
		peptide p2 = all_peptides[i + 1];

		if (p1.scan_number + 1 != p2.scan_number)
			continue;

		if (!(has_small_evalue(p1.scan_number) && has_small_evalue(p2.scan_number)))
			continue;

		if (pept_by_scan[p1.scan_number] != pept_by_scan[p2.scan_number])
			continue;

		//else we have pair of scans with good evalue and same peptides

		//amino sequence
		string pept(pept_by_scan[p1.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);

		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		// checking that jth position annotated
		// pref[i] = mass of prefix of length i + 1;
		// suff[i] = mass of suffix of length i + 1;
		vector<DetectionEntry> v1; // annotated positions in p1;
		vector<int> vv1;
		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;

			double eps = (p1.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p1.actual_masses) {
				if (abs(mass - pref_m) < eps) {
					v1.push_back({ j, false, 0 });
					vv1.push_back(j);
				}
				else if (abs(mass - pref_m_h2o) < eps) {
					v1.push_back({ j, false, 1 });
					vv1.push_back(j);
				}
				else if (abs(mass - pref_m_nh3) < eps) {
					v1.push_back({ j, false, 2 });
					vv1.push_back(j);
				}
				else if (abs(mass - pref_m_h2o_nh3) < eps) {
					v1.push_back({ j, false, 3 });
					vv1.push_back(j);
				}
				else if (abs(mass - suf_m) < eps) {
					v1.push_back({ j, true, 0 });
					vv1.push_back(j);
				}
				else if (abs(mass - suf_m_h2o) < eps) {
					v1.push_back({ j, true, 1 });
					vv1.push_back(j);
				}
				else if (abs(mass - suf_m_nh3) < eps) {
					v1.push_back({ j, true, 2 });
					vv1.push_back(j);
				}
				else if (abs(mass - suf_m_h2o_nh3) < eps) {
					v1.push_back({ j, true, 3 });
					vv1.push_back(j);
				}
			}
		}

		vector<DetectionEntry> v2; // annotated positions in p1;
		vector<int> vv2;
		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;

			double eps = (p2.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p2.actual_masses) {
				if (abs(mass - pref_m) < eps) {
					v2.push_back({ j, false, 0 });
					vv2.push_back(j);
				}
				else if (abs(mass - pref_m_h2o) < eps) {
					v2.push_back({ j, false, 1 });
					vv2.push_back(j);
				}
				else if (abs(mass - pref_m_nh3) < eps) {
					v2.push_back({ j, false, 2 });
					vv2.push_back(j);
				}
				else if (abs(mass - pref_m_h2o_nh3) < eps) {
					v2.push_back({ j, false, 3 });
					vv2.push_back(j);
				}
				else if (abs(mass - suf_m) < eps) {
					v2.push_back({ j, true, 0 });
					vv2.push_back(j);
				}
				else if (abs(mass - suf_m_h2o) < eps) {
					v2.push_back({ j, true, 1 });
					vv2.push_back(j);
				}
				else if (abs(mass - suf_m_nh3) < eps) {
					v2.push_back({ j, true, 2 });
					vv2.push_back(j);
				}
				else if (abs(mass - suf_m_h2o_nh3) < eps) {
					v2.push_back({ j, true, 3 });
					vv2.push_back(j);
				}
			}
		}

		ofs << "Pair of scans for peptide " << pept << "\n";
		if (p1.frag_method == "HCD" && p2.frag_method == "CID") {
			int h = std::unique(vv1.begin(), vv1.end()) - vv1.begin();
			int c = std::unique(vv2.begin(), vv2.end()) - vv2.begin();
			ofs << "For HCD precantage of annotated cuts = " << 100.0 * h / (clear_pept.length() - 1) << "%\n";
			ofs << "For CID precantage of annotated cuts = " << 100.0 * c / (clear_pept.length() - 1) << "%\n";
			make_picture_with_ions(pept, v2, v1, p2.scan_number, p1.scan_number);
		}
		else if (p2.frag_method == "HCD" && p1.frag_method == "CID") {
			int c = std::unique(vv1.begin(), vv1.end()) - vv1.begin();
			int h = std::unique(vv2.begin(), vv2.end()) - vv2.begin();
			ofs << "For HCD precantage of annotated cuts = " << 100.0 * h / (clear_pept.length() - 1) << "%\n";
			ofs << "For CID precantage of annotated cuts = " << 100.0 * c / (clear_pept.length() - 1) << "%\n";
			make_picture_with_ions(pept, v1, v2, p1.scan_number, p2.scan_number);
		}
		ofs << endl;

		i++;
	}
	ofs.close();
}

void calculate_percentage_of_detection_in_total_no_losses(ofstream &ofs);
void calculate_percentage_of_detection_in_total_losses_h2o(ofstream &ofs);
void calculate_percentage_of_detection_in_total_losses_nh3(ofstream &ofs);
void calculate_percentage_of_detection_in_total_losses_h2o_nh3(ofstream &ofs);
void calculate_percentage_of_detection_in_total(ofstream &ofs);

void calculate_percentage_of_detection_in_total_for_each(ofstream &ofs);
void calculate_percentage_of_detection_in_total_for_each_rel(ofstream &ofs);

int main() {
	freopen("piq.msalign", "r", stdin);
	freopen("annotations2.txt", "w", stdout);
	
	init();
	parse();
	fill_pept_by_scan();
	resume_annotated();
	build_annotated_positions();
	
	ofstream ofs("total_annotated.txt", ofstream::out);
	calculate_percentage_of_detection_in_total_no_losses(ofs);
	calculate_percentage_of_detection_in_total_losses_h2o(ofs);
	calculate_percentage_of_detection_in_total_losses_nh3(ofs);
	calculate_percentage_of_detection_in_total_losses_h2o_nh3(ofs);
	calculate_percentage_of_detection_in_total(ofs);
	ofs.close();
	
	ofstream ofs1("total_annotated_for_each_div_total.txt", ofstream::out);
	calculate_percentage_of_detection_in_total_for_each(ofs1);
	ofs1.close();
	
	ofstream ofs2("total_annotated_for_each_div_detected.txt", ofstream::out);
	calculate_percentage_of_detection_in_total_for_each_rel(ofs2);
	ofs2.close();
	
	build_gists_equal_length();
	
	build_gists_close_length(7, 9);
	build_gists_close_length(10, 13);
	build_gists_close_length(14, 17);
	build_gists_close_length(18, 21);
	build_gists_close_length(22, 24);
	build_gists_close_length(27, 30);
	build_gists_close_length(32, 35);
	build_gists_close_length(37, 39);
	return 0;
}













void calculate_percentage_of_detection_in_total_for_each_rel(ofstream &ofs) {
	map<string, set<int>> positions;
	map<string, set<int>> positions_cid;
	map<string, set<int>> positions_hcd;
	for (peptide p : all_peptides) {
		if (!(has_small_evalue(p.scan_number)))
			continue;
		string pept(pept_by_scan[p.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);

		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			positions[clear_pept];
			positions_hcd[clear_pept];
			positions_cid[clear_pept];
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;
			double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p.actual_masses) {
				if (abs(mass - pref_m) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - pref_m_h2o) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - pref_m_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - pref_m_h2o_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m_h2o) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m_h2o_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
			}
		}
	}

	ofs << "Assuming all avaliable detections!\n\n";

	for (pair<string, set<int>> pss : positions_hcd) {
		set<int> sss = positions[pss.first];
		ofs << "Total percentage of detected cuts using HCD for peptide " << pss.first << ": " << (double)pss.second.size() / (sss.size()) * 100 << "% (" << pss.second.size() << "/" << sss.size() << ")\n";
		set<int> & ss = positions_cid[pss.first];
		ofs << "Total percentage of detected cuts using CID for peptide " << pss.first << ": " << (double)ss.size() / (sss.size()) * 100 << "% (" << ss.size() << "/" << sss.size() << ")\n";
		ofs << endl;
	}
}

void calculate_percentage_of_detection_in_total_for_each(ofstream &ofs) {
	map<string, set<int>> positions;
	map<string, set<int>> positions_cid;
	map<string, set<int>> positions_hcd;
	for (peptide p : all_peptides) {
		if (!(has_small_evalue(p.scan_number)))
			continue;
		string pept(pept_by_scan[p.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);

		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			positions[clear_pept];
			positions_hcd[clear_pept];
			positions_cid[clear_pept];
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;
			double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p.actual_masses) {
				if (abs(mass - pref_m) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - pref_m_h2o) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - pref_m_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - pref_m_h2o_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m_h2o) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
				else if (abs(mass - suf_m_h2o_nh3) < eps) {
					positions[clear_pept].insert(j);
					if (p.frag_method == "HCD") {
						positions_hcd[clear_pept].insert(j);
					}
					else { //CID
						positions_cid[clear_pept].insert(j);
					}
				}
			}
		}
	}
	
	ofs << "Assuming all avaliable detections!\n\n";

	vector<int> vhcd, vcid;
	for (pair<string, set<int>> pss : positions_hcd) {
		ofs << "Total percentage of detected cuts using HCD for peptide " << pss.first << ": "<< (double)pss.second.size() / (pss.first.length() - 1) * 100 << "% (" << pss.second.size() << "/" << pss.first.length() - 1 << ")\n";
		set<int> & ss = positions_cid[pss.first];
		ofs << "Total percentage of detected cuts using CID for peptide " << pss.first << ": " << (double)ss.size() / (pss.first.length() - 1) * 100 << "% (" << ss.size() << "/" << pss.first.length() - 1 << ")\n";
		ofs << endl;

		vhcd.push_back(pss.second.size() * 100 / (pss.first.length() - 1));
		vcid.push_back(ss.size() * 100 / (pss.first.length() - 1));
	}
	make_gist(vhcd.size(), vhcd, vcid, "percentage_for_each_div_total");
}





typedef map<string, vector<int>> msvi;
typedef map<string, set<int>> mssi;
vector<long double> conclusion(ofstream &ofs, msvi& positions, msvi& positions_cid, msvi& positions_hcd, mssi& spositions, mssi& spositions_cid, mssi& spositions_hcd, int sz) {
	int tot = 0;
	//int sz = 0;
	for (pair<string, vector<int>> pss : positions) {
		tot += pss.second.size();
		//sz += pss.first.length() - 1;
	}

	vector<long double> hist;
	ofs << "Total percentage of detected cuts using both HCD and CID: " << (double)tot / sz * 100 << "% (" << tot << "/" << sz << ")\n";
	hist.push_back((double)tot / sz * 100);

	tot = 0;
	for (pair<string, vector<int>> pss : positions_cid) tot += pss.second.size();

	ofs << "Total percentage of detected cuts using CID: " << (double)tot / sz * 100 << "% (" << tot << "/" << sz << ")\n";
	hist.push_back((double)tot / sz * 100);

	tot = 0;
	for (pair<string, vector<int>> pss : positions_hcd) tot += pss.second.size();

	ofs << "Total percentage of detected cuts using HCD: " << (double)tot / sz * 100 << "% (" << tot << "/" << sz << ")\n";
	hist.push_back((double)tot / sz * 100);

	int tot_both = 0;
	int tot_cid = 0;
	int all = 0;
	for (pair<string, set<int>> ps : spositions_cid) {
		for (int x : ps.second)
			if (spositions_hcd[ps.first].find(x) == spositions_hcd[ps.first].end()) {
				tot_cid++;
			}
			else {
				tot_both++;
			}
	}
	for (pair<string, set<int>> ps : spositions) {
		all += ps.second.size();
	}

	ofs << "Total percentage of detected cuts using BOTH methods: " << (double)(tot_both) / all * 100 << "% (" << tot_both << "/" << all << ")\n";
	hist.push_back((double)(tot_both) / all * 100);

	ofs << "Total percentage of detected cuts with CID, but not HCD methods: " << (double)(tot_cid) / all * 100 << "% (" << tot_cid << "/" << all << ")\n";
	hist.push_back((double)(tot_cid) / all * 100);

	ofs << "Total percentage of detected cuts with HCD, but not CID methods: " << (double)(all - tot_cid - tot_both) / all * 100 << "% (" << all - tot_cid - tot_both << "/" << all << ")\n";
	hist.push_back((double)(all - tot_cid - tot_both) / all * 100);


	ofs << "\n\n";

	return hist;
}

void insert_all(msvi& positions, msvi& positions_cid, msvi& positions_hcd, mssi& spositions, mssi& spositions_cid, mssi& spositions_hcd, 
	long double mass, long double pref_m, long double suf_m, peptide &p, double eps, string clear_pept, int j) {
	if (abs(mass - pref_m) < eps) {
		positions[clear_pept].push_back(j);
		spositions[clear_pept].insert(j);
		if (p.frag_method == "HCD") {
			positions_hcd[clear_pept].push_back(j);
			spositions_hcd[clear_pept].insert(j);
		}
		else { //CID
			positions_cid[clear_pept].push_back(j);
			spositions_cid[clear_pept].insert(j);
		}
	}
	if (abs(mass - suf_m) < eps) {
		positions[clear_pept].push_back(j);
		spositions[clear_pept].insert(j);
		if (p.frag_method == "HCD") {
			positions_hcd[clear_pept].push_back(j);
			spositions_hcd[clear_pept].insert(j);
		}
		else { //CID
			positions_cid[clear_pept].push_back(j);
			spositions_cid[clear_pept].insert(j);
		}
	}
}

void init(msvi& positions, msvi& positions_cid, msvi& positions_hcd, mssi& spositions, mssi& spositions_cid, mssi& spositions_hcd, string clear_pept) {
	positions[clear_pept];
	positions_hcd[clear_pept];
	positions_cid[clear_pept];
	spositions[clear_pept];
	spositions_hcd[clear_pept];
	spositions_cid[clear_pept];
}

void calculate_percentage_of_detection_in_total_no_losses(ofstream &ofs) {
	msvi positions, positions_cid, positions_hcd;
	mssi spositions, spositions_cid, spositions_hcd;

	int sz = 0;
	for (peptide p : all_peptides) {
		if (!(has_small_evalue(p.scan_number))) continue;
		string pept(pept_by_scan[p.scan_number]);
		string clear_pept = remove_non_amino(pept);
		sz += clear_pept.length() - 1;
		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			init(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, clear_pept);
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;

			double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p.actual_masses) {
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m, suf_m, p, eps, clear_pept, j);
			}
		}
	}
	
	ofs << "Assuming no losses of H2O and NH3 was during experiment!\n\n";
	vector<long double> concl = conclusion(ofs, positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, sz);
	vector<long double> c1, c2;
	c1.push_back(concl[0]);
	c1.push_back(concl[1]);
	c1.push_back(concl[2]);
	c2.push_back(concl[3]);
	c2.push_back(concl[4]);
	c2.push_back(concl[5]);
	make_gist_hex_percentage1(c1, "hhex_no_loses_1");
	make_gist_hex_percentage2(c2, "hhex_no_loses_2");
}

void calculate_percentage_of_detection_in_total_losses_h2o(ofstream &ofs) {
	msvi positions, positions_cid, positions_hcd;
	mssi spositions, spositions_cid, spositions_hcd;
	int sz = 0;
	for (peptide p : all_peptides) {
		if (!(has_small_evalue(p.scan_number))) continue;
		string pept(pept_by_scan[p.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);
		sz += clear_pept.length() - 1;
		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			init(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, clear_pept);
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;

			long double suf_m_h2o = suf_m - h2o;
			double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p.actual_masses) {
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m_h2o, suf_m_h2o, p, eps, clear_pept, j);
			}
		}
	}

	ofs << "Assuming H2O was lost!\n\n";
	vector<long double> concl = conclusion(ofs, positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, sz);
	vector<long double> c1, c2;
	c1.push_back(concl[0]);
	c1.push_back(concl[1]);
	c1.push_back(concl[2]);
	c2.push_back(concl[3]);
	c2.push_back(concl[4]);
	c2.push_back(concl[5]);
	make_gist_hex_percentage1(c1, "hhex_loses_h2o_1");
	make_gist_hex_percentage2(c2, "hhex_loses_h2o_2");
}

void calculate_percentage_of_detection_in_total_losses_nh3(ofstream &ofs) {
	msvi positions, positions_cid, positions_hcd;
	mssi spositions, spositions_cid, spositions_hcd;
	int sz = 0;
	for (peptide p : all_peptides) {
		if (!(has_small_evalue(p.scan_number))) continue;
		string pept(pept_by_scan[p.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);
		sz += clear_pept.length() - 1;
		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			init(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, clear_pept);
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_nh3 = pref_m - nh3;
			long double suf_m_nh3 = suf_m - nh3;
			double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p.actual_masses) {
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m_nh3, suf_m_nh3, p, eps, clear_pept, j);
			}
		}
	}

	ofs << "Assuming NH3 was lost!\n\n";
	vector<long double> concl = conclusion(ofs, positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, sz);
	vector<long double> c1, c2;
	c1.push_back(concl[0]);
	c1.push_back(concl[1]);
	c1.push_back(concl[2]);
	c2.push_back(concl[3]);
	c2.push_back(concl[4]);
	c2.push_back(concl[5]);
	make_gist_hex_percentage1(c1, "hhex_loses_nh3_1");
	make_gist_hex_percentage2(c2, "hhex_loses_nh3_2");
}

void calculate_percentage_of_detection_in_total_losses_h2o_nh3(ofstream &ofs) {
	msvi positions, positions_cid, positions_hcd;
	mssi spositions, spositions_cid, spositions_hcd;
	int sz = 0;
	for (peptide p : all_peptides) {
		if (!(has_small_evalue(p.scan_number)))
			continue;
		string pept(pept_by_scan[p.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);
		sz += clear_pept.length() - 1;
		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			init(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, clear_pept);
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;
			double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p.actual_masses) {
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m_h2o_nh3, suf_m_h2o_nh3, p, eps, clear_pept, j);
			}
		}
	}

	ofs << "Assuming H2O and NH3 was lost!\n\n";
	vector<long double> concl = conclusion(ofs, positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, sz);
	vector<long double> c1, c2;
	c1.push_back(concl[0]);
	c1.push_back(concl[1]);
	c1.push_back(concl[2]);
	c2.push_back(concl[3]);
	c2.push_back(concl[4]);
	c2.push_back(concl[5]);
	make_gist_hex_percentage1(c1, "hhex_loses_nh3_h2o_1");
	make_gist_hex_percentage2(c2, "hhex_loses_nh3_h2o_2");
}

void calculate_percentage_of_detection_in_total(ofstream &ofs) {
	msvi positions, positions_cid, positions_hcd;
	mssi spositions, spositions_cid, spositions_hcd;
	int sz = 0;
	for (peptide p : all_peptides) {
		if (!(has_small_evalue(p.scan_number)))
			continue;
		string pept(pept_by_scan[p.scan_number]);

		//peptide without +57.021
		string clear_pept = remove_non_amino(pept);
		sz += clear_pept.length() - 1;
		//all perf sums without losses h2o / nh3
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);

		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			init(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, clear_pept);
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;
			double eps = (p.frag_method == "HCD") ? eps_hcd : eps_cid;
			for (long double mass : p.actual_masses) {
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m, suf_m, p, eps, clear_pept, j);
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m_h2o, suf_m_h2o, p, eps, clear_pept, j);
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m_nh3, suf_m_nh3, p, eps, clear_pept, j);
				insert_all(positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, mass, pref_m_h2o_nh3, suf_m_h2o_nh3, p, eps, clear_pept, j);
			}
		}
	}

	ofs << "Assuming all avaliable detections!\n\n";
	vector<long double> concl = conclusion(ofs, positions, positions_cid, positions_hcd, spositions, spositions_cid, spositions_hcd, sz);
	vector<long double> c1, c2;
	c1.push_back(concl[0]);
	c1.push_back(concl[1]);
	c1.push_back(concl[2]);
	c2.push_back(concl[3]);
	c2.push_back(concl[4]);
	c2.push_back(concl[5]);
	make_gist_hex_percentage1(c1, "hhex_total_1");
	make_gist_hex_percentage2(c2, "hhex_total_2");
}
