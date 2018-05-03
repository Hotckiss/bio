#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include "painter.h"

using namespace std;

class peptide {
public:
	int scan_number;
	string frag_method;
	long double precursor_mass, precursor_charge;
	vector<long double> actual_masses;
};

const long double eps = 2e-1;
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
	cout << "SELECTED DETECTION ERROR = " << eps << "\n\n";
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
		vector<int> v1; // annotated positions in p1;
		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;
			
			bool annotated = false;
			for (long double mass : p1.actual_masses) {
				if ((abs(mass - pref_m) < eps) || (abs(mass - pref_m_h2o) < eps) || (abs(mass - pref_m_nh3) < eps) || (abs(mass - pref_m_h2o_nh3) < eps) ||
					(abs(mass - suf_m) < eps) || (abs(mass - suf_m_h2o) < eps) || (abs(mass - suf_m_nh3) < eps) || (abs(mass - suf_m_h2o_nh3) < eps)) {
					annotated = true;
					break;
				}
			}

			if (annotated)
				v1.push_back(j);
		}

		vector<int> v2; // annotated positions in p1;
		for (int j = 0; j < clear_pept.length() - 1; j++) { // cut after j + 1 positions
			long double pref_m = pref_sums[j]; // pref, L = j + 1;
			long double suf_m = suf_sums[clear_pept.length() - 2 - j]; // suff, L = n - j - 1;
			long double pref_m_h2o = pref_m - h2o;
			long double pref_m_nh3 = pref_m - nh3;
			long double pref_m_h2o_nh3 = pref_m - h2o - nh3;

			long double suf_m_h2o = suf_m - h2o;
			long double suf_m_nh3 = suf_m - nh3;
			long double suf_m_h2o_nh3 = suf_m - h2o - nh3;

			bool annotated = false;
			for (long double mass : p2.actual_masses) {
				if ((abs(mass - pref_m) < eps) || (abs(mass - pref_m_h2o) < eps) || (abs(mass - pref_m_nh3) < eps) || (abs(mass - pref_m_h2o_nh3) < eps) ||
					(abs(mass - suf_m) < eps) || (abs(mass - suf_m_h2o) < eps) || (abs(mass - suf_m_nh3) < eps) || (abs(mass - suf_m_h2o_nh3) < eps)) {
					annotated = true;
					break;
				}
			}

			if (annotated)
				v2.push_back(j);
		}

		ofs << "Pair of scans for peptide " << pept << "\n";
		if (p1.frag_method == "HCD" && p2.frag_method == "CID") {
			int h = std::unique(v1.begin(), v1.end()) - v1.begin();
			int c = std::unique(v2.begin(), v2.end()) - v2.begin();
			ofs << "For HCD precantage of annotated cuts = " << 100.0 * h / (clear_pept.length() - 1) << "%\n";
			ofs << "For CID precantage of annotated cuts = " << 100.0 * c / (clear_pept.length() - 1) << "%\n";
			make_picture(pept, v2, v1, p2.scan_number, p1.scan_number);
		}
		else if (p2.frag_method == "HCD" && p1.frag_method == "CID") {
			int c = std::unique(v1.begin(), v1.end()) - v1.begin();
			int h = std::unique(v2.begin(), v2.end()) - v2.begin();
			ofs << "For HCD precantage of annotated cuts = " << 100.0 * h / (clear_pept.length() - 1) << "%\n";
			ofs << "For CID precantage of annotated cuts = " << 100.0 * c / (clear_pept.length() - 1) << "%\n";
			make_picture(pept, v1, v2, p1.scan_number, p2.scan_number);
		}
		ofs << endl;

		i++;
	}
	ofs.close();
}
int main() {
	freopen("piq.msalign", "r", stdin);
	freopen("annotations2.txt", "w", stdout);

	init();
	parse();
	fill_pept_by_scan();
	//resume_annotated();
	build_annotated_positions();

	return 0;
}
