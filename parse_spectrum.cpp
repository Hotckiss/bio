#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

class peptide {
public:
	int scan_number;
	string frag_method;
	long double precursor_mass, precursor_charge;
	vector<long double> actual_masses;
};

const long double eps = 2e-2;
const long double val = 57.021;

vector<peptide> all_peptides;
map<int, string> pept_by_scan;
map<char, long double> amino_mass;

map<string, vector<long double>> pref_sums;
map<string, vector<long double>> suf_sums;
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

	return ans;
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
int main() {
	freopen("piq.msalign", "r", stdin);
	freopen("annotations.txt", "w", stdout);
	init();
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

	ifstream input ("table.txt", ifstream::in);

	int n;
	input >> n;

	for (int i = 0; i < n; i++) {
		int sn;
		string pept;
		
		input >> sn >> pept;

		pept_by_scan[sn] = pept;
	}

	input.close();

	cout << "SELECTED DETECTION ERROR = " << eps << "\n\n";
	for (peptide &p : all_peptides) {
		string pept(pept_by_scan[p.scan_number]);
		if (pept.empty())
			continue;
		vector<long double> pref_sums = get_prefix_sums(pept);
		vector<long double> suf_sums = get_suffix_sums(pept);
		vector<int> detected_pref;
		vector<int> detected_suf;
		for (long double mass : p.actual_masses) {
			for (int i = 0; i < pref_sums.size(); i++)
				if (abs(mass - pref_sums[i]) < eps)
					detected_pref.push_back(i);

			for (int i = 0; i < suf_sums.size(); i++)
				if (abs(mass - suf_sums[i]) < eps)
					detected_suf.push_back(i);
		}
		if (detected_pref.size() + detected_suf.size() > 0) {
			cout << "Detected for peptide " << pept << ":\n" << "SCAN = " << p.scan_number << "; FRAG METHOD = " << p.frag_method << "; PRECURSOR MASS = " << p.precursor_mass << "; CHARGE = " << p.precursor_charge << "\n\n";

			size_t pos = pept.find("+57.021");
			if (pos != string::npos) {
				pept = pept.substr(0, pos) + pept.substr(pos + 7, pept.length());
			}

			for (int i = 0; i < detected_pref.size(); i++) {
				cout << "Prefix " << pept.substr(0, detected_pref[i] + 1) << "; MASS = " << pref_sums[detected_pref[i]] << "\n";
			}

			for (int i = 0; i < detected_suf.size(); i++) {
				cout << "Suffix " << pept.substr(pept.length() - 1 - detected_suf[i], pept.length()) << "; MASS = " << suf_sums[detected_suf[i]] << "\n";
			}

			cout << "------------------------\n\n\n";
		}
	}

	return 0;
}
