#include <iostream>  
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
	string peptid;
	long double evalue, precursor_mass, precursor_charge;
	vector<long double> masses;
};

const long double eps = 1e-2;

vector<peptide> all_peptides;
int main() {
	freopen("piq.msalign", "r", stdin);

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

				p.masses.push_back(stod(pch));
			}
			getline(cin, str, '\n');

			all_peptides.push_back(p);
		}
		else {
			break;
		}
	}

	cerr << all_peptides.size() << "\n";
	return 0;
}
