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

	bool operator<(const peptide& o) const {
		return scan_number < o.scan_number;
	}

};

map<string, int> col_id;
map<int, string> id_col;
vector<peptide> all_scans;

int main() {
	freopen("db.tsv", "r", stdin);
	freopen("table.txt", "w", stdout);
	string str;

	getline(cin, str);

	char * pch = strtok(const_cast<char*>(str.c_str()), "\t");

	int id = 0;
	while (pch != NULL) {
		col_id[pch] = id;
		id_col[id] = pch;
		id++;
		pch = strtok(NULL, "\t");
	}

	while (true) {
		getline(cin, str);
		if (str.empty())
			break;

		pch = strtok(const_cast<char*>(str.c_str()), "\t");

		id = 0;

		peptide p;
		while (pch != NULL) {

			if (id_col[id] == "ScanNum") {
				p.scan_number = atoi(pch);
			}

			if (id_col[id] == "Peptide") {
				p.peptid = pch;
			}

			id++;
			pch = strtok(NULL, "\t");
		}

		all_scans.push_back(p);
	}
	sort(all_scans.begin(), all_scans.end());

	for (peptide &p : all_scans) {
		cout << p.scan_number << " " << p.peptid << "\n";
	}
	return 0;
}
