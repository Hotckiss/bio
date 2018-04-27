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
	long double evalue;
	long double precursor;

};

const long double eps = 1e-10; // 10^-10, 10^-5

const int graph_len = 70;
map<string, int> col_id;
map<int, string> id_col;
vector<peptide> scans_eval, scans_number, all_scans;
vector<peptide> pair_hcd, pair_cid, no_pair;
map<string, vector<peptide>> groups_by_peptide;
int cid_exp[100];
int hcd_exp[100];
bool cmp_id(peptide &p1, peptide &p2) {
	return p1.scan_number < p2.scan_number;
}

int digit_number(long long num) {
	if (num == 0) {
		return 1;
	}
	int ans = 0;
	while (num > 0) {
		ans++;
		num /= 10;
	}
	return max(ans - 1, 0);
}

int exp_value(long double num) {
	int ans = 0;

	while (num < 1) {
		ans++;
		num *= 10;
	}

	return ans;
}

int d_exp_hcd_cid(long double hcd, long double cid) {
	return exp_value(hcd) - exp_value(cid);
}

void out2d(int num) {
	if (num < 10) {
		cout << "0" << num;
	}
	else cout << num;
}
int main() {
	freopen("db.tsv", "r", stdin);
	freopen("otchet.txt", "w", stdout);

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
	int maxlen = -1;
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

			if (id_col[id] == "FragMethod") {
				p.frag_method = pch;
			}

			if (id_col[id] == "EValue") {
				p.evalue = stod(pch);
			}

			if (id_col[id] == "Peptide") {
				p.peptid = pch;
			}

			if (id_col[id] == "Precursor") {
				p.precursor = stod(pch); // по Х порядок отличия по У колво пар с таким порядком, по тем где много пар 4-8. +
			}

			id++;
			pch = strtok(NULL, "\t");    // для частых по времени отличие порядка
		}                                // посмотреть какое было евалюе для тех пар, у которых отсеялось пара по теочности и посмотреть отличие точности


		if (p.evalue < eps) {
			scans_eval.push_back(p);
			scans_number.push_back(p);
		}
		all_scans.push_back(p);
	}

	for (peptide& pep : scans_eval) {
		int ll = pep.peptid.length();
		if (ll > maxlen) {
			maxlen = ll;
		}
		groups_by_peptide[pep.peptid].push_back(pep);
	}

	for (pair<string, vector<peptide>> vec : groups_by_peptide) {
		sort(vec.second.begin(), vec.second.end(), cmp_id);
	}

	sort(scans_number.begin(), scans_number.end(), cmp_id);
	sort(all_scans.begin(), all_scans.end(), cmp_id);

	int hcd_better = 0, cid_better = 0;
	double herror = 0, cerror = 0;
	int hnum = 0, cnum = 0;
	int tr = (int)scans_number.size() / 3;
	int n1 = 0, t1 = 0, n2 = 0, t2 = 0, n3 = 0, t3 = 0;
	for (int i = 0; i < scans_number.size() - 1; i++) {
		peptide& p1 = scans_number[i];
		peptide& p2 = scans_number[i + 1];

		if (p1.scan_number + 1 == p2.scan_number) {
			cout << "Парное измерение: #" << p1.scan_number << " #" << p2.scan_number << "\n";
			if (p1.frag_method == "HCD") {
				cout << "HCD Evalue = " << p1.evalue << "; CID EValue = " << p2.evalue << "; HCD Precursor = " << p1.precursor << "; CID Precursor = " << p2.precursor;
				cout << "\n";
				if (p1.evalue < p2.evalue) {
					int delta = d_exp_hcd_cid(p1.evalue, p2.evalue);
					cout << "HCD лучше, порядок отличия Evalue = " << delta; //digit_number((long long)(p2.evalue / p1.evalue));
																			 //hcd_exp[digit_number((long long)(p2.evalue / p1.evalue))]++;
					hcd_exp[delta]++;
				}
				else {
					int delta = d_exp_hcd_cid(p1.evalue, p2.evalue);
					cout << "CID лучше, порядок отличия Evalue = " << -delta;//digit_number((long long)(p1.evalue / p2.evalue));
																			 //cid_exp[digit_number((long long)(p1.evalue / p2.evalue))]++;
					cid_exp[-delta]++;
				}
				cout << "\n------------------\n\n";
				pair_hcd.push_back(p1);
				pair_cid.push_back(p2);
				if (p1.evalue < p2.evalue) {
					hcd_better++;
					herror += p1.evalue;
					hnum++;
					if (i < tr) {
						n1++, t1++;
					}
					else if (i < tr * 2) {
						n2++, t2++;
					}
					else {
						n3++, t3++;
					}
				}
				else {
					cid_better++;
					cerror += p2.evalue;
					cnum++;
					if (i < tr) {
						t1++;
					}
					else if (i < tr * 2) {
						t2++;
					}
					else {
						t3++;
					}
				}
			}
			else {
				cout << "HCD Evalue = " << p2.evalue << "; CID EValue = " << p1.evalue << "; HCD Precursor = " << p2.precursor << "; CID Precursor = " << p1.precursor;

				cout << "\n";

				if (p2.evalue < p1.evalue) {
					cout << "HCD лучше, порядок отличия Evalue = " << digit_number((long long)(p1.evalue / p2.evalue));
					hcd_exp[digit_number((long long)(p1.evalue / p2.evalue))]++;
				}
				else {
					cout << "CID лучше, порядок отличия Evalue = " << digit_number((long long)(p2.evalue / p1.evalue));
					cid_exp[digit_number((long long)(p2.evalue / p1.evalue))]++;
				}
				cout << "\n------------------\n\n";
				pair_hcd.push_back(p2);
				pair_cid.push_back(p1);
				if (p2.evalue < p1.evalue) {
					hcd_better++;
					herror += p2.evalue;
					hnum++;
					if (i < tr) {
						n1++, t1++;
					}
					else if (i < tr * 2) {
						n2++, t2++;
					}
					else {
						n3++, t3++;
					}
				}
				else {
					cid_better++;
					cerror += p1.evalue;
					cnum++;
					if (i < tr) {
						t1++;
					}
					else if (i < tr * 2) {
						t2++;
					}
					else {
						t3++;
					}
				}
			}
			i++;
		}
		else {
			no_pair.push_back(p1);
		}
	}

	cout << "Для " << pair_cid.size() << " пар измерений результаты сравнения погрешности оказались следующие:\n";

	int num_cid = (double)cid_better / (cid_better + hcd_better) * graph_len;
	int num_hcd = graph_len - num_cid;

	cout << "CID -- " << cid_better << " измерений из " << cid_better + hcd_better << " оказались точнее\n";
	cout << "HCD -- " << hcd_better << " измерений из " << cid_better + hcd_better << " оказались точнее.\n\n";

	cout << "Средняя ошибка CID по тем измерениям, где он лучше -- " << (cerror / cnum) << "\n";
	cout << "Средняя ошибка HCD по тем измерениям, где он лучше -- " << (herror / hnum) << "\n\n";

	cout << "CID ";
	for (int i = 0; i < num_cid; i++) {
		cout << "*";
	}
	cout << "\n";

	cout << "HCD ";
	for (int i = 0; i < num_hcd; i++) {
		cout << "*";
	}
	cout << "\n";

	int nohcd = 0;
	for (int i = 0; i < no_pair.size(); i++) {
		if (no_pair[i].frag_method == "HCD")
			nohcd++;
	}
	cout << "\nБез пары оказалось " << no_pair.size() << " измерений, " << nohcd << " из них это HCD, а " << no_pair.size() - nohcd << " это CID.\n--------------------------\n\n";

	cout << "Среди пептидов, для которых было 8 и более измерений:\n\n";

	int tot = 0, hcdb = 0;
	for (pair<string, vector<peptide>> vec : groups_by_peptide) {
		double erh = 0, erc = 0;
		int nc = 0, nh = 0;
		if (vec.second.size() >= 8) {
			for (peptide p1 : vec.second) {
				if (p1.frag_method == "HCD") {
					nh++;
					erh += p1.evalue;
				}
				else {
					nc++;
					erc += p1.evalue;
				}
			}
			string pept = vec.first;
			int l = pept.length();

			for (int k = l; k < maxlen; k++)
				pept += " ";

			cout << pept << " средняя ошибка HCD = " << erh / nh << ", CID = " << erc / nc;
			if (erh / nh < erc / nc) {
				cout << ", HCD точнее\n";
				hcdb++;
			}
			else {
				cout << ", CID точнее\n";
			}
			tot++;
		}


	}

	cout << "----------------------------\n\nCID был точнее " << tot - hcdb << " раз из " << tot << "\n\n";

	int num_cid1 = (double)(tot - hcdb) / tot * graph_len;
	cout << "CID ";
	for (int i = 0; i < num_cid1; i++) {
		cout << "*";
	}
	cout << "\n";

	cout << "HCD ";
	for (int i = 0; i < graph_len - num_cid1; i++) {
		cout << "*";
	}
	cout << "\n--------------------------\n\n";

	cout << "Среди 1й трети измерений по времени HCD был точнее " << n1 << " раз из " << t1 << "\n";
	cout << "Среди 2й трети измерений по времени HCD был точнее " << n2 << " раз из " << t2 << "\n";
	cout << "Среди 3й трети измерений по времени HCD был точнее " << n3 << " раз из " << t3 << "\n";

	cout << "\n";

	cout << "Распределение порядков, где лучше HCD:\n";
	for (int i = 0; i < 20; i++) {
		cout << hcd_exp[i] << " ";
	}
	cout << "\n\n";
	cout << "Распределение порядков, где лучше CID:\n";

	for (int i = 0; i < 20; i++) {
		cout << cid_exp[i] << " ";
	}

	cout << "\n--------------------------\n\n";
	cout << "Гистограммы, которые для пептидов с большим количеством встретившихся пар(хотя бы 6) показывают,\nв скольких парах Evalue было меньше у HCD, и на сколько порядков, аналогично CID.";
	cout << "\nПо горизонтали откладывается количество порядков, по вертикали - для такого порядка количества пар, где точнее HCD или CID.\n\n";
	for (pair<string, vector<peptide>> vec : groups_by_peptide) {
		double erh = 0, erc = 0;
		int nc = 0, nh = 0;
		if (vec.second.size() >= 12) {
			sort(vec.second.begin(), vec.second.end(), cmp_id);
			cout << "Пептид: " << vec.first << "\n";
			vector<peptide> & vp = vec.second;
			vector<int> exponent_cid(20, 0), exponent_hcd(20, 0);

			for (int i = 0; i < vp.size() - 1; i++) {
				peptide &p1 = vp[i];
				peptide &p2 = vp[i + 1];
				if (p1.scan_number + 1 == p2.scan_number) {
					if (p1.frag_method == "HCD") {
						if (p1.evalue < p2.evalue) {
							int d = d_exp_hcd_cid(p1.evalue, p2.evalue);
							exponent_hcd[d]++;
						}
						else {
							int d = d_exp_hcd_cid(p1.evalue, p2.evalue);
							exponent_cid[-d]++;
						}
					}
					else {
						if (p2.evalue < p1.evalue) {
							int d = d_exp_hcd_cid(p2.evalue, p1.evalue);
							exponent_hcd[d]++;
						}
						else {
							int d = d_exp_hcd_cid(p2.evalue, p1.evalue);
							exponent_cid[-d]++;
						}
					}
				}
			}

			char gist1[20][20]; // FIRST = EXP SECOND = NUMBER; gist for hcd;
			for (int i = 0; i < 20; i++)
				for (int j = 0; j < 20; j++)
					gist1[i][j] = ' ';

			for (int i = 0; i < 20; i++) {
				for (int j = 0; j < exponent_hcd[i]; j++)
					gist1[i][j] = 'H';
			}

			char gist2[20][20]; // FIRST = EXP SECOND = NUMBER; gist for hcd;
			for (int i = 0; i < 20; i++)
				for (int j = 0; j < 20; j++)
					gist2[i][j] = ' ';

			for (int i = 0; i < 20; i++) {
				for (int j = 0; j < exponent_cid[i]; j++)
					gist2[i][j] = 'C';
			}

			for (int i = 14; i < 20; i++) {
				for (int j = 0; j < 17; j++) {
					cout << gist2[j][19 - i] << gist1[j][19 - i] << " ";
					//cout << gist2[j][19 - i] << gist1[j][19 - i] << " ";
				}
				cout << "\n";
			}
			for (int i = 0; i < 17; i++) {
				out2d(i);
				cout << " ";
			}
			cout << "\n\n";
		}
	}

	cout << "\n\n";

	cout << "Теперь хочется для пептидов, у которых есть хотя бы 4 пары измерений,\nи если эти пары располагаются недалеко по времени, хотим посмотреть порядки отличия Evalue на этих парах.";
	cout << "\nСкажем, что пары расположены близко по времени, если между ними не более 300 измерений.\n\n";

	for (pair<string, vector<peptide>> vec : groups_by_peptide) {
		double erh = 0, erc = 0;
		int nc = 0, nh = 0;
		if (vec.second.size() >= 12) {
			//!! if deleted above
			sort(vec.second.begin(), vec.second.end(), cmp_id);
			vector<peptide> & vp = vec.second;
			vector<pair<peptide, peptide>> pairs;
			vector<int> times;
			for (int i = 0; i < vp.size() - 1; i++) {
				peptide &p1 = vp[i];
				peptide &p2 = vp[i + 1];

				if (p1.scan_number + 1 == p2.scan_number) {
					if (p1.frag_method == "HCD") {
						pairs.push_back({ p1, p2 });
					}
					else {
						pairs.push_back({ p2, p1 });
					}

					times.push_back(p2.scan_number);
				}
			}
			int mxx = -1;
			for (int i = 1; i < times.size(); i++)
				if (times[i] - times[i - 1] > mxx)
					mxx = times[i] - times[i - 1];

			if (mxx > 300)
				continue;
			cout << "Пептид: " << vec.first << "\n";
			cout << "Последовательность чисел, которые показывают на сколько порядков точнее(по Evalue) HCD чем CID, (отрицательно, если CID точнее)\n\n";

			for (pair<peptide, peptide> ppp : pairs) {
				cout << d_exp_hcd_cid(ppp.first.evalue, ppp.second.evalue) << " ";
			}
			cout << "\n\n";
		}
	}

	cout << "\n";

	int delc = 0, delh = 0;
	for (int i = 0; i < all_scans.size() - 1; i++) {
		peptide& p1 = all_scans[i];
		peptide& p2 = all_scans[i + 1];

		if (p1.scan_number + 1 == p2.scan_number) {
			if (p1.evalue < eps && p2.evalue >= eps) {
				cout << "Пара, где одно из измерений показало крайне низкую точность Evalue.\nОтличие измерений на " << d_exp_hcd_cid(p1.evalue, p2.evalue) << " порядков";
				//cerr << p2.evalue << endl;
				if (p1.frag_method == "HCD") {
					cout << ", отсеялось CID.\n\n";
					delc++;
				}
				else {
					cout << ", отсеялось HCD.\n\n";
					delh++;
				}
			}

			if (p1.evalue >= eps && p2.evalue < eps) {
				cout << "Пара, где одно из измерений показало крайне низкую точность Evalue.\nОтличие измерений на " << d_exp_hcd_cid(p2.evalue, p1.evalue) << " порядков";
				//cerr << p1.evalue << endl;
				if (p1.frag_method == "HCD") {
					cout << ", отсеялось HCD.\n\n";
					delh++;
				}
				else {
					cout << ", отсеялось CID.\n\n";
					delc++;
				}
			}
		}
	}

	cout << "Всего отсеялось по точности " << delc << " CID и " << delh << " HCD\n";
	return 0;
}
