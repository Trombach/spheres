#include <string>
#include <vector>
#include "structure.h"
#include "iop.h"
#include "timer.h"
#include "lina.h"

using namespace std;

int main (int argc, char *argv[]) 
{
	timer T;

	cout << endl;

	//argument processing
	if (argc < 3)
	{
		cerr << "Please provide 2 coordinate files to compare" << endl;
		return 1;
	}


	string file1 = argv[1], file2 = argv[2];

	if (fexists(file1) && fexists(file2))
	{
		cout << "\tInput ok." << endl;
	}
	else
	{
		cerr << "Input incorrect." << endl;
		return 1;
	}

	cout << endl;

	//read structures
	vector<structure> optKS1 = readallstruct(file1), optKS2 = readallstruct(file2);

	if (optKS1.size() != optKS2.size())
	{
		cerr << "Number of structures differs." << endl;
		return 1;
	}

	cout << "\tReading structures" << endl;
#pragma omp parallel for	
	for (vector<structure>::size_type i = 0; i < optKS1.size(); i++)
	{


		//inter-particle distance
		vector<double> interPartDist1;
		vector<double> interPartDist2;
		vector<coord3d> currentCoord1 = optKS1[i].getCoordinates();
		vector<coord3d> currentCoord2 = optKS2[i].getCoordinates();
		for (vector<coord3d>::size_type j=0; j<currentCoord1.size(); j++) {
			for (vector<coord3d>::size_type k=j+1; k<currentCoord1.size(); k++) {
				interPartDist1.push_back(coord3d::dist(currentCoord1[j], currentCoord1[k]));	
				interPartDist2.push_back(coord3d::dist(currentCoord2[j], currentCoord2[k]));	
			}
		}
		sort(interPartDist1.begin(), interPartDist1.end());
		sort(interPartDist2.begin(), interPartDist2.end());
		optKS1[i].setInterPartDist(interPartDist1);
		optKS2[i].setInterPartDist(interPartDist2);
	}

	cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
	cout << "\tComparison of distance vectors" << endl;


	//compare function for vectors of doubles
	auto compare_vector_double = [&] (vector<double> a, vector<double> b) {
		assert (a.size() == b.size());
		double eps = 1e-5;
		vector<double> diff;
		for (vector<double>::size_type i = 0; i < a.size(); i++) {
			diff.push_back(abs(a[i]-b[i]));
		}
		auto max = max_element (begin(diff), end(diff));
		if (*max < eps) return true;
		return false;
	};


	vector< pair<structure, structure> > matchingKS;
	unsigned int fail(0);
	for (vector<structure>::size_type i = 1; i < optKS1.size(); i++)
	{
		bool matched(false);
		for (vector<structure>::size_type j = 0; j < optKS2.size(); j++)
		{
			matched = compare_vector_double (optKS1[i].getInterPartDist(), optKS2[j].getInterPartDist());	
			if (matched)
			{
				matchingKS.push_back(make_pair(optKS1[i], optKS2[j]));
				break;
			}
		}
		if (matched == false)
		{
			cout << "\t\tNo match found for structure " << optKS1[i].getNumber() << " in file " << file1 << endl;
			fail++;
		}
	}

	if (fail > 0)
	{
		cout << "\tMismatches: " << fail << endl;
	} else
	{
		cout << "\tIdentical." << endl;
	}

	return 0;
}
