#include <string>
#include <vector>
#include <iomanip>
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
	if (argc != 3)
	{
		cerr << "Please provide 2 coordinate files to match" << endl;
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
	cout << "\tReading structures" << endl;
	vector<structure> optKS1 = readallstruct(file1), optKS2 = readallstruct(file2), set1 = optKS1, set2 = optKS2;
    cout << endl << "\tFile " << file1 << ": " << optKS1.size() << endl;
    cout << "\tFile " << file2 << ": " << optKS2.size() << endl;
    
	if (optKS1.size() != optKS2.size())
	{
		cerr << "Warning: Number of structures differs." << endl;

        if (optKS2.size() > optKS1.size())
        {
            set1 = optKS2;
            set2 = optKS1;
        }
        cout << "\tSet 1: " << set1.size() << endl;
        cout << "\tSet 2: " << set2.size() << endl;
	}



	cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
	cout << "\tComparison of distance vectors" << endl;


	//compare function for vectors of doubles
	auto compare_vector_double = [&] (vector<double> a, vector<double> b) 
    {
        //cout << a.size() << " " << b.size() << endl;
		assert (a.size() == b.size());
		double eps = 1e-5;
		vector<double> diff;
		for (vector<double>::size_type i = 0; i < a.size(); i++) 
        {
			diff.push_back(abs(a[i]-b[i]));
		}
		auto max = max_element (begin(diff), end(diff));
		if (*max < eps) return true;
		return false;
	};


	vector< pair<structure, structure> > matchingKS;
	unsigned int fail(0);
	for (vector<structure>::size_type i = 0; i < set1.size(); i++)
	{
		bool matched(false);
		for (vector<structure>::size_type j = 0; j < set2.size(); j++)
		{
			matched = compare_vector_double (set1[i].getInterPartDist(), set2[j].getInterPartDist());	
			if (matched)
			{
				matchingKS.push_back(make_pair(set1[i], set2[j]));
				break;
			}
		}
		if (matched == false)
		{
			cout << "\t\tNo match found for structure " << set1[i].getNumber() << " in set1." << endl;
			fail++;
		}
	}

	if (fail > 0)
	{
		cout << "\tMismatches: " << fail << endl;
        cout << "\t\t" << setw(4) << "set1" << setw(4) << "set2" << endl;
        for (vector< pair<structure, structure> >::size_type i = 0; i < matchingKS.size(); i++)
        {
            cout << "\t\t" << setw(4) << matchingKS[i].first.getNumber() << setw(4) << matchingKS[i].second.getNumber() << endl;
        }
        cout << "\tMatches: " << matchingKS.size() << endl;
	} else
	{
		cout << "\tIdentical." << endl;
	}
    

	return 0;
}
