#include <string>
#include <vector>
#include <iomanip>
#include <map>
#include "structure.h"
#include "iop.h"
#include "timer.h"
#include "lina.h"

using namespace std;

#define container_output(container) \
template <typename T> ostream& operator<<(ostream& s, const container<T>& v) \
	{ \
	s << "{"; \
	for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
		s << *x; \
		if(++x!=v.end()) s << ","; \
	} \
	s << "}"; \
	return s; \
	}
container_output(vector);

struct matchingStructure
{
    matchingStructure() : matched(false) {}
    matchingStructure(structure x) : matched(false) {KS = x;}
    bool matched;
    structure KS;
};

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

	//read structures and sort sets according to size
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

    //put sets in data structures
    vector<matchingStructure> set1_match, set2_match;
    for (vector<structure>::size_type i = 0; i < set1.size(); i++)
    {
        matchingStructure KS(set1[i]);
        set1_match.push_back(KS);
    }

    for (vector<structure>::size_type i = 0; i < set2.size(); i++)
    {
        matchingStructure KS(set2[i]);
        set2_match.push_back(KS);
    }

	cout << "\t\tTiming: " << T.timing() << " s" << endl << endl;
	cout << "\tComparison of distance vectors" << endl;


	//compare function for vectors of doubles
	auto compare_vector_double = [&] (vector<double> a, vector<double> b) 
    {
		assert (a.size() == b.size());
		double eps = 1e-4;
		vector<double> diff;
		for (vector<double>::size_type i = 0; i < a.size(); i++) 
        {
			diff.push_back(abs(a[i]-b[i]));
		}
		auto max = max_element (begin(diff), end(diff));
		if (*max < eps) return true;
		return false;
	};


    //matching
	vector< pair<matchingStructure, matchingStructure> > matchingKS;
	unsigned int fail(0);
	for (vector<matchingStructure>::size_type i = 0; i < set1_match.size(); i++)
	{
		bool matched(false);
		for (vector<matchingStructure>::size_type j = 0; j < set2_match.size(); j++)
		{
            if (set2_match[j].matched) continue;
			matched = compare_vector_double (set1_match[i].KS.getInterPartDist(), set2_match[j].KS.getInterPartDist());	
			if (matched)
			{
                set1_match[i].matched = true;
                set2_match[j].matched = true;
				matchingKS.push_back(make_pair(set1_match[i], set2_match[j]));
				break;
			}
		}
		if (matched == false)
		{
			//cout << "\t\tNo match found for structure " << set1[i].getNumber() << " in set1." << endl;
			fail++;
		}
	}

	if (fail > 0)
	{
		cout << "\tUnmatched: " << fail << endl;
        cout << "\t\tset1:" << endl;
        for (vector<matchingStructure>::size_type i = 0; i < set1_match.size(); i++)
        {
            if (!set1_match[i].matched) cout << "\t\t" << set1_match[i].KS.getNumber() << " " << set1_match[i].KS.getMomentOfInertia() << endl;
        }

        cout << "\t\tset2:" << endl;
        for (vector<matchingStructure>::size_type i = 0; i < set2_match.size(); i++)
        {
            if (!set2_match[i].matched) cout << "\t\t" << set2_match[i].KS.getNumber() << " " << set2_match[i].KS.getMomentOfInertia() << endl;
        }

        cout << "\tMatched: " << matchingKS.size() << endl;
        cout << "\t\t" << setw(4) << "set1" << setw(4) << "set2" << endl;
        for (vector< pair<matchingStructure, matchingStructure> >::size_type i = 0; i < matchingKS.size(); i++)
        {
            cout << "\t\t" << setw(4) << matchingKS[i].first.KS.getNumber() << setw(4) << matchingKS[i].second.KS.getNumber() << endl;
        }
	} else
	{
		cout << "\tIdentical." << endl;
	}
    

	return 0;
}
