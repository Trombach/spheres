#ifndef GLOBALS
#define GLOBALS

#include <dlib/matrix.h>

/*
Global definition of column_vector (dlib library).

Precompiler command that allows std::vector to be printed without for-loop:
cout << vector << endl;
*/

typedef dlib::matrix<double,0,1> column_vector;

#define container_output(container) \
template <typename T> std::ostream& operator<<(std::ostream& s, const container<T>& v) \
    { \
    s << "{"; \
    for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
        s << *x; \
        if(++x!=v.end()) s << ","; \
    } \
    s << "}"; \
    return s; \
    }
container_output(std::vector);

#endif
