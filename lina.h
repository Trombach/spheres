#ifndef LINA
#define LINA

#include <vector>
#include "geometry.h"

std::vector<double> diag (std::vector< std::vector<double> > &matrix);
std::vector<std::pair<double, std::vector<double> > > diagv (std::vector< std::vector<double> > &matrix);
matrix3d m3d_diagv (matrix3d &matrix);

#endif
