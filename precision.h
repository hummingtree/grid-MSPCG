/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/ConjugateGradient.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef JKS_PRECISION_H
#define JKS_PRECISION_H

#include <Grid/Grid.h>

#include <qlat/qlat.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class double_type, class float_type>
void D2F(const double_type& in, float_type& out){
	const double* dp = (const double*)&in;
	float* fp = (float*)&out;
	for(int i = 0; i < sizeof(double_type)/sizeof(double); i++){
		fp[i] = dp[i];
	}
}

template<class float_type, class double_type>
void F2D(const float_type& in, double_type& out){
	double* dp = (double*)&out;
	const float* fp = (const float*)&in;
	for(int i = 0; i < sizeof(double_type)/sizeof(double); i++){
		dp[i] = fp[i];
	}
}

#endif
