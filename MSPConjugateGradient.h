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
#ifndef GRID_MSP_CONJUGATE_GRADIENT_H
#define GRID_MSP_CONJUGATE_GRADIENT_H

namespace Grid {

template<class sc>
void expandLatticeFermion(const Lattice<sc>& in, Lattice<sc>& out){
	
	typedef typename sc::scalar_object sobj;
	
	GridBase* gbin = in._grid;
	GridBase* gbout = out._grid;
	assert(gbin._ndimension == gbout._ndimension);
	assert(gbin._ndimension == 5);
	for(int d = 1; d < 5; d++){ // Ls is the 0th index here.
		assert(gbin._ldimensions[d]+4 == gbout._ldimensions[d]);
	}

	for(int index = 0; index < out._grid->lSites(); index++){
		std::vector<int> lcoor_out, lcoor_in;
		out._grid->LocalIndexToLocalCoor(index, lcoor_out);
		
		lcoor_in[0] = lcoor_out[0];
		for(int d = 1; d < 5; d++){
			lcoor_in[d] = lcoor_out[d]-2;
		}
		
		if( lcoor_in[1]>=0 and lcoor_in[1]<gbin._ldimensions[1] and
			lcoor_in[2]>=0 and lcoor_in[2]<gbin._ldimensions[2] and
			lcoor_in[3]>=0 and lcoor_in[3]<gbin._ldimensions[3] and
			lcoor_in[4]>=0 and lcoor_in[4]<gbin._ldimensions[4]){

			sobj s;
			peekLocalSite(s, in, lcoor_in);
			pokeLocalSite(s, out, lcoor_out);
			
		}else{
			pokeLocalSite(0., out, lcoor_out);
		}
	}
	// Simply expand and then copy/merge.

}

template<class gf>
void expandLatticeGaugeField(const Lattice<gf>& in, Lattice<gf>& out){
	
	typedef typename gf::scalar_object sobj;
	
	GridBase* gbin = in._grid;
	GridBase* gbout = out._grid;
	assert(gbin._ndimension == gbout._ndimension);
	assert(gbin._ndimension == 5);
	for(int d = 1; d < 5; d++){ // Ls is the 0th index here.
		assert(gbin._ldimensions[d]+4 == gbout._ldimensions[d]);
	}

	for(int index = 0; index < out._grid->lSites(); index++){
		std::vector<int> lcoor_out, lcoor_in;
		out._grid->LocalIndexToLocalCoor(index, lcoor_out);
		
		lcoor_in[0] = lcoor_out[0];
		for(int d = 1; d < 5; d++){
			lcoor_in[d] = lcoor_out[d]-2;
		}
		
		if( lcoor_in[1]>=0 and lcoor_in[1]<gbin._ldimensions[1] and
			lcoor_in[2]>=0 and lcoor_in[2]<gbin._ldimensions[2] and
			lcoor_in[3]>=0 and lcoor_in[3]<gbin._ldimensions[3] and
			lcoor_in[4]>=0 and lcoor_in[4]<gbin._ldimensions[4]){

			sobj s;
			peekLocalSite(s, in, lcoor_in);
			pokeLocalSite(s, out, lcoor_out);
			
		}else{
			pokeLocalSite(0., out, lcoor_out);
		}
	}
	// Simply expand and then copy/merge.

}

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

class MSPConjugateGradient : public OperatorFunction<LatticeFermion> {
 public:
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  
  ConjugateGradient(RealD tol, Integer maxit, bool err_on_no_conv = true)
      : Tolerance(tol),
        MaxIterations(maxit),
        ErrorOnNoConverge(err_on_no_conv){};

  void operator()(CayleyFermion5D& M, const LatticeFermion& src,
                  LatticeFermion& psi) {
    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    RealD cp, c, a, d, b, ssq, qq, b_pred;

    Field p(src);
    Field mmp(src);
    Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    
    Linop.HermOpAndNorm(psi, mmp, d, b);
    

    r = src - mmp;
    p = r;

    a = norm2(p);
    cp = a;
    ssq = norm2(src);

    std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradient: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradient:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradient:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradient:   mmp " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradient:  cp,r " << cp << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradient:     p " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient: k=0 residual " << cp << " target " << rsq << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {
      c = cp;

      MatrixTimer.Start();
      Linop.HermOpAndNorm(p, mmp, d, qq);
      MatrixTimer.Stop();

      LinalgTimer.Start();
      //  RealD    qqck = norm2(mmp);
      //  ComplexD dck  = innerProduct(p,mmp);

      a = c / d;
      b_pred = a * (a * qq - d) / c;

      cp = axpy_norm(r, -a, mmp, r);
      b = cp / c;

      // Fuse these loops ; should be really easy
      psi = a * p + psi;
      p = p * b + r;

      LinalgTimer.Stop();

      std::cout << GridLogIterative << "ConjugateGradient: Iteration " << k
                << " residual " << cp << " target " << rsq << std::endl;
      std::cout << GridLogDebug << "a = "<< a << " b_pred = "<< b_pred << "  b = "<< b << std::endl;
      std::cout << GridLogDebug << "qq = "<< qq << " d = "<< d << "  c = "<< c << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();
        Linop.HermOpAndNorm(psi, mmp, d, qq);
        p = mmp - src;

        RealD srcnorm = sqrt(norm2(src));
        RealD resnorm = sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage << "ConjugateGradient Converged on iteration " << k << std::endl;
        std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq)<<std::endl;
		std::cout << GridLogMessage << "\tTrue residual " << true_residual<<std::endl;
		std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

        std::cout << GridLogMessage << "Time breakdown "<<std::endl;
		std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() <<std::endl;
		std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
		std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

		IterationsToComplete = k;	

        return;
      }
    }
    std::cout << GridLogMessage << "ConjugateGradient did NOT converge"
              << std::endl;

    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = k;

  }
};
}
#endif
