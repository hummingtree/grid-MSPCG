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

#include <Grid/Grid.h>

#include <qlat/qlat.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

class ExpandGrid{
public:
	GridCartesian* coarse_gauge_grid;
	GridCartesian* fine_gauge_grid;
	GridCartesian* coarse_fermion_grid;
	GridCartesian* fine_fermion_grid; //
	GridRedBlackCartesian* coarse_gauge_rb_grid; //
	GridRedBlackCartesian* fine_gauge_rb_grid; //
	GridRedBlackCartesian* coarse_fermion_rb_grid;
	GridRedBlackCartesian* fine_fermion_rb_grid;
	int expansion;

	std::vector<std::vector<int>> fine_fermion_rb_grid_in_list;
	std::vector<std::vector<int>> coarse_fermion_rb_grid_in_list;
	std::vector<std::vector<int>> fine_fermion_rb_grid_out_list;

	std::vector<std::vector<int>> fine_gauge_grid_in_list;
	std::vector<std::vector<int>> coarse_gauge_grid_in_list;

	inline void init(GridCartesian* _coarse_gauge_grid, int _expansion, int Ls){

		expansion = _expansion;		
		coarse_gauge_grid = _coarse_gauge_grid;
		coarse_gauge_rb_grid = SpaceTimeGrid::makeFourDimRedBlackGrid(coarse_gauge_grid);
		coarse_fermion_grid = SpaceTimeGrid::makeFiveDimGrid(Ls, coarse_gauge_grid);
		coarse_fermion_rb_grid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, coarse_gauge_grid);

		std::vector<int> simd_layout = coarse_gauge_grid->_simd_layout;
    	std::vector<int> mpi_layout  = coarse_gauge_grid->_processors;
    	std::vector<int> latt_size   = coarse_gauge_grid->_fdimensions;
    	
		std::vector<int> expanded_latt_size = latt_size;
    	for(int i = 0; i < 4; i++){
       		expanded_latt_size[i] = latt_size[i] + 2*expansion*mpi_layout[i];
    	}

		fine_gauge_grid = SpaceTimeGrid::makeFourDimGrid(expanded_latt_size, simd_layout, mpi_layout);
		fine_gauge_rb_grid = SpaceTimeGrid::makeFourDimRedBlackGrid(fine_gauge_grid);
		fine_fermion_grid = SpaceTimeGrid::makeFiveDimGrid(Ls, fine_gauge_grid);
		fine_fermion_rb_grid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, fine_gauge_grid);

		fine_fermion_rb_grid_in_list.resize(0);
		fine_fermion_rb_grid_out_list.resize(0);

// fermion copy list. Maybe by directly get index for memory we will get a better performance?
		for(int index = 0; index < fine_fermion_rb_grid->lSites(); index++){
			std::vector<int> lcoor_out(5), lcoor_in(5), full_lcoor_out(5), full_lcoor_in(5);
			fine_fermion_rb_grid->LocalIndexToLocalCoor(index, lcoor_out);

			lcoor_in[0] = lcoor_out[0]; // s: Ls
			lcoor_in[1] = lcoor_out[1] - expansion/2; // x is the _checker_dim
			for(int d = 2; d < 5; d++){
				lcoor_in[d] = lcoor_out[d] - expansion;
			}

			for(int cb = 0; cb < 2; cb++){

				for(int d = 0; d < 5; d++){ // Ls is the 0th index here.
					full_lcoor_out[d] = lcoor_out[d];
					full_lcoor_in[d] = lcoor_in[d];
				}

				full_lcoor_out[1] = lcoor_out[1]*2+cb;
				full_lcoor_in[1] = lcoor_in[1]*2+cb;

				if(1 != fine_fermion_rb_grid->CheckerBoard(full_lcoor_out)) continue;


				if( lcoor_in[1] >= 0 and lcoor_in[1] < coarse_fermion_rb_grid->_ldimensions[1] and
					lcoor_in[2] >= 0 and lcoor_in[2] < coarse_fermion_rb_grid->_ldimensions[2] and
					lcoor_in[3] >= 0 and lcoor_in[3] < coarse_fermion_rb_grid->_ldimensions[3] and
					lcoor_in[4] >= 0 and lcoor_in[4] < coarse_fermion_rb_grid->_ldimensions[4]){
					
					fine_fermion_rb_grid_in_list.push_back(full_lcoor_out);
					coarse_fermion_rb_grid_in_list.push_back(full_lcoor_in);

				}else{
					fine_fermion_rb_grid_out_list.push_back(full_lcoor_out);
				}

			}
		}

// gauge copy list.
		fine_gauge_grid_in_list.resize(0);
		coarse_gauge_grid_in_list.resize(0);
		for(int index = 0; index < fine_gauge_grid->lSites(); index++){
			std::vector<int> lcoor_out(4), lcoor_in(4);
			fine_gauge_grid->LocalIndexToLocalCoor(index, lcoor_out);

			for(int d = 0; d < 4; d++){
				lcoor_in[d] = lcoor_out[d] - expansion;
			}

			if( lcoor_in[0] >= 0 and lcoor_in[0] < coarse_gauge_grid->_ldimensions[0] and
				lcoor_in[1] >= 0 and lcoor_in[1] < coarse_gauge_grid->_ldimensions[1] and
				lcoor_in[2] >= 0 and lcoor_in[2] < coarse_gauge_grid->_ldimensions[2] and
				lcoor_in[3] >= 0 and lcoor_in[3] < coarse_gauge_grid->_ldimensions[3]){
				// local gauge links
				fine_gauge_grid_in_list.push_back(lcoor_out);
				coarse_gauge_grid_in_list.push_back(lcoor_in);
			}
		}
	}
};

// Double inner product
template<class vobj>
inline ComplexD local_inner_product(const Lattice<vobj> &left,const Lattice<vobj> &right) 
{
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_typeD vector_type;
  scalar_type  nrm;
  
  GridBase *grid = left._grid;
  
  std::vector<vector_type,alignedAllocator<vector_type> > sumarray(grid->SumArraySize());
  
  parallel_for(int thr=0;thr<grid->SumArraySize();thr++){
    int nwork, mywork, myoff;
    GridThread::GetWork(left._grid->oSites(),thr,mywork,myoff);
    
    decltype(innerProductD(left._odata[0],right._odata[0])) vnrm=zero; // private to thread; sub summation
    for(int ss=myoff;ss<mywork+myoff; ss++){
      vnrm = vnrm + innerProductD(left._odata[ss],right._odata[ss]);
    }
    sumarray[thr]=TensorRemove(vnrm) ;
  }
  
  vector_type vvnrm; vvnrm=zero;  // sum across threads
  for(int i=0;i<grid->SumArraySize();i++){
    vvnrm = vvnrm+sumarray[i];
  } 
  nrm = Reduce(vvnrm);// sum across simd
  return nrm;
}

template<class vobj> 
inline RealD local_norm_sqr(const Lattice<vobj> &arg){
  ComplexD nrm = local_inner_product(arg,arg);
  return std::real(nrm); 
}

template<class sc>
void zero_boundary_fermion(ExpandGrid& eg, Lattice<sc>& in){
	
	typedef typename sc::scalar_object sobj;

	parallel_for(size_t index = 0; index < eg.fine_fermion_rb_grid_out_list.size(); index++){
		sobj s; memset(&s, 0, sizeof(sobj));
		pokeLocalSite(s, in, eg.fine_fermion_rb_grid_out_list[index]);
	}
}

template<class F>
int local_conjugate_gradient_MdagM(ExpandGrid& eg, LinearOperatorBase<F> &D, const F &src, F &sol, int iterations){

	TIMER("local_conjugate_gradient_MdagM");
	sol.checkerboard = src.checkerboard;
//	conformable(src, sol);    

	sol = zero;

	F mp(src);
	F mmp(src);
	F r(src);
	F p(src);

	r = src;
	p = src;

	RealD rk2 = local_norm_sqr(r);
	RealD Mpk2, MdagMpk2, alpha, beta, rkp12;
	
	WilsonFermion5DStatic::dirichlet = true;

	for(int local_loop_count = 0; local_loop_count < iterations; local_loop_count++){

		D.Op(p, mp);
        D.AdjOp(mp, mmp);

		zero_boundary_fermion(eg, mmp);
		zero_boundary_fermion(eg, p); // not necessary?
        Mpk2 = std::real(local_inner_product(p, mmp));
        MdagMpk2 = local_norm_sqr(mmp); // Dag yes, (Mdag * M * p_k, Mdag * M * p_k)

        alpha = rk2 / Mpk2; 

		sol = sol + alpha * p;
		r = r - alpha * mmp;
		zero_boundary_fermion(eg, r);
		rkp12 = local_norm_sqr(r);

		beta = rkp12 / rk2;
		rk2 = rkp12;

		p = r + beta * p;

        if ( true ){
			zero_boundary_fermion(eg, sol);
            RealD sol2 = local_norm_sqr(sol);
            if( src._grid->ThisRank() == 0 ){
                printf("local_conjugate_gradient_MdagM: local iter. #%04d on NODE #%04d: r2 = %.4E psi2 = %.4E alpha = %.4E beta = %.4E \n",
                        local_loop_count, src._grid->ThisRank(), rk2, sol2, alpha, beta);
            }
        }
 		
	}
	
	WilsonFermion5DStatic::dirichlet = false;
    
	return 0;	
}

template<class F>
int MSP_conjugate_gradient(ExpandGrid& eg, LinearOperatorBase<F>& cD, LinearOperatorBase<F>& fD, const F& c_src, F& c_sol, RealD percent, int f_iter)
{

	RealD alpha, beta, rkzk, pkApk, zkp1rkp1, zkrk;

	Fermion_t p   = threadedAllocFermion(mem_fast); 
	Fermion_t tmp = threadedAllocFermion(mem_fast); 
	Fermion_t mp  = threadedAllocFermion(mem_fast); 
	Fermion_t mmp = threadedAllocFermion(mem_fast); 
	Fermion_t r   = threadedAllocFermion(mem_fast); 
	
	Fermion_t z   = threadedAllocFermion(mem_fast); 

	F p(c_src);
	F mp(c_src);
	F mmp(c_src);
	F r(c_src);

	cD.Op(c_sol, mp);
	cD.AdjOp(mp, mmp);

	r = src - mmp;
	

	d=Mprec(psi, mp, tmp, dag);
	b=Mprec(mp, mmp, tmp, ndag);

	axpy(r, mmp, src, -1.0);
	MSinv(z, r, dag); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
	copy(p, z); // p0 = z0
	// axpy(p, mmp, src, -1.0);

//	a = norm(p);
	cp = norm(r);

	Float ssq =  norm(src);
	ThreadBossDebug("CGNE_prec gues %le \n",guess);
	ThreadBossDebug("CGNE_prec src  %le \n",ssq);
	ThreadBossDebug("CGNE_prec  mp  %le \n",d);
	ThreadBossDebug("CGNE_prec  mmp %le \n",b);
	ThreadBossDebug("CGNE_prec   r  %le \n",cp);

	Float rsq =  residual * residual * ssq;

	//Check if guess is really REALLY good :)
	if ( cp <= rsq ) {
		ThreadBossMessage("CGNE_prec_MSP k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
		threadedFreeFermion(tmp);
		threadedFreeFermion(p);
		threadedFreeFermion(mp);
		threadedFreeFermion(mmp);
		threadedFreeFermion(r);
		threadedFreeFermion(z);
		if ( this->isBoss() && (!me) ) { 
			this->InverterExit();
		}
		return 0;
	}

	ThreadBossMessage("CGNE_prec_MSP: k=0 residual %le rsq %le time=0.0\n",cp,rsq);

	if ( this->watchfile ) {
		ThreadBossDebug("CGNE_prec watching file \"%s\"\n",this->watchfile);
	}
	struct timeval start,stop;
	gettimeofday(&start,NULL);

	int k;
	uint64_t t0,t1;
	t0 = GetTimeBase();
	for (k=1;k<=max_iter;k++){

		this->iter=k;
		uint64_t t_iter_1=GetTimeBase();

		// c=cp;
		rkzk = dot(r, z).real();

		uint64_t t_mprec_1=GetTimeBase();
		pkApk = Mprec(p, mp, tmp, dag, 1);// Dag no
		uint64_t t_mprec_2=GetTimeBase();
		alpha = rkzk / pkApk; // alpha_k

		uint64_t t_mprec_3=GetTimeBase();
		double qq=Mprec(mp, mmp, tmp, ndag); // Dag yes
		uint64_t t_mprec_4=GetTimeBase();
		
		axpy(psi, p, psi, alpha); // x_k+1 = x_k + alpha * p_k
		axpy(r, mmp, r, -alpha); // r_k+1 = r_k - alpha * Ap_k

		MSinv(z, r, dag); // z_k+1 = M^-1 * r_k+1

		zkp1rkp1 = dot(z, r).real();	
		
		beta = zkp1rkp1 / rkzk;
		axpy(p, p, z, beta); // p_k+1 = z_k+1 + beta * p_k
	
		cp = norm(r);

		uint64_t tpsi2=GetTimeBase();
		uint64_t t_iter_2=GetTimeBase();

		if ( ((k%1 == 0) && (verbose!=0)) || (verbose > 10) ){
			t1 = GetTimeBase();
			ThreadBossMessage("CGNE_prec_MSP/iter.count/r**2/time/: %05d %8.4e %8.4e %4.2f\n",k,cp,sqrt(cp/ssq),(t1-t0)/MHz()*1.0e-06);
		}

		// Stopping condition
		if ( cp <= rsq ) { 

			gettimeofday(&stop,NULL);
			struct timeval diff;
			timersub(&stop,&start,&diff);

			ThreadBossMessage("CGNE_prec converged in %d iterations\n",k);
			ThreadBossMessage("CGNE_prec converged in %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);


			double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
			flops = flops * k;

			double t = diff.tv_sec*1.0E6 + diff.tv_usec;
			ThreadBossMessage("CGNE_prec: %d mprec flops/site\n",mprecFlopsPerSite());
			ThreadBossMessage("CGNE_prec: %le flops\n",flops);
			ThreadBossMessage("CGNE_prec: %le mflops per node\n",flops/t);

			Mprec(psi,mp,tmp,dag);
			Mprec(mp,mmp,tmp,ndag); 
			axpy(tmp,src,mmp,-1.0);

			double  mpnorm = sqrt(norm(mp));
			double mmpnorm = sqrt(norm(mmp));
			double psinorm = sqrt(norm(psi));
			double srcnorm = sqrt(norm(src));
			double tmpnorm = sqrt(norm(tmp));
			double true_residual = tmpnorm/srcnorm;
			ThreadBossMessage("CGNE_prec: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
			ThreadBossMessage("CGNE_prec: target residual was %le \n",residual);
			ThreadBossMessage("CGNE_prec: mp %le, mmp %le\n",mpnorm,mmpnorm);

			threadedFreeFermion(tmp);
			threadedFreeFermion(p);
			threadedFreeFermion(mp);
			threadedFreeFermion(mmp);
			threadedFreeFermion(r);
			if ( this->isBoss() && (!me) ) { 
				this->InverterExit();
			}
			return k;

		}

	}
	ThreadBossMessage("CGNE_prec: CG not converged after %d iter resid %le\n",k,sqrt(cp/ssq));
	threadedFreeFermion(tmp);
	threadedFreeFermion(p);
	threadedFreeFermion(mp);
	threadedFreeFermion(mmp);
	threadedFreeFermion(r);
	threadedFreeFermion(z);
	if ( this->isBoss() && (!me) ) { 
		this->InverterExit();
	}

	return -1;

}

#endif
