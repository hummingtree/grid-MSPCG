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

	std::vector<size_t> fine_fermion_rb_grid_in_olist;
	std::vector<size_t> fine_fermion_rb_grid_in_ilist;
	
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
		fine_fermion_rb_grid_in_olist.resize(0);
		fine_fermion_rb_grid_in_ilist.resize(0);
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
			
					size_t odx, idx;
					odx = ((GridCartesian*)fine_fermion_rb_grid)->oIndex(full_lcoor_out);
					idx = ((GridCartesian*)fine_fermion_rb_grid)->iIndex(full_lcoor_out);
					// i/o lists
					fine_fermion_rb_grid_in_olist.push_back(odx);
					fine_fermion_rb_grid_in_ilist.push_back(idx);

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

template<class sc>
void expand_fermion(ExpandGrid& eg, const Lattice<sc>& in, Lattice<sc>& out){

	TIMER("expand_fermion");

//	conformable(eg.coarse_fermion_rb_grid, in._grid);
//	conformable(eg.fine_fermion_rb_grid, out._grid);

	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename sc::scalar_object sobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.fine_fermion_rb_grid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.coarse_fermion_rb_grid_in_list[index]);
        pokeLocalSite(s, out, eg.fine_fermion_rb_grid_in_list[index]);
    }

	parallel_for(size_t index = 0; index < eg.fine_fermion_rb_grid_out_list.size(); index++){
		sobj s; memset(&s, 0, sizeof(sobj));
		pokeLocalSite(s, out, eg.fine_fermion_rb_grid_out_list[index]);
	}

}

template<class sc>
void expand_fermion_qlat(ExpandGrid& eg, const Lattice<sc>& in, Lattice<sc>& out){
	
	TIMER("expand_fermion_qlat");

//	conformable(eg.coarse_fermion_rb_grid, in._grid);
//	conformable(eg.fine_fermion_rb_grid, out._grid);

	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename sc::scalar_object sobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.fine_fermion_rb_grid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.coarse_fermion_rb_grid_in_list[index]);
        pokeLocalSite(s, out, eg.fine_fermion_rb_grid_in_list[index]);
    }

	parallel_for(size_t index = 0; index < eg.fine_fermion_rb_grid_out_list.size(); index++){
		sobj s; memset(&s, 1, sizeof(sobj));
		pokeLocalSite(s, out, eg.fine_fermion_rb_grid_out_list[index]);
	}

	std::vector<sobj> out_lex(out._grid->lSites());
  	unvectorizeToLexOrdArray(out_lex, out);

	qlat::Coordinate global_size(in._grid->_gdimensions[1], in._grid->_gdimensions[2], in._grid->_gdimensions[3], in._grid->_gdimensions[4]);
	qlat::Geometry geo; geo.init(global_size, in._grid->_gdimensions[0]); // multiplicity = in._grid->_gdimensions[0]

	qlat::Coordinate expanse(eg.expansion/2, eg.expansion, eg.expansion, eg.expansion); // x direction is the checkerboard dimension so the /2.
	geo.resize(expanse, expanse);

	qlat::Field<sobj> f; f.init(geo);
	assert(f.field.size() == out._grid->lSites()); 
	memcpy(f.field.data(), out_lex.data(), f.field.size()*sizeof(sobj));
	// DO comm.
	refresh_expanded(f);
	
	memcpy(out_lex.data(), f.field.data(), f.field.size()*sizeof(sobj));
	vectorizeFromLexOrdArray(out_lex, out);
}

template<class sc>
void shrink_fermion(ExpandGrid& eg, const Lattice<sc>& in, Lattice<sc>& out){

	TIMER("shrink_fermion");

//	conformable(eg.coarse_fermion_rb_grid, in._grid);
//	conformable(eg.fine_fermion_rb_grid, out._grid);

	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename sc::scalar_object sobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.coarse_fermion_rb_grid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.fine_fermion_rb_grid_in_list[index]);
        pokeLocalSite(s, out, eg.coarse_fermion_rb_grid_in_list[index]);
    }

}

template<class vobj>
void expand_gauge_field_qlat(ExpandGrid& eg, const Lattice<vobj>& in, Lattice<vobj>& out){
	
	// Simply expand and then copy/merge.
	typedef typename vobj::scalar_object sobj;
	
	parallel_for(size_t index = 0; index < eg.fine_gauge_grid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.coarse_gauge_grid_in_list[index]);
        pokeLocalSite(s, out, eg.fine_gauge_grid_in_list[index]);
    }
	
	std::vector<sobj> out_lex(out._grid->lSites());
  	unvectorizeToLexOrdArray(out_lex, out);

	qlat::Coordinate global_size(in._grid->_gdimensions[0], in._grid->_gdimensions[1], in._grid->_gdimensions[2], in._grid->_gdimensions[3]);
	qlat::Geometry geo; geo.init(global_size, 1); // packing gauge links in four directions together, therefore multiplicity = 1

	qlat::Coordinate expanse(eg.expansion, eg.expansion, eg.expansion, eg.expansion);
	geo.resize(expanse, expanse);

	qlat::Field<sobj> f; f.init(geo);
	assert(f.field.size()*1 == out._grid->lSites()); // 1 for multiplicity
	memcpy(f.field.data(), out_lex.data(), f.field.size()*sizeof(sobj));
	// DO comm.
	refresh_expanded(f);
	
	memcpy(out_lex.data(), f.field.data(), f.field.size()*sizeof(sobj));
	vectorizeFromLexOrdArray(out_lex, out);

}

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

template<class vobj>
inline ComplexD local_inner_product_center(ExpandGrid& eg, const Lattice<vobj> &left, const Lattice<vobj> &right) 
{
	TIMER("local_inner_product_center");
	
	typedef typename vobj::scalar_type scalar_type;
	typedef typename vobj::vector_type vector_type;

	int Nsimd = left._grid->Nsimd();
	int num_threads;
#pragma omp parallel
{
	num_threads = omp_get_num_threads();
}	
	std::vector<ComplexD> ssum(num_threads, 0.);

#pragma omp parallel for
	for(size_t index = 0; index < eg.fine_fermion_rb_grid_in_olist.size(); index++){
		int thread_num = omp_get_thread_num();
		scalar_type* lp = (scalar_type*)&left._odata[eg.fine_fermion_rb_grid_in_olist[index]];
		scalar_type* rp = (scalar_type*)&right._odata[eg.fine_fermion_rb_grid_in_olist[index]];
		for(size_t w = 0; w < sizeof(vobj)/sizeof(vector_type); w++){
			size_t offset = eg.fine_fermion_rb_grid_in_ilist[index] + w * Nsimd;
	//		size_t offset = eg.fine_fermion_rb_grid_in_ilist[index];
			ssum[thread_num] += innerProduct(lp[offset], rp[offset]);
		}
	}

	ComplexD ret = 0.;
	for(int i = 0; i < num_threads; i++){
		ret += ssum[i];
	}

//	right._grid->GlobalSum(ret);

	return ret;
}

template<class vobj> 
inline RealD local_norm_sqr_center(ExpandGrid& eg, const Lattice<vobj> &arg){
  ComplexD nrm = local_inner_product_center(eg, arg, arg);
  return std::real(nrm); 
}

template<class sc>
void zero_boundary_fermion(ExpandGrid& eg, Lattice<sc>& in){

	TIMER("zero_boundary_fermion");

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

	static F mp(src);
	static F mmp(src);
	static F r(src);
	static F p(src);

	r = src;
	p = src;

	RealD rk2 = local_norm_sqr(r);
	RealD Mpk2, MdagMpk2, alpha, beta, rkp12;
	
	WilsonFermion5DStatic::dirichlet = true;

	for(int local_loop_count = 0; local_loop_count < iterations; local_loop_count++){
		{
		TIMER("local_conjugate_gradient/dslash");
		D.Op(p, mp);
        D.AdjOp(mp, mmp);
		}

		{
		TIMER("local_conjugate_gradient/linalg");
//		zero_boundary_fermion(eg, mmp);
//		zero_boundary_fermion(eg, p); // not necessary?
        Mpk2 = std::real(local_inner_product(p, mmp));
        MdagMpk2 = local_norm_sqr(mmp); // Dag yes, (Mdag * M * p_k, Mdag * M * p_k)

        alpha = rk2 / Mpk2; 

		sol = sol + alpha * p;
		r = r - alpha * mmp;
//		zero_boundary_fermion(eg, r);
		rkp12 = local_norm_sqr(r);

		beta = rkp12 / rk2;
		rk2 = rkp12;

		p = r + beta * p;
		}
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
int local_conjugate_gradient_MdagM_variant(ExpandGrid& eg, LinearOperatorBase<F> &D, const F &src, F &sol, int iterations){

	TIMER("local_CG_variant");
	sol.checkerboard = src.checkerboard;
//	conformable(src, sol);    

	sol = zero;

	static F mp(src);
	static F mmp(src);
	static F r(src);
	static F p(src);

	r = src;
	p = src;

//	RealD rk2 = local_norm_sqr_center(eg, r);
	RealD rk2 = local_norm_sqr(r);
	RealD Mpk2, MdagMpk2, alpha, beta, rkp12;
	
	WilsonFermion5DStatic::dirichlet = true;

	for(int local_loop_count = 0; local_loop_count < iterations; local_loop_count++){
		{
		TIMER("local_CG_variant/dslash");
		D.Op(p, mp);
        D.AdjOp(mp, mmp);
		}

		{
		TIMER("local_CG_variant/linalg");
//		zero_boundary_fermion(eg, mmp);
//		zero_boundary_fermion(eg, p); // not necessary?
//		Mpk2 = std::real(local_inner_product_center(eg, p, mmp));
        Mpk2 = std::real(local_inner_product(p, mmp));
//      MdagMpk2 = local_norm_sqr_center(eg, mmp); // Dag yes, (Mdag * M * p_k, Mdag * M * p_k)
        MdagMpk2 = local_norm_sqr(mmp); // Dag yes, (Mdag * M * p_k, Mdag * M * p_k)

        alpha = rk2 / Mpk2; 

		sol = sol + alpha * p;
		r = r - alpha * mmp;
//		zero_boundary_fermion(eg, r);
//		rkp12 = local_norm_sqr_center(eg, r);
		rkp12 = local_norm_sqr(r);

		beta = rkp12 / rk2;
		rk2 = rkp12;

		p = r + beta * p;
		}
        if ( true ){
			zero_boundary_fermion(eg, sol);
//			RealD sol2 = local_norm_sqr_center(eg, sol);
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
int MSP_conjugate_gradient(ExpandGrid& eg, LinearOperatorBase<F>& D, LinearOperatorBase<F>& fD, const F& src, F& sol, RealD percent, int f_iter, size_t max_iter=50000)
{

	RealD alpha, beta, rkzk, pkApk, zkp1rkp1, zkrk;
	RealD r2, src_sqr, tgt_r2;

	F p(src);
	F mp(src);
	F mmp(src);
	F r(src);
	F z(src);

	F f_r(eg.fine_fermion_rb_grid);
	F f_z(eg.fine_fermion_rb_grid);

	D.Op(sol, mp);
	D.AdjOp(mp, mmp);

	r = src - mmp;
	r2 = norm2(r);
	
	src_sqr = norm2(src);
	tgt_r2  = src_sqr * percent * percent;

	if(r2 < tgt_r2){
		if(not src._grid->ThisRank()){
			printf("MSP_conjugate_gradient() CONVERGED at iteration 0;\n");
		}
		return 0;
	}

	expand_fermion(eg, r, f_r);
	local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
	shrink_fermion(eg, f_z, z);

	p = z; // p0 = z0

	for(size_t k = 0; k < max_iter; k++){

		TIMER("MSPCG_iteration");

		rkzk = std::real(innerProduct(r, z));

		D.Op(p, mp);
		D.AdjOp(mp, mmp);
		pkApk = std::real(innerProduct(p, mmp));
		alpha = rkzk / pkApk; // alpha_k

		sol = sol + alpha * p; // x_k+1 = x_k + alpha * p_k
		r = r - alpha * mmp; // r_k+1 = r_k - alpha * Ap_k
    
//		expand_fermion_qlat(eg, r, f_r);
		expand_fermion(eg, r, f_r);
	    local_conjugate_gradient_MdagM_variant(eg, fD, f_r, f_z, f_iter); // z_k+1 = M^-1 * r_k+1
		shrink_fermion(eg, f_z, z);
		
		zkp1rkp1 = std::real(innerProduct(z, r));	
		
		beta = zkp1rkp1 / rkzk;
		p = z + beta * p; // p_k+1 = z_k+1 + beta * p_k
	
		r2 = norm2(r);

		if ( true ){
		    if(not src._grid->ThisRank()){
				printf("MSPCG/iter.count/r2/target r2/: %05d %8.4e %8.4e\n", k, r2, tgt_r2);
			}
		}

		// Stopping condition
		if ( r2 < tgt_r2 ) { 
		    
			D.Op(sol, mp);
		    D.AdjOp(mp, mmp);
		    r = src - mmp;
		    r2 = norm2(r);	
            
			if(not src._grid->ThisRank()){
                printf("MSPCG CONVERGED at iteration %05d with true r2 = %8.4e: target r2 = %8.4e\n", k, r2, tgt_r2);
            }
			
			return k;

		}

	}
    if(not src._grid->ThisRank()){
        printf("MSPCG has NOT CONVERGED after iteration %05d\n", max_iter);
    }
	
	return -1;
}

#endif
