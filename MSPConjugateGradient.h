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

#include "precision.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

class ExpandGrid{
public:
	GridCartesian* cUGrid;
	GridCartesian* fUGrid;
	GridCartesian* cFGrid;
	GridCartesian* fFGrid; //
	GridRedBlackCartesian* cUrbGrid; //
	GridRedBlackCartesian* fUrbGrid; //
	GridRedBlackCartesian* cFrbGrid;
	GridRedBlackCartesian* fFrbGrid;
	
	GridCartesian* fUGrid_F;
	GridCartesian* fFGrid_F;
	GridRedBlackCartesian* fUrbGrid_F;
	GridRedBlackCartesian* fFrbGrid_F;
	
	int expansion;

	std::vector<std::vector<int>> fFrbGrid_in_list;
	std::vector<std::vector<int>> cFrbGrid_in_list;
	std::vector<std::vector<int>> fFrbGrid_out_list;

	std::vector<size_t> fFrbGrid_in_olist;
	std::vector<size_t> fFrbGrid_in_ilist;
	
	std::vector<std::vector<int>> fUGrid_in_list;
	std::vector<std::vector<int>> cUGrid_in_list;

	inline void init(GridCartesian* _cUGrid, int _expansion, int Ls){

		expansion = _expansion;		
		cUGrid = _cUGrid;
		cUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(cUGrid);
		cFGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, cUGrid);
		cFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, cUGrid);

		std::vector<int> simd_layout = cUGrid->_simd_layout;
    	std::vector<int> mpi_layout  = cUGrid->_processors;
    	std::vector<int> mpi_layout_no_comm(4, 1); // no comm
    	std::vector<int> latt_size   = cUGrid->_fdimensions;
    	
		std::vector<int> expanded_latt_size = latt_size;
    	for(int i = 0; i < 4; i++){
			expanded_latt_size[i] = latt_size[i]/mpi_layout[i] + 2*expansion;
//   		expanded_latt_size[i] = latt_size[i] + 2*expansion*mpi_layout[i];
		}

		fUGrid = new GridCartesian(expanded_latt_size, simd_layout, mpi_layout_no_comm, *cUGrid);
//		fUGrid = SpaceTimeGrid::makeFourDimGrid(expanded_latt_size, simd_layout, mpi_layout);
		fUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(fUGrid);
		fFGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, fUGrid);
		fFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, fUGrid);
		
		fUGrid_F = new GridCartesian(expanded_latt_size, GridDefaultSimd(Nd,vComplexF::Nsimd()), mpi_layout_no_comm, *cUGrid);
//		fUGrid_F = SpaceTimeGrid::makeFourDimGrid(expanded_latt_size, GridDefaultSimd(Nd,vComplexF::Nsimd()), mpi_layout);
		fUrbGrid_F = SpaceTimeGrid::makeFourDimRedBlackGrid(fUGrid_F);
		fFGrid_F = SpaceTimeGrid::makeFiveDimGrid(Ls, fUGrid_F);
		fFrbGrid_F = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, fUGrid_F);

		fFrbGrid_in_list.resize(0);
		cFrbGrid_in_list.resize(0);
		fFrbGrid_in_olist.resize(0);
		fFrbGrid_in_ilist.resize(0);
		fFrbGrid_out_list.resize(0);

// fermion copy list. Maybe by directly get index for memory we will get a better performance?
		for(int index = 0; index < fFrbGrid->lSites(); index++){
			std::vector<int> lcoor_out(5), lcoor_in(5), full_lcoor_out(5), full_lcoor_in(5);
			fFrbGrid->LocalIndexToLocalCoor(index, lcoor_out);

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

				if(1 != fFrbGrid->CheckerBoard(full_lcoor_out)) continue;


				if( lcoor_in[1] >= 0 and lcoor_in[1] < cFrbGrid->_ldimensions[1] and
					lcoor_in[2] >= 0 and lcoor_in[2] < cFrbGrid->_ldimensions[2] and
					lcoor_in[3] >= 0 and lcoor_in[3] < cFrbGrid->_ldimensions[3] and
					lcoor_in[4] >= 0 and lcoor_in[4] < cFrbGrid->_ldimensions[4]){
					
					fFrbGrid_in_list.push_back(full_lcoor_out);
					cFrbGrid_in_list.push_back(full_lcoor_in);
			
					size_t odx, idx;
					odx = ((GridCartesian*)fFrbGrid)->oIndex(full_lcoor_out);
					idx = ((GridCartesian*)fFrbGrid)->iIndex(full_lcoor_out);
					// i/o lists
					fFrbGrid_in_olist.push_back(odx);
					fFrbGrid_in_ilist.push_back(idx);

				}else{
					fFrbGrid_out_list.push_back(full_lcoor_out);
				}

			}
		}

// gauge copy list.
		fUGrid_in_list.resize(0);
		cUGrid_in_list.resize(0);
		for(int index = 0; index < fUGrid->lSites(); index++){
			std::vector<int> lcoor_out(4), lcoor_in(4);
			fUGrid->LocalIndexToLocalCoor(index, lcoor_out);

			for(int d = 0; d < 4; d++){
				lcoor_in[d] = lcoor_out[d] - expansion;
			}

			if( lcoor_in[0] >= 0 and lcoor_in[0] < cUGrid->_ldimensions[0] and
				lcoor_in[1] >= 0 and lcoor_in[1] < cUGrid->_ldimensions[1] and
				lcoor_in[2] >= 0 and lcoor_in[2] < cUGrid->_ldimensions[2] and
				lcoor_in[3] >= 0 and lcoor_in[3] < cUGrid->_ldimensions[3]){
				// local gauge links
				fUGrid_in_list.push_back(lcoor_out);
				cUGrid_in_list.push_back(lcoor_in);
			}
		}
	}
};

template<class sc>
void expand_fermion(ExpandGrid& eg, const Lattice<sc>& in, Lattice<sc>& out){

	TIMER("expand_fermion");

//	conformable(eg.cFrbGrid, in._grid);
//	conformable(eg.fFrbGrid, out._grid);

	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename sc::scalar_object sobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.fFrbGrid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.cFrbGrid_in_list[index]);
		pokeLocalSite(s, out, eg.fFrbGrid_in_list[index]);
    }

	parallel_for(size_t index = 0; index < eg.fFrbGrid_out_list.size(); index++){
		sobj s; memset(&s, 0, sizeof(sobj));
		pokeLocalSite(s, out, eg.fFrbGrid_out_list[index]);
	}

}

template<class double_type, class float_type>
void expand_fermion_D2F(ExpandGrid& eg, const Lattice<double_type>& in, Lattice<float_type>& out){

	TIMER("expand_fermion_D2F");

//	conformable(eg.cFrbGrid, in._grid);
//	conformable(eg.fFrbGrid, out._grid);

	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename double_type::scalar_object Dsobj;
	typedef typename float_type::scalar_object Fsobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.fFrbGrid_in_list.size(); index++){
        Dsobj ds;
        Fsobj fs;
        peekLocalSite(ds, in, eg.cFrbGrid_in_list[index]);
		D2F(ds, fs);
		pokeLocalSite(fs, out, eg.fFrbGrid_in_list[index]);
    }

	parallel_for(size_t index = 0; index < eg.fFrbGrid_out_list.size(); index++){
		Fsobj fs; memset(&fs, 0, sizeof(Fsobj));
		pokeLocalSite(fs, out, eg.fFrbGrid_out_list[index]);
	}

}

template<class sc>
void expand_fermion_qlat(ExpandGrid& eg, const Lattice<sc>& in, Lattice<sc>& out){
	
	TIMER("expand_fermion_qlat");

//	conformable(eg.cFrbGrid, in._grid);
//	conformable(eg.fFrbGrid, out._grid);

	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename sc::scalar_object sobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.fFrbGrid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.cFrbGrid_in_list[index]);
        pokeLocalSite(s, out, eg.fFrbGrid_in_list[index]);
    }

	parallel_for(size_t index = 0; index < eg.fFrbGrid_out_list.size(); index++){
		sobj s; memset(&s, 1, sizeof(sobj));
		pokeLocalSite(s, out, eg.fFrbGrid_out_list[index]);
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

//	conformable(eg.cfrbgrid, in._grid);
//	conformable(eg.ffrbgrid, out._grid);

	// simply expand and then copy/merge.
	// set the boundary sites to zero.
	typedef typename sc::scalar_object sobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.cFrbGrid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.fFrbGrid_in_list[index]);
		pokeLocalSite(s, out, eg.cFrbGrid_in_list[index]);
    }

}

template<class float_type, class double_type>
void shrink_fermion_F2D(ExpandGrid& eg, const Lattice<float_type>& in, Lattice<double_type>& out){

	TIMER("shrink_fermion_F2D");

//	conformable(eg.cfrbgrid, in._grid);
//	conformable(eg.ffrbgrid, out._grid);

	// simply expand and then copy/merge.
	// set the boundary sites to zero.
	typedef typename float_type::scalar_object Fsobj;
	typedef typename double_type::scalar_object Dsobj;
	
	out.checkerboard = in.checkerboard;

	parallel_for(size_t index = 0; index < eg.cFrbGrid_in_list.size(); index++){
        Fsobj fs;
        Dsobj ds;
        peekLocalSite(fs, in, eg.fFrbGrid_in_list[index]);
		F2D(fs, ds);
		pokeLocalSite(ds, out, eg.cFrbGrid_in_list[index]);
    }

}

template<class vobj>
void expand_gauge_field_qlat(ExpandGrid& eg, const Lattice<vobj>& in, Lattice<vobj>& out){
	
	// Simply expand and then copy/merge.
	typedef typename vobj::scalar_object sobj;
	
	parallel_for(size_t index = 0; index < eg.fUGrid_in_list.size(); index++){
        sobj s;
        peekLocalSite(s, in, eg.cUGrid_in_list[index]);
        pokeLocalSite(s, out, eg.fUGrid_in_list[index]);
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
	
	size_t u_size = sizeof(sobj)/4;

	size_t count = 0;
	// zero (plus) all boundary gauge links.
	for(size_t offset = 0; offset < f.field.size(); offset++){
		qlat::Coordinate x = geo.coordinate_from_offset(offset);
		void* ptr = &(f.field[offset]);
		for(int mu = 0; mu < 4; mu++){
			if( x[mu] == geo.node_site_expanded[mu]-geo.expansion_left[mu]-1 ) {
				count++;
				memset(ptr+u_size*mu, 0, u_size);
			}
		}
	}

	printf("count = %d\n", count);

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
	for(size_t index = 0; index < eg.fFrbGrid_in_olist.size(); index++){
		int thread_num = omp_get_thread_num();
		scalar_type* lp = (scalar_type*)&left._odata[eg.fFrbGrid_in_olist[index]];
		scalar_type* rp = (scalar_type*)&right._odata[eg.fFrbGrid_in_olist[index]];
		for(size_t w = 0; w < sizeof(vobj)/sizeof(vector_type); w++){
			size_t offset = eg.fFrbGrid_in_ilist[index] + w * Nsimd;
	//		size_t offset = eg.fFrbGrid_in_ilist[index];
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

	parallel_for(size_t index = 0; index < eg.fFrbGrid_out_list.size(); index++){
		sobj s; memset(&s, 0, sizeof(sobj));
		pokeLocalSite(s, in, eg.fFrbGrid_out_list[index]);
	}
}

template<class F>
int local_conjugate_gradient_MdagM(ExpandGrid& eg, LinearOperatorBase<F> &D, const F &src, F &sol, int iterations){

	TIMER("local_CG");
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
	
    if( eg.cUGrid->IsBoss() ){
        printf("local_CG_MdagM: BEFORE starting on NODE #%04d: r2 = %8.4e \n",
                src._grid->ThisRank(), rk2);
    }

//	WilsonFermion5DStatic::dirichlet = true;

	for(int local_loop_count = 0; local_loop_count < iterations; local_loop_count++){
		D.Op(p, mp);
        D.AdjOp(mp, mmp);
	//	mmp = mmp + e * p;

		zero_boundary_fermion(eg, mmp);
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
        
		if ( true ){
			zero_boundary_fermion(eg, sol);
            RealD sol2 = local_norm_sqr(sol);
            if( eg.cUGrid->IsBoss() ){
                printf("local_CG_MdagM: l.i. #%04d on NODE #%04d: r2 = %8.4e psi2 = %8.4e alpha = %8.4e beta = %8.4e Mpk2 = %8.4e \n",
                        local_loop_count, src._grid->ThisRank(), rk2, sol2, alpha, beta, Mpk2);
            }
        }
 		
	}
	
//	WilsonFermion5DStatic::dirichlet = false;
    
	return 0;	
}

template<class F>
int local_CG_MdagM_even_odd(ExpandGrid& eg, LinearOperatorBase<F> &D, const F &src, F &sol, int iterations, int eo)
{

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
		D.Op(p, mp);
        D.AdjOp(mp, mmp);

		zero_boundary_fermion(eg, mmp);
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
	if( (src._grid->ThisProcessorCoor()[0]+src._grid->ThisProcessorCoor()[1]+src._grid->ThisProcessorCoor()[2]+src._grid->ThisProcessorCoor()[3]) % 2 != eo)
		sol = zero;
   	
	norm2(sol);

	return 0;	
}

template<class F>
int MSP_conjugate_gradient(ExpandGrid& eg, LinearOperatorBase<F>& D, LinearOperatorBase<F>& fD, const F& src, F& sol, 
								RealD percent, int f_iter, size_t max_iter = 50000, RealD e = 0.)
{

	RealD alpha, beta, rkzk, pkApk, zkp1rkp1, zkrk;
	RealD r2, src_sqr, tgt_r2;

	F p(src);
	F mp(src);
	F mmp(src);
	F r(src);
	F z(src);

	F f_r(eg.fFrbGrid);
	F f_z(eg.fFrbGrid);

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
	local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter, e); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
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
	    local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter, e); // z_k+1 = M^-1 * r_k+1
//	    local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter); // z_k+1 = M^-1 * r_k+1
		shrink_fermion(eg, f_z, z);
		
		zkp1rkp1 = std::real(innerProduct(z, r));	
//		zkp1rkp1 = -alpha*std::real(innerProduct(z, mmp));	
		
		beta = zkp1rkp1 / rkzk;
		p = z + beta * p; // p_k+1 = z_k+1 + beta * p_k
	
		r2 = norm2(r);

		if ( true ){
		    if(not src._grid->ThisRank()){
				printf("MSPCG/iter.count/r2/target_r2/%%/target_%%: %05d %8.4e %8.4e %8.4e %8.4e\n", k, r2, tgt_r2, std::sqrt(r2/src_sqr), percent);
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

template<class F, class G>
int MSPCG_half(ExpandGrid& eg, LinearOperatorBase<G>& D, LinearOperatorBase<F>& fD, const G& src, G& sol, 
								RealD percent, int f_iter, size_t max_iter = 50000)
{

	RealD alpha, beta, rkzk, pkApk, zkp1rkp1, zkrk;
	RealD r2, src_sqr, tgt_r2;

	G p(src);
	G mp(src);
	G mmp(src);
	G r(src);
	G z(src);

	G f_r(eg.fFrbGrid);
	G f_z(eg.fFrbGrid);
	
	F f_r_F(eg.fFrbGrid_F);
	F f_z_F(eg.fFrbGrid_F);

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

	expand_fermion_D2F(eg, r, f_r_F);
//	precisionChange(f_r_F, f_r);
	local_conjugate_gradient_MdagM(eg, fD, f_r_F, f_z_F, f_iter); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
//	precisionChange(f_z, f_z_F);
	shrink_fermion_F2D(eg, f_z_F, z);

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
   
		{
		TIMER("local_CG/padding")
		expand_fermion_D2F(eg, r, f_r_F);
//		precisionChange(f_r_F, f_r);
		local_conjugate_gradient_MdagM(eg, fD, f_r_F, f_z_F, f_iter); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
//		precisionChange(f_z, f_z_F);
		shrink_fermion_F2D(eg, f_z_F, z);
		}

		zkp1rkp1 = std::real(innerProduct(z, r));	
//		zkp1rkp1 = -alpha*std::real(innerProduct(z, mmp));	
		
		beta = zkp1rkp1 / rkzk;
		p = z + beta * p; // p_k+1 = z_k+1 + beta * p_k
	
		r2 = norm2(r);

		if ( true ){
		    if(not src._grid->ThisRank()){
				printf("MSPCG/iter.count/r2/target_r2/%%/target_%%: %05d %8.4e %8.4e %8.4e %8.4e\n", k, r2, tgt_r2, std::sqrt(r2/src_sqr), percent);
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

template<class F>
int DD_CG(ExpandGrid& eg, LinearOperatorBase<F>& D, LinearOperatorBase<F>& fD, const F& src, F& sol, 
					RealD percent, int f_iter, size_t max_iter = 50000)
{

	RealD alpha, beta, rkzk, pkApk, zkp1rkp1, zkrk;
	RealD r2, src_sqr, tgt_r2;

	F p(src);
	F mp(src);
	F mmp(src);
	F r(src);
	F z(src);

	F f_r(eg.fFrbGrid);
	F f_z(eg.fFrbGrid);

// DD specifically
	F f_z2(eg.fFrbGrid);
	
	F f_tmp(eg.fFrbGrid);
	F f_tmp2(eg.fFrbGrid);
	
	F r_tmp(src);
	F r_tmp2(src);

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
	
//	local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method

// --- DD
		
//	local_conjugate_gradient_MdagM(eg, fD, f_r, f_tmp, f_iter);
	local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter);
//	shrink_fermion(eg, f_tmp, r_tmp);
	
//	WilsonFermion5DStatic::dirichlet = true;
//	fD.Op(f_tmp, f_tmp2);
//	fD.AdjOp(f_tmp2, f_tmp);
//	WilsonFermion5DStatic::dirichlet = false;
	
//	D.Op(r_tmp, r_tmp2);
//	D.AdjOp(r_tmp2, r_tmp);
//
//	expand_fermion(eg, r_tmp, f_tmp2);
//	
//	f_r = f_r - f_tmp2;
//	zero_boundary_fermion(eg, f_r);
//	
//	local_CG_MdagM_even_odd(eg, fD, f_r, f_tmp2, f_iter, 1);
//
//	f_z = f_tmp + f_tmp2;

// --- DD

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
//	    local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter); // z_k+1 = M^-1 * r_k+1

// --- DD
//		local_CG_MdagM_even_odd(eg, fD, f_r, f_tmp, f_iter, 0);
		local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter);
//		shrink_fermion(eg, f_tmp, r_tmp);
		
//		WilsonFermion5DStatic::dirichlet = true;
//		fD.Op(f_tmp, f_tmp2);
//		fD.AdjOp(f_tmp2, f_tmp);
//		WilsonFermion5DStatic::dirichlet = false;
	
//		D.Op(r_tmp, r_tmp2);
//		D.AdjOp(r_tmp2, r_tmp);
//
//		expand_fermion(eg, r_tmp, f_tmp2);
//	
//		f_r = f_r - f_tmp2;
//		zero_boundary_fermion(eg, f_r);
//	
//		local_CG_MdagM_even_odd(eg, fD, f_r, f_tmp2, f_iter, 1);
//
//		f_z = f_tmp + f_tmp2;
// --- DD

//	    local_conjugate_gradient_MdagM(eg, fD, f_r, f_z, f_iter); // z_k+1 = M^-1 * r_k+1
		shrink_fermion(eg, f_z, z);
		
//		zkp1rkp1 = std::real(innerProduct(z, r));	
		zkp1rkp1 = -alpha*std::real(innerProduct(z, mmp));	
		
		beta = zkp1rkp1 / rkzk;
		p = z + beta * p; // p_k+1 = z_k+1 + beta * p_k
	
		r2 = norm2(r);

		if ( true ){
		    if(not src._grid->ThisRank()){
				printf("MSPCG/iter.count/r2/target_r2/%%/target_%%: %05d %8.4e %8.4e %8.4e %8.4e\n", k, r2, tgt_r2, std::sqrt(r2/src_sqr), percent);
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
