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
#ifndef GRID_EXPAND_GRID_H
#define GRID_EXPAND_GRID_H

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
	
	std::array<int, 4> expansion_left;
	std::array<int, 4> expansion_right;
	std::array<int, 4> inner_expansion_left;
	std::array<int, 4> inner_expansion_right;

	std::vector<std::vector<int>> fFrbGrid_in_list;
	std::vector<std::vector<int>> cFrbGrid_in_list;
	std::vector<std::vector<int>> fFrbGrid_out_list;

	std::vector<size_t> fFrbGrid_in_olist;
	std::vector<size_t> fFrbGrid_in_ilist;
	
	std::vector<std::vector<int>> fUGrid_in_list;
	std::vector<std::vector<int>> cUGrid_in_list;
	
	std::vector<std::vector<int>> fFrbGrid_inner_out_list;

	inline void init(	GridCartesian* _cUGrid, 
						std::array<int, 4>& _expansion_left, 
						std::array<int, 4>& _expansion_right,
						std::array<int, 4>& _inner_expansion_left,
						std::array<int, 4>& _inner_expansion_right, 
						int Ls
					){

		expansion_left = _expansion_left;		
		expansion_right = _expansion_right;		
		inner_expansion_left = _inner_expansion_left;		
		inner_expansion_right = _inner_expansion_right;		
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
			expanded_latt_size[i] = latt_size[i]/mpi_layout[i] + expansion_left[i] + expansion_right[i];
//   		expanded_latt_size[i] = latt_size[i] + 2*expansion*mpi_layout[i];
			assert(latt_size[i] > expanded_latt_size[i]); // qlat does NOT support local volume larger than global volume
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

		fFrbGrid_inner_out_list.resize(0);

		size_t count_inner = 0;

// fermion copy list. Maybe by directly get index for memory we will get a better performance?
		for(int index = 0; index < fFrbGrid->lSites(); index++){
			std::vector<int> lcoor_out(5), lcoor_in(5), full_lcoor_out(5), full_lcoor_in(5);
			fFrbGrid->LocalIndexToLocalCoor(index, lcoor_out);

			lcoor_in[0] = lcoor_out[0]; // s: Ls
			lcoor_in[1] = lcoor_out[1] - expansion_left[0]/2; // x is the _checker_dim
			for(int d = 2; d < 5; d++){
				lcoor_in[d] = lcoor_out[d] - expansion_left[d-1];
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
				
				
				if( lcoor_in[1] < -inner_expansion_left[0]/2 or lcoor_in[1] >= cFrbGrid->_ldimensions[1]+inner_expansion_right[0]/2 or 
					lcoor_in[2] < -inner_expansion_left[1] or lcoor_in[2] >= cFrbGrid->_ldimensions[2]+inner_expansion_right[1] or 
					lcoor_in[3] < -inner_expansion_left[2] or lcoor_in[3] >= cFrbGrid->_ldimensions[3]+inner_expansion_right[2] or 
					lcoor_in[4] < -inner_expansion_left[3] or lcoor_in[4] >= cFrbGrid->_ldimensions[4]+inner_expansion_right[3]){
					
					fFrbGrid_inner_out_list.push_back(full_lcoor_out);
					count_inner++;
				}

			}
		}

		printf("count_inner = %d\n", count_inner);

// gauge copy list.
		fUGrid_in_list.resize(0);
		cUGrid_in_list.resize(0);
		for(int index = 0; index < fUGrid->lSites(); index++){
			std::vector<int> lcoor_out(4), lcoor_in(4);
			fUGrid->LocalIndexToLocalCoor(index, lcoor_out);

			for(int d = 0; d < 4; d++){
				lcoor_in[d] = lcoor_out[d] - expansion_left[d];
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

template<class double_type, class float_type>
void expand_fermion_D2F_qlat(ExpandGrid& eg, const Lattice<double_type>& in, Lattice<float_type>& out){

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
	
	std::vector<Fsobj> out_lex(out._grid->lSites());
  	unvectorizeToLexOrdArray(out_lex, out);

	qlat::Coordinate global_size(in._grid->_gdimensions[1], in._grid->_gdimensions[2], in._grid->_gdimensions[3], in._grid->_gdimensions[4]);
	qlat::Geometry geo; geo.init(global_size, in._grid->_gdimensions[0]); // multiplicity = in._grid->_gdimensions[0]

	qlat::Coordinate e_left(eg.expansion_left[0]/2, eg.expansion_left[1], eg.expansion_left[2], eg.expansion_left[3]); // x direction is the checkerboard dimension so the /2.
	qlat::Coordinate e_right(eg.expansion_right[0]/2, eg.expansion_right[1], eg.expansion_right[2], eg.expansion_right[3]); // x direction is the checkerboard dimension so the /2.
	geo.resize(e_left, e_right);

	qlat::Field<Fsobj> f; f.init(geo);
	assert(f.field.size() == out._grid->lSites()); 
	memcpy(f.field.data(), out_lex.data(), f.field.size()*sizeof(Fsobj));
	// DO comm.
//	refresh_expanded_m2(f);
	refresh_expanded(f);
	
	memcpy(out_lex.data(), f.field.data(), f.field.size()*sizeof(Fsobj));
	vectorizeFromLexOrdArray(out_lex, out);

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

	qlat::Coordinate e_left(eg.expansion_left[0]/2, eg.expansion_left[1], eg.expansion_left[2], eg.expansion_left[3]); // x direction is the checkerboard dimension so the /2.
	qlat::Coordinate e_right(eg.expansion_right[0]/2, eg.expansion_right[1], eg.expansion_right[2], eg.expansion_right[3]); // x direction is the checkerboard dimension so the /2.
	geo.resize(e_left, e_right);

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
	
	qlat::Coordinate e_left(eg.expansion_left[0], eg.expansion_left[1], eg.expansion_left[2], eg.expansion_left[3]);
	qlat::Coordinate e_right(eg.expansion_right[0], eg.expansion_right[1], eg.expansion_right[2], eg.expansion_right[3]);
	geo.resize(e_left, e_right);

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
//			// TODO: !!! for debug ONLY
//			
//			if( x[mu] == geo.node_site_expanded[mu]-geo.expansion_left[mu]-2 ) {
//				count++;
//				memset(ptr+u_size*mu, 0, u_size);
//			}
//			
//			if( x[mu] == -2 ) {
//				count++;
//				memset(ptr+u_size*mu, 0, u_size);
//			}
//			
//			if( x[mu] == -1 ) {
//				count++;
//				memset(ptr+u_size*mu, 0, u_size);
//			}
		}
	}
	
	if(in._grid->IsBoss()){
		printf("A total of %d boundary gauge links are zeroed.\n", count);
	}
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

template<class sc>
void zero_boundary_fermion_inner(ExpandGrid& eg, Lattice<sc>& in){

	TIMER("zero_boundary_fermion_inner");

	typedef typename sc::scalar_object sobj;

	parallel_for(size_t index = 0; index < eg.fFrbGrid_inner_out_list.size(); index++){
		sobj s; memset(&s, 0, sizeof(sobj));
		pokeLocalSite(s, in, eg.fFrbGrid_inner_out_list[index]);
	}
}

#endif
