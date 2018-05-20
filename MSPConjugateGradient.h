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
#include "expand-grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class F>
int local_conjugate_gradient_MdagM(ExpandGrid& eg, LinearOperatorBase<F> &D, const F &src, F &sol, int iterations, RealD e, RealD percent = 1e-3){

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
	RealD Mpk2_imag;
	RealD tgt_r2 = rk2 * percent * percent;

    if( eg.cUGrid->IsBoss() ){
        printf("local_CG_MdagM: BEFORE starting on NODE #%04d: r2 = %8.4e \n",
                src._grid->ThisRank(), rk2);
    }

//	WilsonFermion5DStatic::dirichlet = true;

	for(int local_loop_count = 0; local_loop_count < iterations; local_loop_count++){
		{
		TIMER("local_dslash");
		D.Op(p, mp);
        D.AdjOp(mp, mmp);
		mmp = mmp + e * p;
		}
		zero_boundary_fermion(eg, mmp);
//		zero_boundary_fermion(eg, p); // not necessary?
        Mpk2 = std::real(local_inner_product(p, mmp));
        Mpk2_imag = std::imag(local_inner_product(p, mmp));
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
                printf("local_CG_MdagM: l.i. #%04d on NODE #%04d: r2 = %8.4e psi2 = %8.4e alpha = %8.4e beta = %8.4e Mpk2 = %8.4e imag = %+8.4e \n",
                        local_loop_count, src._grid->ThisRank(), rk2, sol2, alpha, beta, Mpk2, Mpk2_imag);
            }
        }

		if( rk2 < tgt_r2 ){
			int max; 
			MPI_Reduce(&local_loop_count, &max, 1, MPI_INT, MPI_MAX, 0, eg.cUGrid->communicator);
			if( eg.cUGrid->IsBoss() ){
                printf("local_CG EXITS after reaching precision %% = %8.4e: maximum iteration count = %04d\n", percent, max);
            }
			return max; 
		}
 		
	}
	
//	WilsonFermion5DStatic::dirichlet = false;
    
	return iterations;	
}

template<class F>
int local_CG_pad(ExpandGrid& eg, LinearOperatorBase<F> &D, const F &src, F &sol, int iterations, RealD e, RealD percent = 1e-3){

	TIMER("local_CG_pad");
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
	RealD Mpk2_imag;
	RealD tgt_r2 = rk2 * percent * percent;

    if( eg.cUGrid->IsBoss() ){
        printf("local_CG_MdagM: BEFORE starting on NODE #%04d: r2 = %8.4e \n",
                src._grid->ThisRank(), rk2);
    }

//	WilsonFermion5DStatic::dirichlet = true;

	for(int local_loop_count = 0; local_loop_count < iterations; local_loop_count++){
		D.Op(p, mp);
        D.AdjOp(mp, mmp);
		mmp = mmp + e * p;

		zero_boundary_fermion_inner(eg, mmp);
//		zero_boundary_fermion(eg, p); // not necessary?
        Mpk2 = std::real(local_inner_product(p, mmp));
        Mpk2_imag = std::imag(local_inner_product(p, mmp));
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
			zero_boundary_fermion_inner(eg, sol);
            RealD sol2 = local_norm_sqr(sol);
            if( eg.cUGrid->IsBoss() ){
                printf("local_CG_MdagM: l.i. #%04d on NODE #%04d: r2 = %8.4e psi2 = %8.4e alpha = %8.4e beta = %8.4e Mpk2 = %8.4e imag = %8.4e \n",
                        local_loop_count, src._grid->ThisRank(), rk2, sol2, alpha, beta, Mpk2, Mpk2_imag);
            }
        }

		if( rk2 < tgt_r2 ){
			int max; 
			MPI_Reduce(&local_loop_count, &max, 1, MPI_INT, MPI_MAX, 0, eg.cUGrid->communicator);
			if( eg.cUGrid->IsBoss() ){
                printf("local_CG EXITS after reaching precision %% = %8.4e: maximum iteration count = %04d\n", percent, max);
            }
			return max; 
		}
 		
	}
	
//	WilsonFermion5DStatic::dirichlet = false;
    
	return iterations;	
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
								RealD percent, int f_iter, size_t max_iter = 50000, RealD e = 0.)
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
	local_conjugate_gradient_MdagM(eg, fD, f_r_F, f_z_F, f_iter, e); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
//	precisionChange(f_z, f_z_F);
	shrink_fermion_F2D(eg, f_z_F, z);

	p = z; // p0 = z0

	for(size_t k = 0; k < max_iter; k++){

		TIMER("MSPCG_iteration");

		rkzk = std::real(innerProduct(r, z));

		{
		TIMER("global_dslash")
		D.Op(p, mp);
		D.AdjOp(mp, mmp);
		}
		pkApk = std::real(innerProduct(p, mmp));
		alpha = rkzk / pkApk; // alpha_k

		sol = sol + alpha * p; // x_k+1 = x_k + alpha * p_k
		r = r - alpha * mmp; // r_k+1 = r_k - alpha * Ap_k
   
		{
		TIMER("local_CG/padding")
		expand_fermion_D2F(eg, r, f_r_F);
//		precisionChange(f_r_F, f_r);
		local_conjugate_gradient_MdagM(eg, fD, f_r_F, f_z_F, f_iter, e); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
//		precisionChange(f_z, f_z_F);
		shrink_fermion_F2D(eg, f_z_F, z);
		}

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

template<class F>
void shift_ff(const F& in, F& out, bool forward){
	if(forward){
		out = Cshift( in, 1, out._grid->_ldimensions[1]*2/2);
		out = Cshift(out, 2, out._grid->_ldimensions[2]/2);
		out = Cshift(out, 3, out._grid->_ldimensions[3]/2);
		out = Cshift(out, 4, out._grid->_ldimensions[4]/2);
		return;
	}else{
		out = Cshift( in, 1, -out._grid->_ldimensions[1]*2/2);
		out = Cshift(out, 2, -out._grid->_ldimensions[2]/2);
		out = Cshift(out, 3, -out._grid->_ldimensions[3]/2);
		out = Cshift(out, 4, -out._grid->_ldimensions[4]/2);
		return;
	}
}

template<class F, class G>
int MSPCG_shift(ExpandGrid& eg, LinearOperatorBase<G>& D, LinearOperatorBase<F>& fD, 
						LinearOperatorBase<G>& shifted_D, LinearOperatorBase<F>& shifted_fD, const G& src, G& sol, 
								RealD percent, int f_iter, size_t max_iter = 50000, RealD e = 0.)
{

	size_t shoushuliangduan = 500;

	G shifted_src(src);
	G shifted_sol(sol);
	shift_ff(src, shifted_src, true);

	int indicator = -1;
	size_t count = 0;

	while(indicator<0){
		MSPCG_half(eg, D, fD, src, sol, percent, f_iter, shoushuliangduan, e);
		shift_ff(sol, shifted_sol, true);
		indicator = MSPCG_half(eg, shifted_D, shifted_fD, shifted_src, shifted_sol, percent, f_iter, shoushuliangduan, e);
		shift_ff(shifted_sol, sol, false);
	
		count += shoushuliangduan*2;
	}
	return count;
}

template<class F, class G>
int MSPCG_half_v2(ExpandGrid& eg, LinearOperatorBase<G>& D, LinearOperatorBase<F>& fD, const G& src, G& sol, 
								RealD percent, int f_iter, size_t max_iter = 50000, RealD e = 0.)
{

	RealD alpha, beta, rkzk, pkApk, zkp1rkp1, zkrk;
	RealD r2, src_sqr, tgt_r2;

	RealD gamma = +0.6;

	if(src._grid->IsBoss()){
    	printf("MSPCG_half_v2: gamma= %+8.2e\n", gamma);
	}	

	G p(src);
	G mp(src);
	G mmp(src);
	G r(src);
	G z(src);
	
	G mmp_old(src);
	G r_cc(src);
	G z_old(src);

	G f_r(eg.fFrbGrid);
	G f_z(eg.fFrbGrid);
	
	F f_r_F(eg.fFrbGrid_F);
	F f_z_F(eg.fFrbGrid_F);

	D.Op(sol, mp);
	D.AdjOp(mp, mmp);

	mmp_old = mmp;

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
	local_conjugate_gradient_MdagM(eg, fD, f_r_F, f_z_F, f_iter, e); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
//	precisionChange(f_z, f_z_F);
	shrink_fermion_F2D(eg, f_z_F, z);

	p = z; // p0 = z0

	for(size_t k = 0; k < max_iter; k++){

		TIMER("MSPCG_iteration");

		rkzk = std::real(innerProduct(r, z));

		{
		TIMER("global_dslash")
		D.Op(p, mp);
		D.AdjOp(mp, mmp);
		}
		pkApk = std::real(innerProduct(p, mmp));
		alpha = rkzk / pkApk; // alpha_k

		sol = sol + alpha * p; // x_k+1 = x_k + alpha * p_k
		r = r - alpha * mmp; // r_k+1 = r_k - alpha * Ap_k
   
		{
		TIMER("local_CG/padding")
		
		if(k>0){
			r_cc = r - gamma*(mmp - beta*mmp_old);
			expand_fermion_D2F(eg, r_cc, f_r_F);
		}else{	
			expand_fermion_D2F(eg, r, f_r_F);
		}
		mmp_old = mmp;

//		expand_fermion_D2F(eg, r, f_r_F);
//		precisionChange(f_r_F, f_r);
		local_conjugate_gradient_MdagM(eg, fD, f_r_F, f_z_F, f_iter, e); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
//		precisionChange(f_z, f_z_F);
		shrink_fermion_F2D(eg, f_z_F, z);
		
		if(k>0){
			z = z + gamma*z_old;
		}

		}
		z_old = z;

//		zkp1rkp1 = std::real(innerProduct(z, r));	
		zkp1rkp1 = -alpha*std::real(innerProduct(z, mmp));	
		
		beta = zkp1rkp1 / rkzk;
		p = z + beta * p; // p_k+1 = z_k+1 + beta * p_k
	
		r2 = norm2(r);
		RealD z2 = norm2(z);
		if ( true ){
		    if(not src._grid->ThisRank()){
				printf("MSPCG/iter.count/r2/target_r2/%%/target_%%: %05d %8.4e %8.4e %8.4e %8.4e \n", k, r2, tgt_r2, std::sqrt(r2/src_sqr), percent);
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
int MSPCG_pad(ExpandGrid& eg, LinearOperatorBase<G>& D, LinearOperatorBase<F>& fD, const G& src, G& sol, 
								RealD percent, int f_iter, size_t max_iter = 50000, RealD e = 0.)
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

	expand_fermion_D2F_qlat(eg, r, f_r_F);
	zero_boundary_fermion_inner(eg, f_r_F);
//	precisionChange(f_r_F, f_r);
	local_CG_pad(eg, fD, f_r_F, f_z_F, f_iter, e); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
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
		expand_fermion_D2F_qlat(eg, r, f_r_F);
		zero_boundary_fermion_inner(eg, f_r_F);
//		precisionChange(f_r_F, f_r);
		local_CG_pad(eg, fD, f_r_F, f_z_F, f_iter, e); // z0 = M^-1 * r0 // notations follow https://en.wikipedia.org/wiki/Conjugate_gradient_method
//		precisionChange(f_z, f_z_F);
		shrink_fermion_F2D(eg, f_z_F, z);
		}

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
