#include <Grid/Grid.h>
#include <cassert>

#include <qlat/qlat.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template <class d>
struct scal {
	d internal;
};

Gamma::Algebra Gmu[] = {Gamma::Algebra::GammaX, Gamma::Algebra::GammaY, 
						Gamma::Algebra::GammaZ, Gamma::Algebra::GammaT};
/*
template<class vobj>
void copy_local_site(const Lattice<vobj>& from, Lattice<vobj>& to,
                    std::vector<int>& from_site, std::vector<int>& to_site){

	typedef typename vobj::scalar_type scalar_type;

    GridBase *from_grid = from._grid;
    GridBase *to_grid = to._grid;
    assert(from_grid->Nsimd() == to_grid->Nsimd());

	// forcing to use the base class version of i/oIndex since the coordinate is the checkerboarded one.
    int from_idx = from_grid->GridBase::iIndex(from_site);
    int to_idx = to_grid->GridBase::iIndex(to_site);
    int from_odx = from_grid->GridBase::oIndex(from_site);
    int to_odx = to_grid->GridBase::oIndex(to_site);

   	scalar_type* to_vp = (scalar_type*)&to._odata[to_odx];
    scalar_type* from_vp = (scalar_type*)&from._odata[from_odx];
	to_vp[to_idx] = from_vp[from_idx];

    return;
};

template<class vobj>
void set_local_site(Lattice<vobj>& to, std::vector<int>& to_site, int ch){

	typedef typename vobj::scalar_type scalar_type;
    
	GridBase *to_grid = to._grid;

	// forcing to use the base class version of i/oIndex since the coordinate is the checkerboarded one.
    int to_idx = to_grid->GridBase::iIndex(to_site);
    int to_odx = to_grid->GridBase::oIndex(to_site);

    scalar_type* to_vp = (scalar_type*)&to._odata[to_odx];
    memset(to_vp+to_idx, ch, sizeof(scalar_type));

    return;
};
*/

template<class sc>
void expandLatticeFermion(const Lattice<sc>& in, Lattice<sc>& out){
	
	GridStopWatch watch;
	watch.Start();

	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename sc::scalar_object sobj;
	
	GridBase* gbin = in._grid;
	GridBase* gbout = out._grid;
	assert(gbin->_ndimension == gbout->_ndimension);
	assert(gbin->_ndimension == 5);
	assert(gbin->_ldimensions[0] == gbout->_ldimensions[0]);
	assert(gbin->_ldimensions[1]+2 == gbout->_ldimensions[1]); // x is the _checker_dim
	for(int d = 2; d < 5; d++){ // Ls is the 0th index here.
		assert(gbin->_ldimensions[d]+4 == gbout->_ldimensions[d]);
	}

	out.checkerboard = in.checkerboard;

//	size_t copy_count = 0;
//	size_t set_count = 0;
	parallel_for(int index = 0; index < out._grid->lSites(); index++){
		std::vector<int> lcoor_out(5), lcoor_in(5), full_lcoor_out(5), full_lcoor_in(5);
		out._grid->LocalIndexToLocalCoor(index, lcoor_out);
	
		lcoor_in[0] = lcoor_out[0]; // s: Ls
		lcoor_in[1] = lcoor_out[1]-1; // x is the _checker_dim
		for(int d = 2; d < 5; d++){
			lcoor_in[d] = lcoor_out[d]-2;
		}
		
		for(int cb = 0; cb < 2; cb++){
		
		for(int d = 0; d < 5; d++){ // Ls is the 0th index here.
			full_lcoor_out[d] = lcoor_out[d];
			full_lcoor_in[d] = lcoor_in[d];
		}
	
		full_lcoor_out[1] = lcoor_out[1]*2+cb;
		full_lcoor_in[1] = lcoor_in[1]*2+cb;

		if(out.checkerboard != gbout->CheckerBoard(full_lcoor_out)) continue;
		

		if( lcoor_in[1]>=0 and lcoor_in[1]<gbin->_ldimensions[1] and
			lcoor_in[2]>=0 and lcoor_in[2]<gbin->_ldimensions[2] and
			lcoor_in[3]>=0 and lcoor_in[3]<gbin->_ldimensions[3] and
			lcoor_in[4]>=0 and lcoor_in[4]<gbin->_ldimensions[4]){
		
			sobj s;
			peekLocalSite(s, in, full_lcoor_in);
//			peek_local_site(s, in, lcoor_in);
			pokeLocalSite(s, out, full_lcoor_out);
//			poke_local_site(s, out, lcoor_out);
//			copy_local_site(in, out, lcoor_in, lcoor_out);			
//			copy_count++;
		}else{
			sobj s; memset(&s, 0, sizeof(sobj));
			pokeLocalSite(s, out, full_lcoor_out);
//			poke_local_site(s, out, lcoor_out);
//			set_local_site(out, lcoor_out, 0);
//			set_count++;
		}

		}
	}
//	std::cout << GridLogMessage << copy_count << "\tset:\t " << set_count << std::endl;
	watch.Stop();

	std::cout << GridLogMessage << "Total fermion expansion time : " << watch.Elapsed() << std::endl;
}

template<class vobj>
void shrinkLatticeFermion(const Lattice<vobj>& in, Lattice<vobj>& out){
	
	// Simply shrink and then copy/merge.
	typedef typename vobj::scalar_object sobj;
	
	GridBase* gb_in = in._grid;
	GridBase* gb_out = out._grid;
	assert(gb_in->_ndimension == gb_out->_ndimension);
	assert(gb_in->_ndimension == 5);
	assert(gb_in->_ldimensions[0] == gb_out->_ldimensions[0]);
	assert(gb_in->_ldimensions[1]-2 == gb_out->_ldimensions[1]); // x is the _checker_dim
	for(int d = 2; d < 5; d++){ // Ls is the 0th index here.
		assert(gb_in->_ldimensions[d]-4 == gb_out->_ldimensions[d]);
	}

	out.checkerboard = in.checkerboard;

	parallel_for(int index = 0; index < out._grid->lSites(); index++){
		std::vector<int> lcoor_out(5), lcoor_in(5), full_lcoor_out(5), full_lcoor_in(5);
		out._grid->LocalIndexToLocalCoor(index, lcoor_out);
	
		lcoor_in[0] = lcoor_out[0]; // s: Ls
		lcoor_in[1] = lcoor_out[1]+1; // x is the _checker_dim
		for(int d = 2; d < 5; d++){
			lcoor_in[d] = lcoor_out[d]+2;
		}
			
		for(int d = 0; d < 5; d++){ // Ls is the 0th index here.
			full_lcoor_out[d] = lcoor_out[d];
			full_lcoor_in[d] = lcoor_in[d];
		}

		for(int cb = 0; cb < 2; cb++){
		
			full_lcoor_out[1] = lcoor_out[1]*2+cb;
			full_lcoor_in[1] = lcoor_in[1]*2+cb;

			if(out.checkerboard != gb_out->CheckerBoard(full_lcoor_out)) continue;
	
			sobj s;
			peekLocalSite(s, in, full_lcoor_in);
			pokeLocalSite(s, out, full_lcoor_out);
		}
	}
}

template<class vobj>
void expandLatticeGaugeField(const Lattice<vobj>& in, Lattice<vobj>& out){
	
	// Simply expand and then copy/merge.
	// Set the Boundary sites to zero.
	typedef typename vobj::scalar_object sobj;
	
	GridBase* gbin = in._grid;
	GridBase* gbout = out._grid;
	assert(in._grid->_ndimension == out._grid->_ndimension);
	assert(in._grid->_ndimension == 4);
	for(int d = 0; d < 4; d++){ // Ls is the 0th index here.
		assert(in._grid->_ldimensions[d]+4 == out._grid->_ldimensions[d]);
	}

	for(int index = 0; index < out._grid->lSites(); index++){
		std::vector<int> lcoor_out(4), lcoor_in(4), gcoor_in(4);
		out._grid->LocalIndexToLocalCoor(index, lcoor_out);
	
		for(int d = 0; d < 4; d++){
			lcoor_in[d] = lcoor_out[d]-2;
		}

		if( lcoor_in[1]>=0 and lcoor_in[1]<gbin->_ldimensions[0] and
			lcoor_in[2]>=0 and lcoor_in[2]<gbin->_ldimensions[1] and
			lcoor_in[3]>=0 and lcoor_in[3]<gbin->_ldimensions[2] and
			lcoor_in[0]>=0 and lcoor_in[0]<gbin->_ldimensions[3]){
			// local gauge links
			sobj s;
			peekLocalSite(s, in, lcoor_in);
			pokeLocalSite(s, out, lcoor_out);
		}else{
			// off-node gauge links
			in._grid->ProcessorCoorLocalCoorToGlobalCoor(in._grid->_processor_coor, lcoor_in, gcoor_in);
			for(int mu = 0; mu < 4; mu++){
				gcoor_in[mu] = (gcoor_in[mu]+in._grid->_gdimensions[mu])%in._grid->_gdimensions[mu];
			}
			sobj s; // memset(&s, 0, sizeof(sobj));
			peekSite(s, in, gcoor_in);
			pokeLocalSite(s, out, lcoor_out);
		}

	}

}

template<class vobj>
void expand_lattice_gauge_field_qlat(const Lattice<vobj>& in, Lattice<vobj>& out){
	
	// Simply expand and then copy/merge.
	typedef typename vobj::scalar_object sobj;
	
	GridBase* gbin = in._grid;
	GridBase* gbout = out._grid;
	assert(in._grid->_ndimension == out._grid->_ndimension);
	assert(in._grid->_ndimension == 4);
	for(int d = 0; d < 4; d++){ // Ls is the 0th index here.
		assert(in._grid->_ldimensions[d]+4 == out._grid->_ldimensions[d]);
	}

	for(int index = 0; index < out._grid->lSites(); index++){
		std::vector<int> lcoor_out(4), lcoor_in(4), gcoor_in(4);
		out._grid->LocalIndexToLocalCoor(index, lcoor_out);
	
		for(int d = 0; d < 4; d++){
			lcoor_in[d] = lcoor_out[d]-2;
		}

		if( lcoor_in[1]>=0 and lcoor_in[1]<gbin->_ldimensions[0] and
			lcoor_in[2]>=0 and lcoor_in[2]<gbin->_ldimensions[1] and
			lcoor_in[3]>=0 and lcoor_in[3]<gbin->_ldimensions[2] and
			lcoor_in[0]>=0 and lcoor_in[0]<gbin->_ldimensions[3]){
			// local gauge links
			sobj s;
			peekLocalSite(s, in, lcoor_in);
			pokeLocalSite(s, out, lcoor_out);
		}else{
			// off-node gauge links
			sobj s; memset(&s, 1, sizeof(sobj));
			pokeLocalSite(s, out, lcoor_out);
		}

	}
	
	std::vector<sobj> out_lex(out._grid->lSites());
  	unvectorizeToLexOrdArray(out_lex, out);

	qlat::Coordinate global_size(in._grid->_gdimensions[0], in._grid->_gdimensions[1], in._grid->_gdimensions[2], in._grid->_gdimensions[3]);
	qlat::Geometry geo; geo.init(global_size, 1); // packing gauge links in four directions together, therefore multiplicity=1

	qlat::Coordinate expansion(2, 2, 2, 2);
	geo.resize(expansion, expansion);

	qlat::Field<sobj> f; f.init(geo);
	assert(f.field.size()*1 == out._grid->lSites()); // 1 for multiplicity
	memcpy(f.field.data(), out_lex.data(), f.field.size()*sizeof(sobj));
	// DO comm.
	refresh_expanded(f);
	
	memcpy(out_lex.data(), f.field.data(), f.field.size()*sizeof(sobj));
	vectorizeFromLexOrdArray(out_lex, out);

}
int main(int argc, char** argv) {
	Grid_init(&argc, &argv);

	GridLogIterative.Active(1);

	const int Ls = 12;

	GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
	GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
	GridCartesian* FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
	GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

	FrbGrid->show_decomposition();

	qlat::Coordinate node_coor(UGrid->ThisProcessorCoor()[0], UGrid->ThisProcessorCoor()[1], UGrid->ThisProcessorCoor()[2], UGrid->ThisProcessorCoor()[3]);
	qlat::Coordinate node_size(GridDefaultMpi()[0], GridDefaultMpi()[1], GridDefaultMpi()[2], GridDefaultMpi()[3]);
	qlat::begin(qlat::index_from_coordinate(node_coor, node_size), node_size);
	printf("Node #%03d(grid): %02dx%02dx%02dx%02d ; #%03d(qlat): %02dx%02dx%02dx%02d\n", UGrid->ThisRank(), 
				UGrid->ThisProcessorCoor()[0], UGrid->ThisProcessorCoor()[1], UGrid->ThisProcessorCoor()[2], UGrid->ThisProcessorCoor()[3], 
				qlat::get_id_node(), qlat::get_coor_node()[0], qlat::get_coor_node()[1], qlat::get_coor_node()[2], qlat::get_coor_node()[3]);
	
	std::vector<int> seeds4({1, 2, 3, 4});
	std::vector<int> seeds5({5, 6, 7, 8});
	GridParallelRNG RNG5(FGrid);
	RNG5.SeedFixedIntegers(seeds5);
	GridParallelRNG RNG4(UGrid);
	RNG4.SeedFixedIntegers(seeds4);

	LatticeFermion src(FGrid);
	random(RNG5, src);
	LatticeFermion result(FGrid);
	result = zero;
	LatticeGaugeField Umu(UGrid);
	
	std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
	std::vector<int> mpi_layout  = GridDefaultMpi();
	std::vector<int> latt_size = GridDefaultLatt();
	std::vector<int> expanded_latt_size = GridDefaultLatt();
	for(int i = 0; i < 4; i++){
		expanded_latt_size[i] = latt_size[i] + 4*mpi_layout[i];
	}

	GridCartesian* ExpandedUGrid = SpaceTimeGrid::makeFourDimGrid(expanded_latt_size, simd_layout, mpi_layout);
	GridRedBlackCartesian* ExpandedUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(ExpandedUGrid);
	GridCartesian* ExpandedFGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, ExpandedUGrid);
	GridRedBlackCartesian* ExpandedFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, ExpandedUGrid);

	ExpandedFrbGrid->show_decomposition();

	int orthodir=3;
	int orthosz =latt_size[orthodir];

	FieldMetaData header;
	std::string file("/global/homes/j/jiquntu/configurations/32x64x12ID_b1.75_mh0.045_ml0.0001/configurations/ckpoint_lat.160");
	NerscIO::readConfiguration(Umu, header, file);

//	SU3::HotConfiguration(RNG4, Umu);

	std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << "   Ls: " << Ls << std::endl;

	RealD mass = 0.0001;
	RealD M5 = 1.8;
	MobiusFermionR DMobius(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, 22./12., 10./12.);
	DMobius.ZeroCounters();

	LatticeFermion src_o(FrbGrid);
	LatticeFermion result_o(FrbGrid);
	pickCheckerboard(Odd, src_o, src);
	result_o = zero;

	LatticeFermion psi_o(ExpandedFrbGrid);
	expandLatticeFermion(src_o, psi_o);
	
	LatticeFermion src_o_check(FrbGrid);
	shrinkLatticeFermion(psi_o, src_o_check);
	std::cout << "Fermion expanded." << std::endl;
	std::cout << GridLogMessage << norm2(src_o) << " \t " << norm2(psi_o) << " \t " << norm2(src_o_check) << std::endl;

	LatticeGaugeField expandedUmu(ExpandedUGrid);
	expand_lattice_gauge_field_qlat(Umu, expandedUmu);
	std::cout << "Gauge field expanded." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << " \t " << WilsonLoops<PeriodicGimplR>::avgPlaquette(expandedUmu) << std::endl;

	MobiusFermionR expandedDMobius(expandedUmu, *ExpandedFGrid, *ExpandedFrbGrid, *ExpandedUGrid, *ExpandedUrbGrid, mass, M5, 22./12., 10./12.);
//	LatticeFermion src_o_dup(FrbGrid);
//	localConvert(src_o, src_o_dup);
//	std::cout << GridLogMessage << norm2(src_o) << " \t " << norm2(src_o_dup) << std::endl;

	GridStopWatch CGTimer;

	SchurDiagMooeeOperator<MobiusFermionR, LatticeFermion> HermOpEO(DMobius);
	SchurDiagMooeeOperator<MobiusFermionR, LatticeFermion> expandedHermOpEO(expandedDMobius);
	ConjugateGradient<LatticeFermion> CG(1.0e-10, 30000, 0);// switch off the assert

	LatticeFermion Mdag_src_o(FrbGrid);
	HermOpEO.AdjOp(src_o, Mdag_src_o);

	LatticeFermion x(FrbGrid);
	pickCheckerboard(Odd, x, src);
	LatticeFermion y(FrbGrid);
	LatticeFermion ydd(FrbGrid);
	if(UGrid->ThisRank() != 0) x = zero;

	LatticeFermion xd(ExpandedFrbGrid);
	LatticeFermion yd(ExpandedFrbGrid);

	expandLatticeFermion(x, xd);

	HermOpEO.AdjOp(x, y);
	
	WilsonFermion5DStatic::dirichlet = true;
	expandedHermOpEO.AdjOp(xd, yd);
	WilsonFermion5DStatic::dirichlet = false;
	shrinkLatticeFermion(yd, ydd);	
	
	if(UGrid->ThisRank() != 0) y = zero;

	std::cout << GridLogMessage << "|y - ydd|**2 : " << norm2(y) << "\t" << norm2(ydd) << std::endl;
	
	CGTimer.Start();
//	CG(HermOpEO, Mdag_src_o, result_o);
	CGTimer.Stop();

	std::cout << GridLogMessage << "Total CG time : " << CGTimer.Elapsed() << std::endl;

	std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
	DMobius.Report();

	Grid_finalize();
}
