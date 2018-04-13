#include <Grid/Grid.h>
#include <cassert>

#include <qlat/qlat.h>

#include "MSPConjugateGradient.h"

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
void expand_fermion(ExpandGrid& eg, const Lattice<sc>& in, Lattice<sc>& out){
	
	GridStopWatch watch;
	watch.Start();

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

	watch.Stop();

	std::cout << GridLogMessage << "Total fermion expansion time : " << watch.Elapsed() << std::endl;
}

template<class sc>
void shrink_fermion(ExpandGrid& eg, const Lattice<sc>& in, Lattice<sc>& out){
	
	GridStopWatch watch;
	watch.Start();

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

	watch.Stop();

	std::cout << GridLogMessage << "Total fermion shrinking time : " << watch.Elapsed() << std::endl;
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
int main(int argc, char** argv) {
	Grid_init(&argc, &argv);

	GridLogIterative.Active(1);

	const int Ls = 12;

	GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
//	GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
//	GridCartesian* FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
//	GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);
//
//	FrbGrid->show_decomposition();

	qlat::Coordinate node_coor(UGrid->ThisProcessorCoor()[0], UGrid->ThisProcessorCoor()[1], UGrid->ThisProcessorCoor()[2], UGrid->ThisProcessorCoor()[3]);
	qlat::Coordinate node_size(GridDefaultMpi()[0], GridDefaultMpi()[1], GridDefaultMpi()[2], GridDefaultMpi()[3]);
	qlat::begin(qlat::index_from_coordinate(node_coor, node_size), node_size);
	printf("Node #%03d(grid): %02dx%02dx%02dx%02d ; #%03d(qlat): %02dx%02dx%02dx%02d\n", UGrid->ThisRank(), 
				UGrid->ThisProcessorCoor()[0], UGrid->ThisProcessorCoor()[1], UGrid->ThisProcessorCoor()[2], UGrid->ThisProcessorCoor()[3], 
				qlat::get_id_node(), qlat::get_coor_node()[0], qlat::get_coor_node()[1], qlat::get_coor_node()[2], qlat::get_coor_node()[3]);
	
	ExpandGrid eg; eg.init(UGrid, 2, Ls);
	
	std::vector<int> seeds4({1, 2, 3, 4});
	std::vector<int> seeds5({5, 6, 7, 8});
	GridParallelRNG RNG5(eg.coarse_fermion_grid);
	RNG5.SeedFixedIntegers(seeds5);
	GridParallelRNG RNG4(UGrid);
	RNG4.SeedFixedIntegers(seeds4);

	LatticeFermion src(eg.coarse_fermion_grid);
	random(RNG5, src);
	LatticeFermion result(eg.coarse_fermion_grid);
	result = zero;
	LatticeGaugeField Umu(UGrid);
	
//	std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
//	std::vector<int> mpi_layout  = GridDefaultMpi();
	std::vector<int> latt_size = GridDefaultLatt();
//	std::vector<int> expanded_latt_size = GridDefaultLatt();
//	for(int i = 0; i < 4; i++){
//		expanded_latt_size[i] = latt_size[i] + 4*mpi_layout[i];
//	}
//
//	GridCartesian* ExpandedUGrid = SpaceTimeGrid::makeFourDimGrid(expanded_latt_size, simd_layout, mpi_layout);
//	GridRedBlackCartesian* ExpandedUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(ExpandedUGrid);
//	GridCartesian* ExpandedFGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, ExpandedUGrid);
//	GridRedBlackCartesian* ExpandedFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, ExpandedUGrid);
//
//	ExpandedFrbGrid->show_decomposition();

	int orthodir=3;
	int orthosz =latt_size[orthodir];

	FieldMetaData header;
	std::string file("/global/homes/j/jiquntu/configurations/32x64x12ID_b1.75_mh0.045_ml0.0001/configurations/ckpoint_lat.160");
	NerscIO::readConfiguration(Umu, header, file);

//	SU3::HotConfiguration(RNG4, Umu);

	std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << "   Ls: " << Ls << std::endl;

	RealD mass = 0.0001;
	RealD M5 = 1.8;
	MobiusFermionR DMobius(Umu, *eg.coarse_fermion_grid, *eg.coarse_fermion_rb_grid, *eg.coarse_gauge_grid, *eg.coarse_gauge_rb_grid, mass, M5, 22./12., 10./12.);
	DMobius.ZeroCounters();
	
	LatticeFermion src_o(eg.coarse_fermion_rb_grid);
	LatticeFermion result_o(eg.coarse_fermion_rb_grid);
	pickCheckerboard(Odd, src_o, src);
	result_o = zero;


	LatticeFermion psi_o(eg.fine_fermion_rb_grid);
	expand_fermion(eg, src_o, psi_o);
	
	LatticeFermion src_o_check(eg.coarse_fermion_rb_grid);
	shrink_fermion(eg, psi_o, src_o_check);
	std::cout << "Fermion expanded." << std::endl;
	std::cout << GridLogMessage << norm2(src_o) << " \t " << norm2(psi_o) << " \t " << norm2(src_o_check) << std::endl;

	LatticeGaugeField expandedUmu(eg.fine_gauge_grid);
	expand_gauge_field_qlat(eg, Umu, expandedUmu);
	std::cout << "Gauge field expanded." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << " \t " << WilsonLoops<PeriodicGimplR>::avgPlaquette(expandedUmu) << std::endl;

	MobiusFermionR expandedDMobius(expandedUmu, *eg.fine_fermion_grid, *eg.fine_fermion_rb_grid, *eg.fine_gauge_grid, *eg.fine_gauge_rb_grid, mass, M5, 22./12., 10./12.);
//	LatticeFermion src_o_dup(FrbGrid);
//	localConvert(src_o, src_o_dup);
//	std::cout << GridLogMessage << norm2(src_o) << " \t " << norm2(src_o_dup) << std::endl;

	GridStopWatch CGTimer;

	SchurDiagMooeeOperator<MobiusFermionR, LatticeFermion> HermOpEO(DMobius);
	SchurDiagMooeeOperator<MobiusFermionR, LatticeFermion> expandedHermOpEO(expandedDMobius);
	ConjugateGradient<LatticeFermion> CG(1.0e-10, 30000, 0);// switch off the assert

	LatticeFermion Mdag_src_o(eg.coarse_fermion_rb_grid);
	HermOpEO.AdjOp(src_o, Mdag_src_o);

	LatticeFermion x(eg.coarse_fermion_rb_grid);
	pickCheckerboard(Odd, x, src);
	LatticeFermion y(eg.coarse_fermion_rb_grid);
	LatticeFermion ydd(eg.coarse_fermion_rb_grid);

	LatticeFermion xd(eg.fine_fermion_rb_grid);
	LatticeFermion yd(eg.fine_fermion_rb_grid);

	expand_fermion(eg, x, xd);

	if(UGrid->ThisRank() != 0) x = zero;
	HermOpEO.Op(x, y);
	HermOpEO.AdjOp(y, x);
	
	expandedDMobius.ZeroCounters();
	WilsonFermion5DStatic::dirichlet = true;
	expandedHermOpEO.Op(xd, yd);
	expandedHermOpEO.AdjOp(yd, xd);
	WilsonFermion5DStatic::dirichlet = false;
	shrink_fermion(eg, xd, ydd);	
	
//	if(UGrid->ThisRank() != 0) x = zero;
//	if(UGrid->ThisRank() != 0) ydd = zero;

	std::cout << GridLogMessage << "|y - ydd|**2 : " << local_norm_sqr(x) << "\t" << local_norm_sqr(ydd) << std::endl;

	local_conjugate_gradient_MdagM(eg, expandedHermOpEO, xd, yd, 20);

	CGTimer.Start();
//	CG(HermOpEO, Mdag_src_o, result_o);
	CGTimer.Stop();

	std::cout << GridLogMessage << "Total CG time : " << CGTimer.Elapsed() << std::endl;

	std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
	DMobius.Report();
	expandedDMobius.Report();

	Grid_finalize();
}
