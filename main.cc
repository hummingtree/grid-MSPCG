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

int main(int argc, char** argv) {
	Grid_init(&argc, &argv);

	GridLogIterative.Active(1);

	const int Ls = 12;

	GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
//	GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
//	GridCartesian* FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
//	GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);
//
	UGrid->show_decomposition();

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
//	std::string file("/global/homes/j/jiquntu/configurations/64x128x10I_mh0.02659_ml0.000661/ckpoint_lat.2850");
	NerscIO::readConfiguration(Umu, header, file);

//	SU3::HotConfiguration(RNG4, Umu);

	std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << "   Ls: " << Ls << std::endl;

	RealD mass = 0.0001;
	RealD fD_mass = mass;
	RealD M5 = 1.8;
//	RealD b = 1.5;
//	RealD c = 0.5;
	RealD b = 22./12.;
	RealD c = 10./12.;
	MobiusFermionR DMobius(Umu, *eg.coarse_fermion_grid, *eg.coarse_fermion_rb_grid, *eg.coarse_gauge_grid, *eg.coarse_gauge_rb_grid, mass, M5, b, c);
	DMobius.ZeroCounters();
	
	LatticeFermionD src_o(eg.coarse_fermion_rb_grid);
	LatticeFermionD result_o(eg.coarse_fermion_rb_grid);
	pickCheckerboard(Odd, src_o, src);
	result_o = zero;


	LatticeFermionD psi_o(eg.fine_fermion_rb_grid);
	expand_fermion_qlat(eg, src_o, psi_o);
	
	LatticeFermionD src_o_check(eg.coarse_fermion_rb_grid);
	shrink_fermion(eg, psi_o, src_o_check);
	std::cout << "Fermion expanded." << std::endl;
	std::cout << GridLogMessage << norm2(src_o) << " \t " << local_norm_sqr_center(eg, psi_o) << " \t " << norm2(src_o_check) << std::endl;

	LatticeGaugeFieldD expandedUmu(eg.fine_gauge_grid);
	expand_gauge_field_qlat(eg, Umu, expandedUmu);
	std::cout << "Gauge field expanded." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << " \t " << WilsonLoops<PeriodicGimplR>::avgPlaquette(expandedUmu) << std::endl;

	LatticeGaugeFieldF expanded_gf_F(eg.fine_gauge_grid_F);
	precisionChange(expanded_gf_F, expandedUmu);
	std::cout << "Changed to single precision Gauge field." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplF>::avgPlaquette(expanded_gf_F) << std::endl;

//	MobiusFermionD expandedDMobius(expandedUmu, *eg.fine_fermion_grid, *eg.fine_fermion_rb_grid, *eg.fine_gauge_grid, *eg.fine_gauge_rb_grid, fD_mass, M5, b, c);
//	expandedDMobius.StencilOdd.zero_comm_recv_buffers();
//	expandedDMobius.StencilEven.zero_comm_recv_buffers();
//	LatticeFermion src_o_dup(FrbGrid);
//	localConvert(src_o, src_o_dup);
//	std::cout << GridLogMessage << norm2(src_o) << " \t " << norm2(src_o_dup) << std::endl;

	MobiusFermionF fDF(expanded_gf_F, *eg.fine_fermion_grid_F, *eg.fine_fermion_rb_grid_F, *eg.fine_gauge_grid_F, *eg.fine_gauge_rb_grid_F, fD_mass, M5, b, c);
	fDF.ZeroCounters();
	
	GridStopWatch CGTimer;
	GridStopWatch MSPCGTimer;

	SchurDiagMooeeOperator<MobiusFermionD, LatticeFermionD> HermOpEO(DMobius);
//	SchurDiagMooeeOperator<MobiusFermionD, LatticeFermionD> expandedHermOpEO(expandedDMobius);
	SchurDiagMooeeOperator<MobiusFermionF, LatticeFermionF> fHermF(fDF);
//	ConjugateGradient<LatticeFermion> CG(1.0e-8, 30000, 0);// switch off the assert

	LatticeFermionD Mdag_src_o(eg.coarse_fermion_rb_grid);
	HermOpEO.AdjOp(src_o, Mdag_src_o);

	Mdag_src_o.checkerboard = Odd;

	LatticeFermionD cXD(eg.coarse_fermion_rb_grid);
	pickCheckerboard(Odd, cXD, src);
	LatticeFermionD cYD(eg.coarse_fermion_rb_grid);
	LatticeFermionD ydd(eg.coarse_fermion_rb_grid);

	cYD.checkerboard = Odd;
	cYD = zero;

	LatticeFermionD fXD(eg.fine_fermion_rb_grid);
	LatticeFermionD fYD(eg.fine_fermion_rb_grid);
	fXD.checkerboard = Odd;
	fYD.checkerboard = Odd;

	LatticeFermionF fXF(eg.fine_fermion_rb_grid_F);
	LatticeFermionF fYF(eg.fine_fermion_rb_grid_F);
	fXF.checkerboard = Odd;
	fYF.checkerboard = Odd;

	expand_fermion_qlat(eg, cXD, fXD);
	precisionChange(fXF, fXD);

	if(UGrid->ThisRank() != 0) cXD = zero;
	HermOpEO.Op(cXD, cYD);
	HermOpEO.AdjOp(cYD, cXD);

//	expandedDMobius.ZeroCounters();
//	WilsonFermion5DStatic::dirichlet = true;
	fHermF.Op(fXF, fYF);
	fHermF.AdjOp(fYF, fXF);
//	WilsonFermion5DStatic::dirichlet = false;
	
	precisionChange(fXD, fXF);
	shrink_fermion(eg, fXD, cYD);	
	
//	if(UGrid->ThisRank() != 0) x = zero;
//	if(UGrid->ThisRank() != 0) ydd = zero;

	std::cout << GridLogMessage << "|cYD - cXD|**2 : " << local_norm_sqr(cXD) << "\t" << local_norm_sqr(cYD) << std::endl;

//	local_conjugate_gradient_MdagM_variant(eg, expandedHermOpEO, xd, yd, 20);

//	CGTimer.Start();
//	CG(HermOpEO, Mdag_src_o, result_o);
//	CGTimer.Stop();
//	std::cout << GridLogMessage << "Total CG time : " << CGTimer.Elapsed() << std::endl;

	int local_iter = 6;
	RealD local_e  = 0.;
	std::cout << GridLogMessage << "MSPCG local iteration : " << local_iter << std::endl;
	std::cout << GridLogMessage << "MSPCG local mass      : " << fD_mass << std::endl;
	std::cout << GridLogMessage << "MSPCG local e         : " << local_e << std::endl;
	
	MSPCGTimer.Start();
//	MSP_conjugate_gradient(eg, HermOpEO, expandedHermOpEO, Mdag_src_o, y, 1e-7, local_iter, 50000, local_e);
//	MSPCG_half(eg, HermOpEO, fHermF, Mdag_src_o, y, 1e-10, local_iter, 50000);
//	DD_CG(eg, HermOpEO, expandedHermOpEO, Mdag_src_o, y, 1e-7, local_iter, 50000);
	MSPCGTimer.Stop();
	std::cout << GridLogMessage << "Total MSPCG time : " << MSPCGTimer.Elapsed() << std::endl;

	std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
//	DMobius.Report();
	DMobius.CayleyReport();
//	expandedDMobius.Report();
	fDF.CayleyReport();

	Timer::display();

	Grid_finalize();
}
