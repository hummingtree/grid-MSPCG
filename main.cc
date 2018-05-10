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

	std::array<int, 4> padding = {2,2,2,2};
	std::array<int, 4> inner_padding = {0,0,0,0};
	ExpandGrid eg; eg.init(UGrid, padding, inner_padding, Ls);
	
	std::vector<int> seeds4({1, 2, 3, 4});
	std::vector<int> seeds5({5, 6, 7, 8});
	GridParallelRNG RNG5(eg.cFGrid);
	RNG5.SeedFixedIntegers(seeds5);
	GridParallelRNG RNG4(UGrid);
	RNG4.SeedFixedIntegers(seeds4);

	LatticeFermion src(eg.cFGrid);
	random(RNG5, src);
	LatticeFermion result(eg.cFGrid);
	result = zero;
	LatticeGaugeField Umu(UGrid);
	LatticeGaugeField shifted_Umu(UGrid);
	
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
	shifted_Umu = Cshift(        Umu, 0, UGrid->_ldimensions[0]/2);
	shifted_Umu = Cshift(shifted_Umu, 1, UGrid->_ldimensions[1]/2);
	shifted_Umu = Cshift(shifted_Umu, 2, UGrid->_ldimensions[2]/2);
	shifted_Umu = Cshift(shifted_Umu, 3, UGrid->_ldimensions[3]/2);

//	SU3::HotConfiguration(RNG4, Umu);

	std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << "   Ls: " << Ls << std::endl;

	RealD mass = 0.0001;
	RealD fD_mass = mass;
	RealD M5 = 1.8;
//	RealD b = 1.5;
//	RealD c = 0.5;
	RealD b = 22./12.;
	RealD c = 10./12.;
	MobiusFermionR DMobius(Umu, *eg.cFGrid, *eg.cFrbGrid, *eg.cUGrid, *eg.cUrbGrid, mass, M5, b, c);
	DMobius.ZeroCounters();

	MobiusFermionR shifted_DMobius(shifted_Umu, *eg.cFGrid, *eg.cFrbGrid, *eg.cUGrid, *eg.cUrbGrid, mass, M5, b, c);
	shifted_DMobius.ZeroCounters();

	LatticeFermionD src_o(eg.cFrbGrid);
	LatticeFermionD result_o(eg.cFrbGrid);
	pickCheckerboard(Odd, src_o, src);
	result_o = zero;


	LatticeFermionD psi_o(eg.fFrbGrid);
	expand_fermion(eg, src_o, psi_o);
	
	LatticeFermionD src_o_check(eg.cFrbGrid);
	shrink_fermion(eg, psi_o, src_o_check);
	std::cout << "Fermion expanded." << std::endl;
	std::cout << GridLogMessage << norm2(src_o) << " \t " << norm2(psi_o) << " \t " << norm2(src_o_check) << std::endl;

	LatticeGaugeFieldD fUField(eg.fUGrid);
	expand_gauge_field_qlat(eg, Umu, fUField);
	std::cout << "Gauge field expanded." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << " \t " << WilsonLoops<PeriodicGimplR>::avgPlaquette(fUField) << std::endl;

	LatticeGaugeFieldF fUField_F(eg.fUGrid_F);
	precisionChange(fUField_F, fUField);
	std::cout << "Changed to single precision Gauge field." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplF>::avgPlaquette(fUField_F) << std::endl;
	
	LatticeGaugeFieldD shifted_fUField(eg.fUGrid);
	expand_gauge_field_qlat(eg, shifted_Umu, shifted_fUField);
	std::cout << "shifted Gauge field expanded." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplR>::avgPlaquette(shifted_Umu) << " \t " << WilsonLoops<PeriodicGimplR>::avgPlaquette(shifted_fUField) << std::endl;

	LatticeGaugeFieldF shifted_fUField_F(eg.fUGrid_F);
	precisionChange(shifted_fUField_F, shifted_fUField);
	std::cout << "shifted Changed to single precision Gauge field." << std::endl;
	std::cout << GridLogMessage << WilsonLoops<PeriodicGimplF>::avgPlaquette(shifted_fUField_F) << std::endl;

//	MobiusFermionD expandedDMobius(expandedUmu, *eg.fFGrid, *eg.fFrbGrid, *eg.fUGrid, *eg.fUrbGrid, fD_mass, M5, b, c);
//	expandedDMobius.StencilOdd.zero_comm_recv_buffers();
//	expandedDMobius.StencilEven.zero_comm_recv_buffers();
//	LatticeFermion src_o_dup(FrbGrid);
//	localConvert(src_o, src_o_dup);
//	std::cout << GridLogMessage << norm2(src_o) << " \t " << norm2(src_o_dup) << std::endl;

	MobiusFermionF fDF(fUField_F, *eg.fFGrid_F, *eg.fFrbGrid_F, *eg.fUGrid_F, *eg.fUrbGrid_F, fD_mass, M5, b, c);
	fDF.ZeroCounters();
	MobiusFermionF shifted_fDF(shifted_fUField_F, *eg.fFGrid_F, *eg.fFrbGrid_F, *eg.fUGrid_F, *eg.fUrbGrid_F, fD_mass, M5, b, c);
	shifted_fDF.ZeroCounters();
	
	GridStopWatch CGTimer;
	GridStopWatch MSPCGTimer;

	SchurDiagMooeeOperator<MobiusFermionD, LatticeFermionD> HermOpEO(DMobius);
	SchurDiagMooeeOperator<MobiusFermionD, LatticeFermionD> shifted_HermOpEO(shifted_DMobius);
//	SchurDiagMooeeOperator<MobiusFermionD, LatticeFermionD> expandedHermOpEO(expandedDMobius);
	SchurDiagMooeeOperator<MobiusFermionF, LatticeFermionF> fHermF(fDF);
	SchurDiagMooeeOperator<MobiusFermionF, LatticeFermionF> shifted_fHermF(shifted_fDF);
//	ConjugateGradient<LatticeFermion> CG(1.0e-8, 30000, 0);// switch off the assert

	LatticeFermionD Mdag_src_o(eg.cFrbGrid);
	HermOpEO.AdjOp(src_o, Mdag_src_o);

	Mdag_src_o.checkerboard = Odd;

	LatticeFermionD cXD(eg.cFrbGrid);
	pickCheckerboard(Odd, cXD, src);
	LatticeFermionD cYD(eg.cFrbGrid);
	LatticeFermionD ydd(eg.cFrbGrid);
	
	LatticeFermionD shifted_cXD(eg.cFrbGrid);
	pickCheckerboard(Odd, cXD, src);
	LatticeFermionD shifted_cYD(eg.cFrbGrid);
	shifted_cXD.checkerboard = Odd;
	shifted_cYD.checkerboard = Odd;

	shift_ff(cXD, shifted_cXD, true);

	cYD.checkerboard = Odd;
	cYD = zero;

	LatticeFermionD fXD(eg.fFrbGrid);
	LatticeFermionD fYD(eg.fFrbGrid);
	fXD.checkerboard = Odd;
	fYD.checkerboard = Odd;

	LatticeFermionF fXF(eg.fFrbGrid_F);
	LatticeFermionF fYF(eg.fFrbGrid_F);
	fXF.checkerboard = Odd;
	fYF.checkerboard = Odd;

	printf("Right before fermion expansion with qlat.\n");
	
	expand_fermion_D2F_qlat(eg, cXD, fXF);
//	expand_fermion_D2F(eg, cXD, fXF);
	zero_boundary_fermion_inner(eg, fXF);
//	precisionChange(fXF, fXD);

//	if(UGrid->ThisRank() != 0) cXD = zero;
	HermOpEO.Op(cXD, cYD);
	HermOpEO.AdjOp(cYD, cXD);

	shifted_HermOpEO.Op(shifted_cXD, shifted_cYD);
	shifted_HermOpEO.AdjOp(shifted_cYD, shifted_cXD);
	shift_ff(shifted_cXD, shifted_cYD, false);

//	expandedDMobius.ZeroCounters();
//	WilsonFermion5DStatic::dirichlet = true;
	fHermF.Op(fXF, fYF);
	fHermF.AdjOp(fYF, fXF);
//	WilsonFermion5DStatic::dirichlet = false;
	
//	precisionChange(fXD, fXF);
	shrink_fermion_F2D(eg, fXF, cYD);	
	
//	if(UGrid->ThisRank() != 0) x = zero;
//	if(UGrid->ThisRank() != 0) ydd = zero;

//	std::cout << GridLogMessage << "|cYD - cXD|**2 : " << local_norm_sqr(cXD) << "\t" << local_norm_sqr(cYD) << std::endl;
//	if(UGrid->IsBoss()){
//		printf("|cXD|**2 = %20.16e, |cYD|**2 = %20.16e\n", local_norm_sqr(cXD), local_norm_sqr(cYD));
		printf("|cXD|**2 = %20.16e, |cYD|**2 = %20.16e, |shifted_cYD|**2 = %20.16e\n", norm2(cXD), norm2(cYD), norm2(shifted_cYD));
//	}

//	local_conjugate_gradient_MdagM_variant(eg, expandedHermOpEO, xd, yd, 20);

//	CGTimer.Start();
//	CG(HermOpEO, Mdag_src_o, result_o);
//	CGTimer.Stop();
//	std::cout << GridLogMessage << "Total CG time : " << CGTimer.Elapsed() << std::endl;

	cYD = zero;

	int local_iter = 6;
	RealD local_e  = 0.;
	std::cout << GridLogMessage << "MSPCG local iteration : " << local_iter << std::endl;
	std::cout << GridLogMessage << "MSPCG local mass      : " << fD_mass << std::endl;
	std::cout << GridLogMessage << "MSPCG local e         : " << local_e << std::endl;
	std::cout << GridLogMessage << "MSPCG local padding   : [" << padding[0] << " " << padding[1] << " " << padding[2] << " " << padding[3] << "]"<< std::endl;
	std::cout << GridLogMessage << "MSPCG inner padding   : [" << inner_padding[0] << " " << inner_padding[1] << " " << inner_padding[2] << " " << inner_padding[3] << "]"<< std::endl;
	
	MSPCGTimer.Start();
//	MSP_conjugate_gradient(eg, HermOpEO, expandedHermOpEO, Mdag_src_o, y, 1e-7, local_iter, 50000, local_e);
//	MSPCG_half(eg, HermOpEO, fHermF, Mdag_src_o, cYD, 1e-10, local_iter, 50000, local_e);
	MSPCG_half(eg, HermOpEO, fHermF, Mdag_src_o, cYD, 1e-10, local_iter, 50000, local_e);
//	MSPCG_shift(eg, HermOpEO, fHermF, shifted_HermOpEO, shifted_fHermF, Mdag_src_o, cYD, 1e-10, local_iter, 50000, local_e);
//	MSPCG_pad(eg, HermOpEO, fHermF, Mdag_src_o, cYD, 1e-10, local_iter, 50000, local_e);
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
