#ifndef TEST_OPERATOR_ELASTIC_H
#define TEST_OPERATOR_ELASTIC_H

#include <AMReX_Vector.H>

namespace Test
{
namespace Operator
{
class Elastic
{
public:
	Elastic() {};

	/// Create a simulation domain where the coarse level has ncells x ncells ( x ncells)
	/// and nlevels.
	/// Each level is identical to the one below but with a quarter of the area and twice
	/// the resolution.
	/// \todo Add option to specify one of multiple refinement geometries
	void Define(int _ncells,
		    int _nlevels);

	/// Impose a set of prescribed values for \f$\mathbf{u}\f$.
	/// Compute right hand side, then use Operator::Reflux to modify the residual \f$\mathbf{r}=\mathbf{b}-D\mathbf{u}\f$.
	/// Residual should always be zero unless something goes wrong with Reflux.
	/// If this test fails, check Reflux.
	int RefluxTest(int);

	/// Compute the exact solution of the governing equation
	///   \f[C_{ijkl} u_{k,jl} + b_i = 0\f]
	/// Where
	///   \f[C_{ijkl} = \alpha \delta_{ik} \delta_{jl}\f]
	/// (see Model::Solid::LinearElastic::Laplacian for implementation of this model)
	///
	/// Consider a simulation domain \f$\Omega = [0, 1]^2\f$
	///
	/// Expand \f$\mathbf{u}\f$ and \f$\mathbf{b}\f$ out as trigonometric series
	///       \f[ u_i = u_i^{mn} \sin(\pi\ x\ m) \sin(\pi\ y\ n) \f] 
	///       \f[ b_i = b_i^{mn} \sin(\pi\ x\ m) \sin(\pi\ y\ n) \f]
	///
	/// Then the second derivative of displacement is
	///       \f[ u_{i,jj} = - u_i^{mn} \pi^2 (m^2 + n^2) \sin(\pi\ x\ m) \sin(\pi\ y\ n) \f]
	/// and the pointwise force is 
	///       \f[ C_{ijkl} u_{k,jl} = \alpha \delta_{ik} \delta_{jl} u_{k,jl} = \alpha u_{i,jj} \f]
	///                     
	/// Substituting into the governing equation gives
	///       \f[-\alpha u_i^{mn} \pi^2 (m^2 + n^2) \sin(\pi\ x\ m) \sin(\pi\ y\ n) 
	///          + b_i^{mn} \sin(\pi\ x\ m)\sin(\pi\ y\ n) 
	///                     = 0\f]
	/// or
	///     \f[\alpha u_i^{mn} \pi^2 (m^2 + n^2) \sin(\pi\ x\ m)\sin(\pi y n) = b_i^{mn} \sin(\pi\ x\ m)\sin(\pi\ y\ n)\f]
	///
	/// Using orthognonality of the eigenfunctions gives
	///     \f[\alpha u_i^{mn} \pi^2 (m^2 + n^2) = b_i^{mn}\f]
	/// or, simplifying
	///     \f[u_i^{mn}  = \frac{b_i^{mn}}{(\alpha\ \pi^2\ (m^2 + n^2))}\f]
	/// 
	///
	int TrigTest(bool verbose,
		     int component = 0,
		     int n = 1,
		     std::string plotfile = "");

	int SPD(bool verbose,
		std::string plotfile = "");

	int SpatiallyVaryingCTest(int);
	
private: // Private member functions

	void WritePlotFile(std::string plotfile)
	{
		const int output_comp = varname.size();
		amrex::Vector<amrex::MultiFab> plotmf(nlevels);
		for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			//if (ilev==1) ngrids[ilev].growHi(1,1);
			plotmf			[ilev].define(ngrids[ilev], dmap[ilev], output_comp, 0);
			amrex::MultiFab::Copy(plotmf	[ilev], solution_exact [ilev], 0, 0*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], solution_numeric[ilev],0, 1*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], solution_error [ilev], 0, 2*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], rhs_prescribed [ilev], 0, 3*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], rhs_exact      [ilev], 0, 4*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], rhs_numeric    [ilev], 0, 5*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], res_exact      [ilev], 0, 6*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], res_numeric    [ilev], 0, 7*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
			amrex::MultiFab::Copy(plotmf	[ilev], ghost_force    [ilev], 0, 8*AMREX_SPACEDIM, AMREX_SPACEDIM, 0); // 
		}


		amrex::WriteMultiLevelPlotfile(plotfile, nlevels, amrex::GetVecOfConstPtrs(plotmf),
					       varname, geom, 0.0, amrex::Vector<int>(nlevels, 0),
					       amrex::Vector<amrex::IntVect>(nlevels, amrex::IntVect{ref_ratio}));
	}		

	
private: // Private member variables
	int nlevels = 1;
	int ncells = 16;
	const int ref_ratio = 2;
 	amrex::Vector<amrex::Geometry> 		   geom;
 	amrex::Vector<amrex::BoxArray> 		   cgrids, ngrids;
 	amrex::Vector<amrex::DistributionMapping>  dmap;
	
 	amrex::Vector<amrex::MultiFab>             solution_exact;
 	amrex::Vector<amrex::MultiFab>             solution_numeric;
 	amrex::Vector<amrex::MultiFab>             solution_error;

 	amrex::Vector<amrex::MultiFab>             rhs_prescribed;
 	amrex::Vector<amrex::MultiFab>             rhs_exact;
 	amrex::Vector<amrex::MultiFab>             rhs_numeric;

 	amrex::Vector<amrex::MultiFab>             res_exact;
 	amrex::Vector<amrex::MultiFab>             res_numeric;
 	amrex::Vector<amrex::MultiFab>             ghost_force;

	amrex::Vector<std::string> varname = {AMREX_D_DECL("solution_exact1", "solution_exact2","solution_exact3"),
					      AMREX_D_DECL("solution_numeric1", "solution_numeric2","solution_numeric3"),
					      AMREX_D_DECL("solution_error1", "solution_error2","solution_error3"),
					      AMREX_D_DECL("rhs_prescribed1", "rhs_prescribed2","rhs_prescribed3"),
					      AMREX_D_DECL("rhs_exact1","rhs_exact2","rhs_exact2"),
					      AMREX_D_DECL("rhs_numeric1","rhs_numeric2","rhs_numeric3"),
					      AMREX_D_DECL("res_exact1","res_exact2","res_exact3"),
					      AMREX_D_DECL("res_numeric1","res_numeric2","res_numeric3"),
					      AMREX_D_DECL("ghost_force1", "ghost_force2","ghost_force2")};
};
}
}





#endif

