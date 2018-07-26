
#include "Operator2.H"
#include <AMReX_MLLinOp_F.H>
#include "AMReX_MG_F.H"
#include <AMReX_MultiFabUtil.H>

namespace Operator {

Operator2::Operator2 ()
{
	m_ixtype = amrex::IntVect::TheCellVector();
}

Operator2::~Operator2 () {}

void
Operator2::define (const amrex::Vector<amrex::Geometry>& a_geom,
		   const amrex::Vector<amrex::BoxArray>& a_grids,
		   const amrex::Vector<amrex::DistributionMapping>& a_dmap,
		   BC::BC& a_bc,
		   const amrex::LPInfo& a_info,
		   const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory)
{
	MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
	defineAuxData();
	defineBC();
	m_bc = &a_bc;
}

void
Operator2::defineAuxData ()
{
	BL_PROFILE("Operator2::defineAuxData()");

	m_undrrelxr.resize(m_num_amr_levels);
	m_maskvals.resize(m_num_amr_levels);
	m_fluxreg.resize(m_num_amr_levels-1);

	const int ncomp = getNComp();

	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_undrrelxr[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			m_undrrelxr[amrlev][mglev].define(m_grids[amrlev][mglev],
							  m_dmap[amrlev][mglev],
							  1, 0, 0, ncomp);
		}
	}
    
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_maskvals[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			for (amrex::OrientationIter oitr; oitr; ++oitr)
			{
				const amrex::Orientation face = oitr();
				const int ngrow = 1;
				const int extent = 1; // extend to corners
				m_maskvals[amrlev][mglev][face].define(m_grids[amrlev][mglev],
								       m_dmap[amrlev][mglev],
								       m_geom[amrlev][mglev],
								       face, 0, ngrow, 0, ncomp, true);
			}
		}
	}

	for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
	{
		const amrex::IntVect ratio{m_amr_ref_ratio[amrlev]};
		m_fluxreg[amrlev].define(m_grids[amrlev+1][0], m_grids[amrlev][0],
					 m_dmap[amrlev+1][0], m_dmap[amrlev][0],
					 m_geom[amrlev+1][0], m_geom[amrlev][0],
					 ratio, amrlev+1, ncomp);
	}
}

void
Operator2::defineBC ()
{
	BL_PROFILE("Operator2::defineBC()");

	const int ncomp = getNComp();

	m_bndry_sol.resize(m_num_amr_levels);
	m_crse_sol_br.resize(m_num_amr_levels);

	m_bndry_cor.resize(m_num_amr_levels);
	m_crse_cor_br.resize(m_num_amr_levels);

	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_bndry_sol[amrlev].reset(new amrex::MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
							       ncomp, m_geom[amrlev][0]));
	}

	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
	{
		const int in_rad = 0;
		const int out_rad = 1;
		const int extent_rad = 2;
		const int crse_ratio = m_amr_ref_ratio[amrlev-1];
		amrex::BoxArray cba = m_grids[amrlev][0];
		cba.coarsen(crse_ratio);
		m_crse_sol_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
								     in_rad, out_rad, extent_rad, ncomp));
	}

	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
	{
		const int in_rad = 0;
		const int out_rad = 1;
		const int extent_rad = 2;
		const int crse_ratio = m_amr_ref_ratio[amrlev-1];
		amrex::BoxArray cba = m_grids[amrlev][0];
		cba.coarsen(crse_ratio);
		m_crse_cor_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
								     in_rad, out_rad, extent_rad, ncomp));
		m_crse_cor_br[amrlev]->setVal(0.0);
	}

	// This has be to done after m_crse_cor_br is defined.
	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_bndry_cor[amrlev].reset(new amrex::MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
							       ncomp, m_geom[amrlev][0]));
		amrex::MultiFab bc_data(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp, 1);
		bc_data.setVal(0.0);

		m_bndry_cor[amrlev]->setBndryValues(*m_crse_cor_br[amrlev], 0, bc_data, 0, 0, ncomp,
						    m_amr_ref_ratio[amrlev-1], amrex::BCRec());

		m_bndry_cor[amrlev]->setLOBndryConds({AMREX_D_DECL(BCType::Dirichlet,
								   BCType::Dirichlet,
								   BCType::Dirichlet)},
			{AMREX_D_DECL(BCType::Dirichlet,
				      BCType::Dirichlet,
				      BCType::Dirichlet)},
			m_amr_ref_ratio[amrlev-1], amrex::RealVect{});
	}

	m_bcondloc.resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_bcondloc[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			m_bcondloc[amrlev][mglev].reset(new BndryCondLoc(m_grids[amrlev][mglev],
										m_dmap[amrlev][mglev]));
		} 
	}
}

void
Operator2::setLevelBC (int amrlev, const amrex::MultiFab* a_levelbcdata)
{
	BL_PROFILE("Operator2::setLevelBC()");

	AMREX_ALWAYS_ASSERT(amrlev >= 0 && amrlev < m_num_amr_levels);

	const int ncomp = getNComp();

	amrex::MultiFab zero;
	if (a_levelbcdata == nullptr) {
		zero.define(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp, 1);
		zero.setVal(0.0);
	} else {
		AMREX_ALWAYS_ASSERT(a_levelbcdata->nGrow() >= 1);
	}
	const amrex::MultiFab& bcdata = (a_levelbcdata == nullptr) ? zero : *a_levelbcdata;

	int br_ref_ratio = -1;

	if (amrlev == 0)
	{
		if (needsCoarseDataForBC())
		{
			br_ref_ratio = m_coarse_data_crse_ratio > 0 ? m_coarse_data_crse_ratio : 2;
			if (m_crse_sol_br[amrlev] == nullptr && br_ref_ratio > 0)
			{
				const int in_rad = 0;
				const int out_rad = 1;
				const int extent_rad = 2;
				const int crse_ratio = br_ref_ratio;
				amrex::BoxArray cba = m_grids[amrlev][0];
				cba.coarsen(crse_ratio);
				m_crse_sol_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
										     in_rad, out_rad,
										     extent_rad, ncomp));
			}
			if (m_coarse_data_for_bc != nullptr) {
				AMREX_ALWAYS_ASSERT(m_coarse_data_crse_ratio > 0);
				const amrex::Box& cbx = amrex::coarsen(m_geom[0][0].Domain(), m_coarse_data_crse_ratio);
				m_crse_sol_br[amrlev]->copyFrom(*m_coarse_data_for_bc, 0, 0, 0, ncomp,
								amrex::Geometry::periodicity(cbx));
			} else {
				m_crse_sol_br[amrlev]->setVal(0.0);
			}
			m_bndry_sol[amrlev]->setBndryValues(*m_crse_sol_br[amrlev], 0,
							    bcdata, 0, 0, ncomp,
							    br_ref_ratio, amrex::BCRec());
			br_ref_ratio = m_coarse_data_crse_ratio;
		}
		else
		{
			m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,ncomp,amrex::BCRec());
			br_ref_ratio = 1;
		}
	}
	else
	{
		m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,ncomp, m_amr_ref_ratio[amrlev-1], amrex::BCRec());
		br_ref_ratio = m_amr_ref_ratio[amrlev-1];
	}

	m_bndry_sol[amrlev]->setLOBndryConds(m_lobc, m_hibc, br_ref_ratio, m_coarse_bc_loc);

	const amrex::Real* dx = m_geom[amrlev][0].CellSize();
	for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
	{
		m_bcondloc[amrlev][mglev]->setLOBndryConds(m_geom[amrlev][mglev], dx,
							   m_lobc, m_hibc,
							   br_ref_ratio, m_coarse_bc_loc);
	}
}

amrex::BoxArray
Operator2::makeNGrids (int grid_size) const
{
	const amrex::Box& dombx = m_geom[0].back().Domain();

	const amrex::BoxArray& old_ba = m_grids[0].back();
	const int N = old_ba.size();
	amrex::Vector<amrex::Box> bv;
	bv.reserve(N);
	for (int i = 0; i < N; ++i)
	{
		amrex::Box b = old_ba[i];
		b.coarsen(grid_size);
		b.refine(grid_size);
		amrex::IntVect sz = b.size();
		const amrex::IntVect nblks {AMREX_D_DECL(sz[0]/grid_size, sz[1]/grid_size, sz[2]/grid_size)};
        
		amrex::IntVect big = b.smallEnd() + grid_size - 1;
		b.setBig(big);

#if (AMREX_SPACEDIM == 3)
		for (int kk = 0; kk < nblks[2]; ++kk) {
#endif
#if (AMREX_SPACEDIM >= 2)
			for (int jj = 0; jj < nblks[1]; ++jj) {
#endif
				for (int ii = 0; ii < nblks[0]; ++ii)
				{
					amrex::IntVect shft{AMREX_D_DECL(ii*grid_size,jj*grid_size,kk*grid_size)};
					amrex::Box bb = amrex::shift(b,shft);
					bb &= dombx;
					bv.push_back(bb);
				}
#if (AMREX_SPACEDIM >= 2)
			}
#endif
#if (AMREX_SPACEDIM == 3)
		}
#endif
	}

	std::sort(bv.begin(), bv.end());
	bv.erase(std::unique(bv.begin(), bv.end()), bv.end());

	amrex::BoxList bl(std::move(bv));

	return amrex::BoxArray{std::move(bl)};
}

void
Operator2::restriction (int, int, amrex::MultiFab& crse, amrex::MultiFab& fine) const
{
	const int ncomp = getNComp();
	amrex::average_down(fine, crse, 0, ncomp, 2);
}

void
Operator2::interpolation (int amrlev, int fmglev,
			  amrex::MultiFab& fine, const amrex::MultiFab& crse) const
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(crse,true); mfi.isValid(); ++mfi)
	{
		const amrex::Box&         bx    = mfi.tilebox();
		const int          ncomp = getNComp();
		const amrex::FArrayBox& cfab    = crse[mfi];
		amrex::FArrayBox&       ffab    = fine[mfi];

		amrex_mg_interp(ffab.dataPtr(),
				AMREX_ARLIM(ffab.loVect()), AMREX_ARLIM(ffab.hiVect()),
				cfab.dataPtr(),
				AMREX_ARLIM(cfab.loVect()), AMREX_ARLIM(cfab.hiVect()),
				bx.loVect(), bx.hiVect(), &ncomp);
	}    
}

void
Operator2::averageDownSolutionRHS (int camrlev, amrex::MultiFab& crse_sol, amrex::MultiFab& crse_rhs,
				   const amrex::MultiFab& fine_sol, const amrex::MultiFab& fine_rhs)
{
	const auto amrrr = AMRRefRatio(camrlev);
	const int ncomp = getNComp();
	amrex::average_down(fine_sol, crse_sol, 0, ncomp, amrrr);
	amrex::average_down(fine_rhs, crse_rhs, 0, ncomp, amrrr);
}

void
Operator2::apply (int amrlev, int mglev, amrex::MultiFab& out, amrex::MultiFab& in, BCMode bc_mode,
		  const amrex::MLMGBndry* bndry) const
{
	BL_PROFILE("Operator2::apply()");
	applyBC(amrlev, mglev, in, bc_mode, bndry);
	Fapply(amrlev, mglev, out, in);
}

void
Operator2::smooth (int amrlev, int mglev, amrex::MultiFab& sol, const amrex::MultiFab& rhs,
		   bool skip_fillboundary) const
{
	BL_PROFILE("Operator2::smooth()");
	for (int redblack = 0; redblack < 2; ++redblack)
	{
		applyBC(amrlev, mglev, sol, BCMode::Homogeneous, nullptr, skip_fillboundary);
		Fsmooth(amrlev, mglev, sol, rhs, redblack);
		skip_fillboundary = false;
	}
}

void
Operator2::updateSolBC (int amrlev, const amrex::MultiFab& crse_bcdata) const
{
	BL_PROFILE("Operator2::updateSolBC()");

	AMREX_ALWAYS_ASSERT(amrlev > 0);
	const int ncomp = getNComp();
	m_crse_sol_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_geom[amrlev-1][0].periodicity());
	m_bndry_sol[amrlev]->updateBndryValues(*m_crse_sol_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
}

void
Operator2::updateCorBC (int amrlev, const amrex::MultiFab& crse_bcdata) const
{
	BL_PROFILE("Operator2::updateCorBC()");
	AMREX_ALWAYS_ASSERT(amrlev > 0);
	const int ncomp = getNComp();
	m_crse_cor_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_geom[amrlev-1][0].periodicity());
	m_bndry_cor[amrlev]->updateBndryValues(*m_crse_cor_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
}

void
Operator2::solutionResidual (int amrlev, amrex::MultiFab& resid, amrex::MultiFab& x,
			     const amrex::MultiFab& b, const amrex::MultiFab* crse_bcdata)
{
	BL_PROFILE("Operator2::solutionResidual()");
	const int ncomp = getNComp();
	if (crse_bcdata != nullptr) {
		updateSolBC(amrlev, *crse_bcdata);
	}
	const int mglev = 0;
	apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

	AMREX_ALWAYS_ASSERT(resid.nComp() == b.nComp());
	amrex::MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 0);
}

void
Operator2::fillSolutionBC (int amrlev, amrex::MultiFab& sol, const amrex::MultiFab* crse_bcdata)
{
	BL_PROFILE("Operator2::fillSolutionBC()");
	if (crse_bcdata != nullptr) {
		updateSolBC(amrlev, *crse_bcdata);
	}
	const int mglev = 0;
	applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());    
}

void
Operator2::correctionResidual (int amrlev, int mglev, amrex::MultiFab& resid, amrex::MultiFab& x,
			       const amrex::MultiFab& b,
			       BCMode bc_mode, const amrex::MultiFab* crse_bcdata)
{
	BL_PROFILE("Operator2::correctionResidual()");
	const int ncomp = getNComp();
	if (bc_mode == BCMode::Inhomogeneous)
	{
		if (crse_bcdata)
		{
			AMREX_ALWAYS_ASSERT(mglev == 0);
			AMREX_ALWAYS_ASSERT(amrlev > 0);
			updateCorBC(amrlev, *crse_bcdata);
		}
		apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, m_bndry_cor[amrlev].get());
	}
	else
	{
		AMREX_ALWAYS_ASSERT(crse_bcdata == nullptr);
		apply(amrlev, mglev, resid, x, BCMode::Homogeneous, nullptr);
	}

	amrex::MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 0);
}

void
Operator2::applyBC (int amrlev, int mglev, amrex::MultiFab& in, BCMode bc_mode,
		    const amrex::MLMGBndry* bndry, bool skip_fillboundary) const
{
	BL_PROFILE("Operator2::applyBC()");
	// No coarsened boundary values, cannot apply inhomog at mglev>0.
	BL_ASSERT(mglev == 0 || bc_mode == BCMode::Homogeneous);
	BL_ASSERT(bndry != nullptr || bc_mode == BCMode::Homogeneous);

	const int ncomp = getNComp();
	const int cross = false;
	if (!skip_fillboundary) {
		in.FillBoundary(0, ncomp, m_geom[amrlev][mglev].periodicity(),cross); 
	}

	int flagbc = (bc_mode == BCMode::Homogeneous) ? 0 : 1;

	const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

	const auto& maskvals = m_maskvals[amrlev][mglev];
	const auto& bcondloc = *m_bcondloc[amrlev][mglev];

	amrex::FArrayBox foo(amrex::Box::TheUnitBox(),ncomp);

#ifdef _OPENMP
#pragma omp parallel
#endif
	
	//m_bc->SetLevel(amrlev);
	// m_bc->define(m_geom[amrlev][mglev]);
	// m_bc->FillBoundary(in,0,0,0.0);
	for (amrex::MFIter mfi(in, amrex::MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
	{
		//continue;
		const amrex::Box& vbx   = mfi.validbox();
		amrex::FArrayBox& iofab = in[mfi];
		const RealTuple & bdl = bcondloc.bndryLocs(mfi);
		const BCTuple   & bdc = bcondloc.bndryConds(mfi);
		for (amrex::OrientationIter oitr; oitr; ++oitr)
		{
			const amrex::Orientation ori = oitr();
			int  cdr = ori;
			amrex::Real bcl = bdl[ori];
			int  bct = bdc[ori];
			foo.setVal(10.0);
			const amrex::FArrayBox& fsfab = (bndry != nullptr) ? bndry->bndryValues(ori)[mfi] : foo;
			const amrex::Mask& m = maskvals[ori][mfi];
			amrex_mllinop_apply_bc(BL_TO_FORTRAN_BOX(vbx),
					       BL_TO_FORTRAN_ANYD(iofab),
					       BL_TO_FORTRAN_ANYD(m),
					       cdr, bct, bcl,
					       BL_TO_FORTRAN_ANYD(fsfab),
					       maxorder, dxinv, flagbc, ncomp, cross);
		}
	}
}

void
Operator2::reflux (int crse_amrlev,
		   amrex::MultiFab& res, const amrex::MultiFab& crse_sol, const amrex::MultiFab&,
		   amrex::MultiFab&, amrex::MultiFab& fine_sol, const amrex::MultiFab&) const
{
	BL_PROFILE("Operator2::reflux()");
	amrex::YAFluxRegister& fluxreg = m_fluxreg[crse_amrlev];
	fluxreg.reset();

	const int ncomp = getNComp();

	const int fine_amrlev = crse_amrlev+1;

	amrex::Real dt = 1.0;
	const amrex::Real* crse_dx = m_geom[crse_amrlev][0].CellSize();
	const amrex::Real* fine_dx = m_geom[fine_amrlev][0].CellSize();

	const int mglev = 0;
	applyBC(fine_amrlev, mglev, fine_sol, BCMode::Inhomogeneous, m_bndry_sol[fine_amrlev].get());

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		std::array<amrex::FArrayBox,AMREX_SPACEDIM> flux;
		std::array<amrex::FArrayBox*,AMREX_SPACEDIM> pflux { AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) };
		std::array<amrex::FArrayBox const*,AMREX_SPACEDIM> cpflux { AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) };

		for (amrex::MFIter mfi(crse_sol, amrex::MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
		{
			if (fluxreg.CrseHasWork(mfi))
			{
				const amrex::Box& tbx = mfi.tilebox();
				AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
					     flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
					     flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
				FFlux(crse_amrlev, mfi, pflux, crse_sol[mfi]);
				fluxreg.CrseAdd(mfi, cpflux, crse_dx, dt);
			}
		}

#ifdef _OPENMP
#pragma omp barrier
#endif

		for (amrex::MFIter mfi(fine_sol, amrex::MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
		{
			if (fluxreg.FineHasWork(mfi))
			{
				const amrex::Box& tbx = mfi.tilebox();
				AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
					     flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
					     flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
				const int face_only = true;
				FFlux(fine_amrlev, mfi, pflux, fine_sol[mfi], face_only);
				fluxreg.FineAdd(mfi, cpflux, fine_dx, dt);            
			}
		}
	}

	fluxreg.Reflux(res);
}

void
Operator2::compFlux (int amrlev, const std::array<amrex::MultiFab*,AMREX_SPACEDIM>& fluxes,
		     amrex::MultiFab& sol) const
{
	BL_PROFILE("Operator2::compFlux()");

	const int mglev = 0;
	const int ncomp = getNComp();
	applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		std::array<amrex::FArrayBox,AMREX_SPACEDIM> flux;
		std::array<amrex::FArrayBox*,AMREX_SPACEDIM> pflux { AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) };
		for (amrex::MFIter mfi(sol, amrex::MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
		{
			const amrex::Box& tbx = mfi.tilebox();
			AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0));,
				     flux[1].resize(amrex::surroundingNodes(tbx,1));,
				     flux[2].resize(amrex::surroundingNodes(tbx,2)););
			FFlux(amrlev, mfi, pflux, sol[mfi]);
			for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
				const amrex::Box& nbx = mfi.nodaltilebox(idim);
				(*fluxes[idim])[mfi].copy(flux[idim], nbx, 0, nbx, 0, ncomp);
			}
		}
	}
}

void
Operator2::compGrad (int amrlev, const std::array<amrex::MultiFab*,AMREX_SPACEDIM>& grad,
		     amrex::MultiFab& sol) const
{
	BL_PROFILE("Operator2::compGrad()");

	if (sol.nComp() > 1)
		amrex::Abort("Operator2::compGrad called, but only works for single-component solves");

	const int mglev = 0;
	applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

	const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(sol, amrex::MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
	{
		AMREX_D_TERM(const amrex::Box& xbx = mfi.nodaltilebox(0);,
			     const amrex::Box& ybx = mfi.nodaltilebox(1);,
			     const amrex::Box& zbx = mfi.nodaltilebox(2););
		amrex_mllinop_grad(AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
						BL_TO_FORTRAN_BOX(ybx),
						BL_TO_FORTRAN_BOX(zbx)),
				   BL_TO_FORTRAN_ANYD(sol[mfi]),
				   AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*grad[0])[mfi]),
						BL_TO_FORTRAN_ANYD((*grad[1])[mfi]),
						BL_TO_FORTRAN_ANYD((*grad[2])[mfi])),
				   dxinv);
	}
}

void
Operator2::prepareForSolve ()
{
	BL_PROFILE("Operator2::prepareForSolve()");

	const int ncomp = getNComp();
	for (int amrlev = 0;  amrlev < m_num_amr_levels; ++amrlev)
	{
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			const auto& bcondloc = *m_bcondloc[amrlev][mglev];
			const auto& maskvals = m_maskvals[amrlev][mglev];
			const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

			amrex::BndryRegister& undrrelxr = m_undrrelxr[amrlev][mglev];
			amrex::MultiFab foo(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], ncomp, 0, amrex::MFInfo().SetAlloc(false));
#ifdef _OPENMP
#pragma omp parallel
#endif
			for (amrex::MFIter mfi(foo, amrex::MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
			{
				const amrex::Box& vbx = mfi.validbox();

				const RealTuple & bdl = bcondloc.bndryLocs(mfi);
				const BCTuple   & bdc = bcondloc.bndryConds(mfi);

				for (amrex::OrientationIter oitr; oitr; ++oitr)
				{
					const amrex::Orientation ori = oitr();
                    
					int  cdr = ori;
					amrex::Real bcl = bdl[ori];
					int  bct = bdc[ori];
                    
					amrex::FArrayBox& ffab = undrrelxr[ori][mfi];
					const amrex::Mask& m   =  maskvals[ori][mfi];

					amrex_mllinop_comp_interp_coef0(BL_TO_FORTRAN_BOX(vbx),
									BL_TO_FORTRAN_ANYD(ffab),
									BL_TO_FORTRAN_ANYD(m),
									cdr, bct, bcl, maxorder, dxinv, ncomp);
				}
			}
		}
	}
}

amrex::Real
Operator2::xdoty (int amrlev, int mglev,
		  const amrex::MultiFab& x, const amrex::MultiFab& y, bool local) const
{
	const int ncomp = getNComp();
	const int nghost = 0;
	amrex::Real result = amrex::MultiFab::Dot(x,0,y,0,ncomp,nghost,true);
	if (!local) {
		amrex::ParallelAllReduce::Sum(result, Communicator(amrlev, mglev));
	}
	return result;
}

Operator2::BndryCondLoc::BndryCondLoc (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
	: bcond(ba, dm),
	  bcloc(ba, dm)
{
}

void
Operator2::BndryCondLoc::setLOBndryConds (const amrex::Geometry& geom, const amrex::Real* dx,
					  const std::array<BCType,AMREX_SPACEDIM>& lobc,
					  const std::array<BCType,AMREX_SPACEDIM>& hibc,
					  int ratio, const amrex::RealVect& a_loc)
{
	const amrex::Box&  domain = geom.Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(bcloc); mfi.isValid(); ++mfi)
	{
		const amrex::Box& bx = mfi.validbox();
		RealTuple & bloc  = bcloc[mfi];
		BCTuple   & bctag = bcond[mfi];

		amrex::MLMGBndry::setBoxBC(bloc, bctag, bx, domain, lobc, hibc, dx, ratio, a_loc);
	}
}

void
Operator2::applyMetricTerm (int, int, amrex::MultiFab&) const
{
	// This is a required method needed only if the geometry is 
	// non-Cartesian. This operator is used for Cartesian coordinates
	// only.
}

void
Operator2::unapplyMetricTerm (int, int, amrex::MultiFab&) const
{
	// This is a required method needed only if the geometry is 
	// non-Cartesian. This operator is used for Cartesian coordinates
	// only.
}

}