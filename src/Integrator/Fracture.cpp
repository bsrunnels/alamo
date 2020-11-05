#include "Fracture.H"

namespace Integrator
{
Fracture::Fracture() :
	Integrator()
{
    //==================================================
    // Problem type
    {
        std::string fracture_problem;
        IO::ParmParse pp("fracture");
        pp.query("problem_type", fracture_problem);

        if (fracture_problem == "brittle") fracture_type = FractureType::Brittle;
        else Util::Abort(INFO, "We are only working with brittle fracture right now");
        
    }
    //==================================================

    //==================================================
    // Crack model
    {
        IO::ParmParse pp_crack("crack");
        std::string crack_type;
        pp_crack.query("type",crack_type);
        pp_crack.query("modulus_scaling_max",crack.scaleModulusMax);
        pp_crack.query("refinement_threshold",crack.refinement_threshold);
        pp_crack.query("df_tol_rel",crack.driving_force_tolerance_rel);
        pp_crack.query("df_tol_abs",crack.driving_force_tolerance_abs);
        
        if(crack_type=="constant")
        {
            Model::Interface::Crack::Constant *tmpbdy = new Model::Interface::Crack::Constant();
            pp_crack.queryclass("constant",*tmpbdy);
            crack.cracktype = tmpbdy;
        }
        else if(crack_type == "sin")
        {
            Model::Interface::Crack::Sin *tmpbdy = new Model::Interface::Crack::Sin();
            pp_crack.queryclass("sin", *tmpbdy);
            crack.cracktype = tmpbdy;
        }
        else
            Util::Abort(INFO,"This crack model hasn't been implemented yet");
    }
    //==================================================

    //==================================================
    // Anistropy
    {
        IO::ParmParse pp_anisotropy("crack.anisotropy");
        pp_anisotropy.query("on", anisotropy.on);
        pp_anisotropy.query("tstart",anisotropy.tstart);
        anisotropy.timestep = timestep;
        pp_anisotropy.query("timestep", anisotropy.timestep);
        anisotropy.plot_int = plot_int;
        pp_anisotropy.query("plot_int", anisotropy.plot_int);
        anisotropy.plot_dt = plot_dt;
        pp_anisotropy.query("plot_dt", anisotropy.plot_dt);
        pp_anisotropy.query("beta", anisotropy.beta);
    }
    //==================================================

    //==================================================
    // ICs
    {
        IO::ParmParse pp("ic"); // Phase-field model parameters
        pp.query("type", crack.ic_type);

        if(crack.ic_type == "ellipsoid")
        {
            IC::Ellipsoid *tmpic = new IC::Ellipsoid(geom);
            pp.queryclass("ellipsoid", *tmpic);
            crack.ic = tmpic;
            crack.is_ic = true;
        }
        else if(crack.ic_type == "notch")
        {
            IC::Notch *tmpic = new IC::Notch(geom);
            pp.queryclass("notch",*tmpic);
            crack.ic = tmpic;
            crack.is_ic = true;
        }
        else 
        {
            Util::Warning(INFO, "No valid IC found. Ignoring ICs");
            crack.is_ic = false;
        }
    }
    //==================================================

    //==================================================
    RegisterNodalFab(crack.field, 1, number_of_ghost_nodes, "crack", true);
    RegisterNodalFab(crack.field_old, 1, number_of_ghost_nodes, "crack_old", true);
    RegisterNodalFab(crack.driving_force, 5, number_of_ghost_nodes, "driving_force", true);
    RegisterIntegratedVariable(&(crack.driving_force_norm),"driving_force_norm");
    
    //==================================================

    //==================================================
    // Material model
    {
        IO::ParmParse pp_material("material");
        pp_material.query("model",material.input_material);

        switch (fracture_type)
        {
            case FractureType::Brittle:
                if(material.input_material == "isotropic") pp_material.queryclass("isotropic",material.brittlemodeltype);
                else Util::Abort(INFO,"This model has not been implemented yet.");
                break;
            default:
                break;
        }

        IO::ParmParse pp_material_ic("material.ic");
        std::string ic_type;
        pp_material_ic.query("type",ic_type);
        if (ic_type == "ellipsoid")
        {
            IC::Ellipsoid *tmpic = new IC::Ellipsoid(geom);
            pp_material_ic.queryclass("ellipsoid", *tmpic);
            material.ic = tmpic;
            material.is_ic = true;
        }
        else if(ic_type == "notch")
        {
            IC::Notch *tmpic = new IC::Notch(geom);
            pp_material_ic.queryclass("notch",*tmpic);
            material.ic = tmpic;
            material.is_ic = true;
        }
        else 
        {
            Util::Warning(INFO, "No valid IC found. Ignoring ICs");
            material.is_ic = false;
        }
    }
    //==================================================

    //==================================================
    // Solver properties
    {
        IO::ParmParse pp_elastic("solver");
        pp_elastic.query("int",				sol.interval);
        pp_elastic.query("type",			sol.type);
        pp_elastic.query("max_iter",		sol.max_iter);
        pp_elastic.query("max_fmg_iter",	sol.max_fmg_iter);
        pp_elastic.query("verbose",			sol.verbose);
        pp_elastic.query("cgverbose",		sol.cgverbose);
        pp_elastic.query("tol_rel",			sol.tol_rel);
        pp_elastic.query("tol_abs",			sol.tol_abs);
        pp_elastic.query("cg_tol_rel",		sol.cg_tol_rel);
        pp_elastic.query("cg_tol_abs",		sol.cg_tol_abs);
        pp_elastic.query("use_fsmooth",		sol.use_fsmooth);
        pp_elastic.query("agglomeration", 	sol.agglomeration);
        pp_elastic.query("consolidation", 	sol.consolidation);

        pp_elastic.query("bottom_solver",       sol.bottom_solver);
        pp_elastic.query("linop_maxorder",      sol.linop_maxorder);
        pp_elastic.query("max_coarsening_level",sol.max_coarsening_level);
        pp_elastic.query("verbose",             sol.verbose);
        pp_elastic.query("cg_verbose",          sol.cgverbose);
        pp_elastic.query("bottom_max_iter",     sol.bottom_max_iter);
        pp_elastic.query("max_fixed_iter",      sol.max_fixed_iter);
        pp_elastic.query("bottom_tol",          sol.bottom_tol);
        pp_elastic.query("pre_smooth", sol.pre_smooth);
        pp_elastic.query("post_smooth", sol.post_smooth);
    }
    //==================================================

    //==================================================
    // Loading conditions and multipliers
    {
        IO::ParmParse pp_elastic("elastic");
        pp_elastic.query("df_mult", elastic.df_mult);

        std::string loading_mode, loading_type;
        IO::ParmParse pp_load("loading");
        if (pp_load.countval("body_force")) pp_load.queryarr("body_force",loading.body_force);

        if( fracture_type == FractureType::Brittle ) pp_elastic.queryclass("bc",elastic.brittlebc);
    }
    //==================================================

    //==================================================
    // Registering fabs now
    {
        nlevels = maxLevel() + 1;
        RegisterNodalFab(elastic.disp,  AMREX_SPACEDIM, number_of_ghost_nodes, "disp", true);
        RegisterNodalFab(elastic.rhs,  AMREX_SPACEDIM, number_of_ghost_nodes, "rhs", true);
        RegisterNodalFab(elastic.residual,  AMREX_SPACEDIM, number_of_ghost_nodes, "res", true);
        RegisterNodalFab(elastic.strain,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strain", true);
        RegisterNodalFab(elastic.stress,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "stress", true);
        RegisterNodalFab(elastic.energy, 1, number_of_ghost_nodes, "energy", true);
        RegisterNodalFab(elastic.energy_pristine, 1, number_of_ghost_nodes, "energy_pristine", true);
        RegisterNodalFab(elastic.energy_pristine_old, 1, number_of_ghost_nodes, "energy_pristine_old", true);

        if (fracture_type == FractureType::Brittle)
        {
            RegisterGeneralFab(material.brittlemodel, 1, number_of_ghost_nodes);
            RegisterNodalFab(material.modulus_field, 1, number_of_ghost_nodes, "modulus_field", true);
            material.brittlemodel.resize(nlevels);
        }
    }
    //==================================================
}

Fracture::~Fracture()
{
}

void
Fracture::Initialize (int ilev)
{
    //==================================================
    // Initialization of crack fields
    {
        if (crack.is_ic)
        {
            crack.ic->Initialize(ilev,crack.field);
            crack.ic->Initialize(ilev,crack.field_old);
        }
        else
        {
            crack.field[ilev]->setVal(1.0);
            crack.field_old[ilev]->setVal(1.0);
        }
        crack.driving_force[ilev]->setVal(0.0);
    }
    //==================================================
	
    //==================================================
    // Initialization of elastic fields
    {
        elastic.disp[ilev]->setVal(0.0);
        elastic.strain[ilev]->setVal(0.0);
        elastic.stress[ilev]->setVal(0.0);
        elastic.rhs[ilev]->setVal(0.0);
        elastic.energy[ilev]->setVal(0.0);
        elastic.residual[ilev]->setVal(0.0);
        elastic.energy_pristine[ilev] -> setVal(0.);
        elastic.energy_pristine_old[ilev] -> setVal(0.);
    }
	//==================================================

    //==================================================
    // Initialization of brittle and ductile specific fields
    {
        if(material.is_ic) material.ic->Initialize(ilev,material.modulus_field);
        else material.modulus_field[ilev]->setVal(1.0);

        if (fracture_type == FractureType::Brittle)
            material.brittlemodel[ilev]->setVal(material.brittlemodeltype);
    }
    //==================================================
}

void
Fracture::TimeStepBegin(amrex::Real time, int iter)
{
    Util::Message(INFO,crack.driving_force_norm," ",crack.driving_force_reference," ",crack.driving_force_tolerance_rel);
    if (crack.driving_force_norm / crack.driving_force_reference < crack.driving_force_tolerance_rel)
        elastic.do_solve_now = true;
    if (crack.driving_force_norm < crack.driving_force_tolerance_abs)
        elastic.do_solve_now = true;

    if (!elastic.do_solve_now) return;

    if (anisotropy.on && time >= anisotropy.tstart)
	{
		SetTimestep(anisotropy.timestep);
		if (anisotropy.elastic_int > 0) 
			if (iter % anisotropy.elastic_int) return;
	}
    //if(iter%sol.interval) return;
    loading.val = loading.init + ((double)loading.step)*loading.rate;
    material.brittlemodeltype.DegradeModulus(0.0);

    for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		std::swap(elastic.energy_pristine_old[ilev], elastic.energy_pristine[ilev]);
        Util::RealFillBoundary(*crack.field[ilev],geom[ilev]);
        if (fracture_type == FractureType::Brittle) 
        {
            material.brittlemodel[ilev]->setVal(material.brittlemodeltype);
            Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
        }
    }

    //==================================================
    // Scaling the modulus
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            Util::RealFillBoundary(*crack.field[ilev],geom[ilev]);
            Util::RealFillBoundary(*crack.field_old[ilev],geom[ilev]);
            Util::RealFillBoundary(*material.modulus_field[ilev],geom[ilev]);

            for (amrex::MFIter mfi(*elastic.disp[ilev],true); mfi.isValid(); ++mfi)
            {
                amrex::Box box = mfi.grownnodaltilebox();
                amrex::Array4<const Set::Scalar> const& c_new = (*crack.field[ilev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[ilev]).array(mfi);

                amrex::Array4<brittle_fracture_model_type> modelfab_b;
                
                if (fracture_type == FractureType::Brittle)
                    modelfab_b = (material.brittlemodel)[ilev]->array(mfi);
                
                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Scalar _temp = 0;
                    
                    if (fracture_type == FractureType::Brittle)
                        _temp = crack.cracktype->g_phi(c_new(i,j,k,0),0);
                    
                    _temp = std::min(_temp, crack.cracktype->g_phi(modbox(i,j,k,0),0.0));
                    
                    if (std::isnan(_temp)) Util::Abort(INFO);
                    if(_temp < 0.0) _temp = 0.;
                    if(_temp > 1.0) _temp = 1.0;
                    _temp = crack.scaleModulusMax + _temp * (1. - crack.scaleModulusMax);

                    if (fracture_type == FractureType::Brittle)
                        modelfab_b(i,j,k,0).DegradeModulus(1.-_temp);
                });
            }
            if (fracture_type == FractureType::Brittle) Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
        }
    }
    //==================================================

    //==================================================
    // Setting the body forces
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            const Real* DX = geom[ilev].CellSize();
            Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

            AMREX_D_TERM(elastic.rhs[ilev]->setVal(loading.body_force(0)*volume,0,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(1)*volume,1,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(2)*volume,2,1););
        }
    }
    //==================================================

    //==================================================
    // Setting the elastic boundary conditions
    {
        if (fracture_type == FractureType::Brittle)
            elastic.brittlebc.Init(elastic.rhs,geom);
    }
    //==================================================

    //==================================================
    // Setting up the solver parameters
    Operator::Elastic<brittle_fracture_model_type::sym> op_b;
    {
        LPInfo info;
        info.setAgglomeration(sol.agglomeration);
        info.setConsolidation(sol.consolidation);
        info.setMaxCoarseningLevel(sol.max_coarsening_level);

        for (int ilev = 0; ilev < nlevels; ilev++) if (elastic.disp[ilev]->contains_nan()) Util::Warning(INFO);

        if (fracture_type == FractureType::Brittle)
        {
            op_b.define(geom, grids, dmap, info);
            op_b.setMaxOrder(sol.linop_maxorder);
            op_b.SetBC(&elastic.brittlebc);

            Solver::Nonlocal::Newton<brittle_fracture_model_type>  solver(op_b);
            solver.setMaxIter(sol.max_iter);
            solver.setMaxFmgIter(sol.max_fmg_iter);
            solver.setFixedIter(sol.max_fixed_iter);
            solver.setVerbose(sol.verbose);
            solver.setBottomVerbose(sol.cgverbose);
            solver.setBottomMaxIter(sol.bottom_max_iter);
            solver.setBottomTolerance(sol.cg_tol_rel) ;
            solver.setBottomToleranceAbs(sol.cg_tol_abs) ;
            if (sol.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
            else if (sol.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
            solver.solve(elastic.disp, elastic.rhs, material.brittlemodel, sol.tol_rel, sol.tol_abs);
            solver.compResidual(elastic.residual,elastic.disp,elastic.rhs,material.brittlemodel);
        }
    }
    //==================================================
    
    //==================================================
    // Computing new stresses, strains and energies
    {
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            if (fracture_type == FractureType::Brittle)
            {
                op_b.Strain(ilev,*elastic.strain[ilev],*elastic.disp[ilev]);
                op_b.Stress(ilev,*elastic.stress[ilev],*elastic.disp[ilev]);
                op_b.Energy(ilev,*elastic.energy[ilev],*elastic.disp[ilev]);
            }
        }
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            Util::RealFillBoundary(*elastic.strain[ilev],geom[ilev]);
            Util::RealFillBoundary(*elastic.energy_pristine_old[ilev],geom[ilev]);
            Util::RealFillBoundary(*material.modulus_field[ilev],geom[ilev]);
            elastic.energy_pristine[ilev]->setVal(0.0);

            for (amrex::MFIter mfi(*elastic.strain[ilev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.grownnodaltilebox();
                amrex::Array4<const Set::Scalar> const& strain_box 	= (*elastic.strain[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_box 		= (*elastic.energy_pristine[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_box_old 	= (*elastic.energy_pristine_old[ilev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[ilev]).array(mfi);

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Matrix eps = Numeric::FieldToMatrix(strain_box,i,j,k);
                    Eigen::SelfAdjointEigenSolver<Set::Matrix> eigensolver(eps);
                    Set::Vector eValues = eigensolver.eigenvalues();
                    Set::Matrix eVectors = eigensolver.eigenvectors();

                    Set::Matrix epsp = Set::Matrix::Zero();
                    Set::Matrix epsn = Set::Matrix::Zero();

                    for (int n = 0; n < AMREX_SPACEDIM; n++)
                    {
                        if(eValues(n) > 0.0) epsp += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                        else epsn += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                    }
                    material.brittlemodeltype.DegradeModulus(crack.scaleModulusMax + crack.cracktype->g_phi(modbox(i,j,k,0),0.0) * (1. - crack.scaleModulusMax));
                    energy_box(i,j,k,0) = material.brittlemodeltype.W(epsp);
                    energy_box(i,j,k,0) = energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);
                    
                    if (std::isnan(energy_box(i,j,k,0))) Util::Abort(INFO, "Nans detected in energy_box. material.brittlemodeltype.W(eps) = ", material.brittlemodeltype.W(eps));
                });
            }
            Util::RealFillBoundary(*elastic.energy_pristine[ilev],geom[ilev]);
        }
    }
    //==================================================

    integrate_variables_before_advance = false;
    integrate_variables_after_advance = true;
}

void 
Fracture::Advance (int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
    std::swap(*crack.field_old[lev], *crack.field[lev]);crack.field_old[lev]->FillBoundary();
    Util::RealFillBoundary(*crack.field[lev],geom[lev]);

	const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain(geom[lev].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());
    const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

    for ( amrex::MFIter mfi(*crack.field[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.validbox();
		amrex::Array4<const Set::Scalar> const& c_old = (*crack.field_old[lev]).array(mfi);
		amrex::Array4<Set::Scalar> const& df = (*crack.driving_force[lev]).array(mfi);
		amrex::Array4<Set::Scalar> const& c_new = (*crack.field[lev]).array(mfi);
		amrex::Array4<const Set::Scalar> const& energy_box = (*elastic.energy_pristine[lev]).array(mfi);
        // amrex::Array4<const Set::Scalar> const& modbox = (*material.modulus_field[lev]).array(mfi);

        amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
#if AMREX_SPACEDIM !=2
            Util::Abort(INFO, "This doesn't work with 1D or 3D yet");
#endif
            if (i == lo.x && j == lo.y) c_new(i,j,k,0) = c_new(i+1,j+1,k,0);
            else if (i == lo.x && j == hi.y) c_new(i,j,k,0) = c_new(i+1,j-1,k,0);
            else if (i == hi.x && j == lo.y) c_new(i,j,k,0) = c_new(i-1,j+1,k,0);
            else if (i == hi.x && j == hi.y) c_new(i,j,k,0) = c_new(i-1,j-1,k,0);
            else if (i == lo.x) c_new(i,j,k,0) = c_new(i+1,j,k,0);
            else if (j == lo.y) c_new(i,j,k,0) = c_new(i,j+1,k,0);
            else if (i == hi.x) c_new(i,j,k,0) = c_new(i-1,j,k,0);
            else if (j == hi.y) c_new(i,j,k,0) = c_new(i,j-1,k,0);
            else
            {
			    Set::Scalar rhs = 0.0;

                // Elastic component of the driving force
                Set::Scalar en_cell = energy_box(i,j,k,0);

                if (std::isnan(en_cell)) Util::Abort(INFO, "Nans detected in en_cell. energy_box(i,j,k,0) = ", energy_box(i,j,k,0));
                if (std::isinf(c_old(i,j,k,0))) Util::Abort(INFO, "Infs detected in c_old");
                if (c_old(i,j,k,0) > 1) Util::Abort(INFO, "Very large values of c_old at (",i,",",j,",",k,") = ", c_old(i,j,k,0));

                df(i,j,k,0) = crack.cracktype->Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;
                rhs += crack.cracktype->Dg_phi(c_old(i,j,k,0),0.0)*en_cell*elastic.df_mult;

                Set::Matrix DDc = Numeric::Hessian(c_old, i, j, k, 0, DX);
                Set::Scalar laplacian = DDc.trace();

                // Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDEta = Numeric::DoubleHessian<AMREX_SPACEDIM>(c_old, i, j, k, 0, DX);
                Set::Scalar bilaplacian = 0.0;
                // for (int p = 0; p < AMREX_SPACEDIM; p++)
                //     for (int q =0; q < AMREX_SPACEDIM; q++)
                //         bilaplacian += DDDDEta(p,p,q,q);
                
                if (anisotropy.on) Util::Abort(INFO, "Anisotropy not implemented yet");

                df(i,j,k,1) = crack.cracktype->Gc(0.0)*crack.cracktype->Dw_phi(c_old(i,j,k,0),0.0)/(4.0*crack.cracktype->Zeta(0.0))*crack.mult_df_Gc;
                df(i,j,k,2) = 2.0*crack.cracktype->Zeta(0.0)*crack.cracktype->Gc(0.0)*laplacian*crack.mult_df_lap;
                df(i,j,k,3) = 0.5*(crack.cracktype->Zeta(0.0)*crack.cracktype->Zeta(0.0)*crack.cracktype->Zeta(0.0))*bilaplacian;

                rhs += crack.cracktype->Gc(0.0)*crack.cracktype->Dw_phi(c_old(i,j,k,0),0.)/(4.0*crack.cracktype->Zeta(0.0))*crack.mult_df_Gc;
                rhs -= 2.0*crack.cracktype->Zeta(0.0)*crack.cracktype->Gc(0.0)*laplacian*crack.mult_df_lap;
                // rhs += 0.5*(crack.cracktype->Zeta(0.0)*crack.cracktype->Zeta(0.0)*crack.cracktype->Zeta(0.0))*bilaplacian;

                // rhs *= crack.cracktype->g_phi(modbox(i,j,k,0),0.0);                
                // rhs *= modbox(i,j,k,0); 

                if (fracture_type == FractureType::Brittle)
                    df(i,j,k,4) = std::max(0.,rhs - crack.cracktype->DrivingForceThreshold(c_old(i,j,k,0)));
            
			    if(std::isnan(rhs)) Util::Abort(INFO, "Dwphi = ", crack.cracktype->Dw_phi(c_old(i,j,k,0),0.0),". c_old(i,j,k,0) = ",c_old(i,j,k,0));
			    c_new(i,j,k,0) = c_old(i,j,k,0) - dt*std::max(0., rhs - crack.cracktype->DrivingForceThreshold(c_old(i,j,k,0)))*crack.cracktype->Mobility(c_old(i,j,k,0));

                if(c_new(i,j,k,0) > 1.0) {/*Util::Message(INFO, "cnew = ", c_new(i,j,k,0) ,", resetting to 1.0");*/ c_new(i,j,k,0) = 1.;}
                if(c_new(i,j,k,0) < 0.0) {/*Util::Message(INFO, "cnew = ", c_new(i,j,k,0) ,", resetting to 0.0");*/ c_new(i,j,k,0) = 0.;}
            }
		});
    }
    Util::RealFillBoundary(*crack.field[lev],geom[lev]);
}

void 
Fracture::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,const amrex::MFIter &mfi, const amrex::Box &box)
{
	const amrex::Real* DX = geom[amrlev].CellSize();
    const Set::Scalar DV = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

	amrex::Array4<const Set::Scalar> const &df = (*crack.driving_force[amrlev]).array(mfi);
    
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
    {
        crack.driving_force_norm += df(i,j,k,4) * DV;
	});
}

void
Fracture::TagCellsForRefinement (int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real *DX = geom[lev].CellSize();
	const Set::Vector dx(DX);
	const Set::Scalar dxnorm = dx.lpNorm<2>();
    amrex::Box domain(geom[lev].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());


	for (amrex::MFIter mfi(*crack.field[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box 							&bx 	= mfi.tilebox();
		amrex::Array4<char> const 					&tags 	= a_tags.array(mfi);
		amrex::Array4<const Set::Scalar> const 		&c_new 	= (*crack.field[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const 		&mod_box 	= (*material.modulus_field[lev]).array(mfi);
		
		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            std::array<Numeric::StencilType,AMREX_SPACEDIM> sten 
                = Numeric::GetStencil(i,j,k,domain);

			Set::Vector grad = Numeric::Gradient(c_new, i, j, k, 0, DX, sten);
            Set::Vector grad_mod = Numeric::Gradient(mod_box, i, j, k, 0, DX, sten);
			if (dxnorm * grad.lpNorm<2>() > crack.refinement_threshold || dxnorm * grad_mod.lpNorm<2>() > crack.refinement_threshold)
				tags(i, j, k) = amrex::TagBox::SET;
		});
	}
}

void
Fracture::TimeStepComplete(amrex::Real /*time*/,int /*iter*/)
{
    if (elastic.do_solve_now)
        crack.driving_force_reference = crack.driving_force_norm;

    elastic.do_solve_now = false;
}

}
