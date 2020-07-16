#ifndef BC_OPERATOR_ELASTIC_CRACK_H
#define BC_OPERATOR_ELASTIC_CRACK_H

// #include "Operator/Elastic.H"
#include "IO/ParmParse.H"
#include "Numeric/Interpolator/Linear.H"
#include "BC/Operator/Elastic/Elastic.H"

namespace BC
{
namespace Operator
{
namespace Elastic
{
class Crack : public Elastic
{
public:

    Crack() 
    {

    }
 	Crack(Set::Scalar disp)
	{
        m_disp = disp;
	};
 	~Crack() {};

	using Elastic::Init;
	virtual void
	Init(amrex::MultiFab * a_rhs,
	     const amrex::Geometry &a_geom,
		 bool a_homogeneous = false) const
	{
			amrex::Box domain(a_geom.Domain());
			domain.convert(amrex::IntVect::TheNodeVector());
			const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
			for (amrex::MFIter mfi(*a_rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
			{
				amrex::Box bx = mfi.tilebox();
				bx.grow(2);
				bx = bx & domain;
				amrex::Array4<amrex::Real> const& rhs       = a_rhs->array(mfi);
				amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    rhs(i,j,k,0) = 0.0;
                    rhs(i,j,k,1) = 0.0;

					if (j==hi.y && !a_homogeneous)  rhs(i,j,k,1) = m_disp(m_time);

			});
		}					
	}

    //
    // ┌┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┐
    // │                    ├
    // │                    ├
    // │                    ├
    // │                    ├
    // └─────────╮          ├
    //           ╰┬┬┬┬┬┬┬┬┬┬┘
    //

 	AMREX_FORCE_INLINE
 	Set::Vector operator () (const Set::Vector &u,
 				 const Set::Matrix &gradu,
 				 const Set::Matrix &sigma,
 				 const int &i, const int &j, const int &k,
				 const amrex::Box &domain)
 	{
		(void)i; (void)j; (void)k; // Suppress "unused variable" warnings
		//Set::Vector f;

		const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
		
		amrex::IntVect m(AMREX_D_DECL(i,j,k));

		// Corners
		#if AMREX_SPACEDIM == 2
		
		if (m[0] == lo.x && m[1] == lo.y)                  return set({Type::Traction,     Type::Traction},         u, gradu, sigma, Set::Vector(-SQRT2INV, -SQRT2INV));
		if (m[0] == lo.x && m[1] == hi.y)                  return set({Type::Traction, Type::Displacement},     u, gradu, sigma, Set::Vector(-SQRT2INV, +SQRT2INV));
		if (m[0] == hi.x && m[1] == lo.y)                  return set({Type::Displacement, Type::Displacement},     u, gradu, sigma, Set::Vector(+SQRT2INV, -SQRT2INV));
		if (m[0] == hi.x && m[1] == hi.y)                  return set({Type::Displacement, Type::Displacement},     u, gradu, sigma, Set::Vector(+SQRT2INV, +SQRT2INV));
		
		if (m[0] == lo.x)                                  return set({Type::Traction,         Type::Traction},         u, gradu, sigma, Set::Vector(-1, 0));
		if (m[0] == hi.x)                                  return set({Type::Displacement,     Type::Traction},         u, gradu, sigma, Set::Vector(+1, 0));
		if (m[1] == hi.y)                                  return set({Type::Traction,     Type::Displacement},     u, gradu, sigma, Set::Vector( 0,+1));
		if (m[1] == lo.y)
        {
            if (m[0] < domain.hiVect()[0]/4)
            {
                //Util::Message(INFO,"Crack",m);
                return set({Type::Traction, Type::Traction},             u, gradu, sigma, Set::Vector( 0,-1));
            }
            else
            {
                //Util::Message(INFO,"No Crack",m);
                return set({Type::Traction, Type::Displacement},         u, gradu, sigma, Set::Vector( 0,-1));
            }
        }
		
		#elif AMREX_SPACEDIM == 3
        Util::Abort(INFO,"BC::Operator::Elastic::Crack currently works in 2D only");
		#endif
		
		Util::Abort(INFO,"Boundary condition error");
		return Set::Vector::Zero();
 	}

	AMREX_FORCE_INLINE
	Set::Vector set(const std::array<Type,AMREX_SPACEDIM> &bc_type, 
					const Set::Vector &u, const Set::Matrix &gradu, const Set::Matrix &sigma, Set::Vector n) const
	{
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
		{
			if      (bc_type[i] == Type::Displacement) 
				f(i) = u(i);
			else if (bc_type[i] == Type::Traction)
				f(i) = (sigma*n)(i);
			else if (bc_type[i] == Type::Neumann)
				f(i) = (gradu*n)(i);
			else if (bc_type[i] == Periodic)
				continue;
		}
		return f;
	}

protected:
	
	#if AMREX_SPACEDIM==2
	static const int m_nfaces = 8;
	#elif AMREX_SPACEDIM==3
	static const int m_nfaces = 26;
	#endif

	//std::array<std::array<Type,                                      AMREX_SPACEDIM>, m_nfaces> m_bc_type; 
	//std::array<std::array<Numeric::Interpolator::Linear<Set::Scalar>,AMREX_SPACEDIM>, m_nfaces> m_bc_val; 

public:
    static void Parse(Crack & value, IO::ParmParse & pp)
    {
        Util::Abort(INFO,"Not implemented yet");
    }

    Numeric::Interpolator::Linear<Set::Scalar> m_disp;
};
}
}
}
#endif
