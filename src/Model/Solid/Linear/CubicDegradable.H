#ifndef MODEL_SOLID_LINEAR_CUBICDEGRADABLE_H_
#define MODEL_SOLID_LINEAR_CUBICDEGRADABLE_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Linear/Cubic.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class CubicDegradable : public Cubic
{
public:
    CubicDegradable() {};
	CubicDegradable(Solid<Set::Sym::MajorMinor> base) : Cubic(base) {};    
	CubicDegradable(Cubic base) : Cubic(base) {};
    //CubicDegradable(Cubic base) : Cubic(base) {};
    CubicDegradable(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
    {
    	Define(C11,C11,C12,C12,C44,C44,phi1,phi1,Phi,Phi,phi2,phi2);
    };


	void Define(Set::Scalar C11, Set::Scalar C11_0, Set::Scalar C12, Set::Scalar C12_0, Set::Scalar C44, Set::Scalar C44_0, Set::Scalar phi1,Set::Scalar phi1_0, Set::Scalar Phi, Set::Scalar Phi_0, Set::Scalar phi2, Set::Scalar phi2_0)
	{
		m_C11_0 = C11_0; m_C12_0 = C12_0; m_C44_0 = C44_0;
		m_phi1 = phi1_0; m_Phi = Phi_0; m_phi2 = phi2_0;
        Cubic::Define(C11,C12,C44, phi1, Phi,phi2);
	}

    virtual ~CubicDegradable() {};

    void DegradeModulus(const Set::Scalar eta)
	{
		if(std::isnan(eta)) Util::Abort(INFO, "Nans in eta");
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta ", eta);
		if(std::isnan(m_C11_0)) Util::Abort(INFO, "Nans in lambda0");
		if(std::isnan(m_C12_0)) Util::Abort(INFO, "Nans in mu0");
		if(std::isnan(m_C44_0)) Util::Abort(INFO, "Nans in mu0");

		Cubic::Define( m_C11_0*(1.-eta), m_C12_0*(1.-eta),m_C44_0*(1.-eta), m_phi1,m_Phi,m_phi2);
	}

    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
	{
		if(eta.size() == 0) return;
		else if(eta.size() == 1)
		{
			if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
			Cubic::Define( m_C11_0*(1.-eta[0]), m_C12_0*(1.-eta[0]),m_C44_0*(1.-eta[0]), m_phi1,m_Phi,m_phi2);
		}
		else if(eta.size() == 2)
		{
			if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
			Cubic::Define( m_C11_0*(1.-eta[0]), m_C12_0*(1.-eta[1]),m_C44_0*(1.-eta[0]), m_phi1,m_Phi,m_phi2);
		}
		else
			return;
	}

private:
    Set::Scalar m_C11_0 = NAN, m_C12_0 = NAN, m_C44_0 = NAN, m_phi1 = NAN, m_Phi = NAN, m_phi2 = NAN;

public:
    static void Parse(CubicDegradable & value, IO::ParmParse & pp)
    {
        Set::Scalar C11, C12, C44;
        pp.query("C11",C11);
        pp.query("C12",C12);
        pp.query("C44",C44);

        if (pp.contains("random"))
        {
            value = Cubic::Random(C11,C12,C44);
            return;
        }

        Set::Scalar phi1, Phi, phi2;
        pp.query("phi1",phi1);
        pp.query("Phi",Phi);
        pp.query("phi2",phi2);
        value.Define(C11,C11,C12,C12,C44,C44,phi1,phi1,Phi,Phi,phi2,phi2);
    }
};

}
}
}

#endif