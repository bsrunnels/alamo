
#include "PFFem.H"
//#include <AMReX_LP_F.H>

namespace amrex {

PFFem::PFFem (const BndryData& bd,
	      Real             _h)
  :
  LinOp(bd,_h) {}

PFFem::~PFFem() {}

Real
PFFem::norm (int nm, int level, const bool local)
{
  // switch ( nm )
  //   {
  //   case 0:
  //     return 8.0/(h[level][0]*h[level][0]);
  //   }
  // amrex::Error("Bad PFFem::norm");
  return -1.0;
}

void
PFFem::compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
		 MultiFab& in, const BC_Mode& bc_mode,
		 int src_comp, int dst_comp, int num_comp, int bnd_comp)
{
//  BL_PROFILE("PFFem::compFlux()");
//
//  const int level    = 0;
//  applyBC(in,src_comp,num_comp,level,bc_mode,bnd_comp);
//
//  const bool tiling = true;
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
//  for (MFIter inmfi(in,tiling); inmfi.isValid(); ++inmfi)
//    {
//      const Box&       tbx  = inmfi.tilebox();
//      FArrayBox&       xfluxfab = xflux[inmfi];
//      FArrayBox&       yfluxfab = yflux[inmfi];
//      const FArrayBox& ufab = in[inmfi];
//
//      for (int p1 = tbx.loVect()[0]; p1<=tbx.hiVect()[0]; p1++)
//	for (int p2 = tbx.loVect()[1]; p2<=tbx.hiVect()[1]; p2++)
//	  for (int n1 = p1-1; n1 <= p1+1; n1++)
//	    for (int n2 = p2-1; n2 <= p2+1; n2++)
//	      {
//		//for (int i = 0; i < BL_SPACEDIM; i++)
//		//for (int k = 0; k < BL_SPACEDIM; k++)
//
//
//		Real K = (mu + lambda)*Phi(p1,p2,n1,n2,i,k); // This is K(p,n,i,k)
//		if (i==k) K += mu*(Phi(p1,p2,n1,n2,0,0)+Phi(p1,p2,n1,n2,1,1));
//
//		      
//
//		bfab(amrex::IntVect(p1,p2),i) = K*ufab(amrex::IntVect(n1,n2),k);
//	      }
//
//
//
//
//      for (int p1 = tbx.loVect()[0]; p1<=tbx.hiVect()[0]; p1++)
//	for (int p2 = tbx.loVect()[1]; p2<=tbx.hiVect()[1]; p2++)
//	  for (int i = 0; i < BL_SPACEDIM; i++)
//	    {
//	      bfab(amrex::IntVect(p1,p2),i) = 0;
//	      for (int n1 = p1-1; n1 <= p1+1; n1++)
//		for (int n2 = p2-1; n2 <= p2+1; n2++)
//		  {
//		    // int k = 0;
//		    // {
//		    //   Real K = (mu + lambda)*Phi(p1,p2,n1,n2,i,k);
//		    //   if (i==k) K += mu*(Phi(p1,p2,n1,n2,0,0)+Phi(p1,p2,n1,n2,1,1));
//		    //   xfluxfab(amrex::IntVect(p1,p2)) = K*ufab(amrex::IntVect(n1,n2),k);
//		    // }
//
//		  }
//	    }
//    }
//
//
}

//void
//PFFem::Fsmooth (MultiFab&       solnL,
//                    const MultiFab& rhsL,
//                    int             level,
//                    int             redBlackFlag)
//{
//    BL_PROFILE("PFFem::Fsmooth()");
//
//    OrientationIter oitr;
//
//    const FabSet& f0 = undrrelxr[level][oitr()]; oitr++;
//    const FabSet& f1 = undrrelxr[level][oitr()]; oitr++;
//    const FabSet& f2 = undrrelxr[level][oitr()]; oitr++;
//    const FabSet& f3 = undrrelxr[level][oitr()]; oitr++;
//#if (BL_SPACEDIM > 2)
//    const FabSet& f4 = undrrelxr[level][oitr()]; oitr++;
//    const FabSet& f5 = undrrelxr[level][oitr()]; oitr++;
//#endif
//
//    oitr.rewind();
//    const MultiMask& mm0 = maskvals[level][oitr()]; oitr++;
//    const MultiMask& mm1 = maskvals[level][oitr()]; oitr++;
//    const MultiMask& mm2 = maskvals[level][oitr()]; oitr++;
//    const MultiMask& mm3 = maskvals[level][oitr()]; oitr++;
//#if (BL_SPACEDIM > 2)
//    const MultiMask& mm4 = maskvals[level][oitr()]; oitr++;
//    const MultiMask& mm5 = maskvals[level][oitr()]; oitr++;
//#endif
//
//    const int nc = rhsL.nComp();
//
//    const bool tiling = true;
//
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
//    for (MFIter solnLmfi(solnL,tiling); solnLmfi.isValid(); ++solnLmfi)
//    {
//	const Mask& m0 = mm0[solnLmfi];
//        const Mask& m1 = mm1[solnLmfi];
//        const Mask& m2 = mm2[solnLmfi];
//        const Mask& m3 = mm3[solnLmfi];
//#if (BL_SPACEDIM > 2)
//        const Mask& m4 = mm4[solnLmfi];
//        const Mask& m5 = mm5[solnLmfi];
//#endif
//
//	const Box&       tbx     = solnLmfi.tilebox();
//        const Box&       vbx     = solnLmfi.validbox();
//        FArrayBox&       solnfab = solnL[solnLmfi];
//        const FArrayBox& rhsfab  = rhsL[solnLmfi];
//        const FArrayBox& f0fab   = f0[solnLmfi];
//        const FArrayBox& f1fab   = f1[solnLmfi];
//        const FArrayBox& f2fab   = f2[solnLmfi];
//        const FArrayBox& f3fab   = f3[solnLmfi];
//#if (BL_SPACEDIM == 3)
//        const FArrayBox& f4fab   = f4[solnLmfi];
//        const FArrayBox& f5fab   = f5[solnLmfi];
//#endif
//
//#if (BL_SPACEDIM == 2)
//        FORT_GSRB(
//            solnfab.dataPtr(), 
//            ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
//            rhsfab.dataPtr(), 
//            ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
//            f0fab.dataPtr(), 
//            ARLIM(f0fab.loVect()), ARLIM(f0fab.hiVect()),
//            m0.dataPtr(), 
//            ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
//            f1fab.dataPtr(), 
//            ARLIM(f1fab.loVect()), ARLIM(f1fab.hiVect()),
//            m1.dataPtr(), 
//            ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
//            f2fab.dataPtr(), 
//            ARLIM(f2fab.loVect()), ARLIM(f2fab.hiVect()),
//            m2.dataPtr(), 
//            ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
//            f3fab.dataPtr(), 
//            ARLIM(f3fab.loVect()), ARLIM(f3fab.hiVect()),
//            m3.dataPtr(), 
//            ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
//	    tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
//            &nc, h[level].data(), &redBlackFlag);
//#endif
//
//#if (BL_SPACEDIM == 3)
//        FORT_GSRB(
//            solnfab.dataPtr(), 
//            ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
//            rhsfab.dataPtr(), 
//            ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
//            f0fab.dataPtr(), 
//            ARLIM(f0fab.loVect()), ARLIM(f0fab.hiVect()),
//            m0.dataPtr(), 
//            ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
//            f1fab.dataPtr(), 
//            ARLIM(f1fab.loVect()), ARLIM(f1fab.hiVect()),
//            m1.dataPtr(), 
//            ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
//            f2fab.dataPtr(), 
//            ARLIM(f2fab.loVect()), ARLIM(f2fab.hiVect()),
//            m2.dataPtr(), 
//            ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
//            f3fab.dataPtr(), 
//            ARLIM(f3fab.loVect()), ARLIM(f3fab.hiVect()),
//            m3.dataPtr(), 
//            ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
//            f4fab.dataPtr(), 
//            ARLIM(f4fab.loVect()), ARLIM(f4fab.hiVect()),
//            m4.dataPtr(), 
//            ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
//            f5fab.dataPtr(), 
//            ARLIM(f5fab.loVect()), ARLIM(f5fab.hiVect()),
//            m5.dataPtr(), 
//            ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
//	    tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
//	    &nc, h[level].data(), &redBlackFlag);
//#endif
//    }
//}

 void
   PFFem::Fsmooth_jacobi (MultiFab&       solnL,
			  const MultiFab& rhsL,
			  int            level)
 {
 }

 void
   PFFem::Fapply (MultiFab&       y,
		  const MultiFab& x,
		  int             level)
 {
   int src_comp = 0;
   int dst_comp = 0;
   int num_comp = 1;
   Fapply(y,dst_comp,x,src_comp,num_comp,level);
 }


 inline amrex::Real PFFem::Phi(int p1, int p2, int n1, int n2, int i, int j)
 {
   if (p1==n1 && p2==n2) // SAME
     if (i==j) return 4./3.;
     else return 0.;
   else if (p2==n2) // EAST/WEST
     if (i==1 && j == 1) return -2./3.;
     else if (i==2 && j==2) return 1./3.;
     else return 0.;
   else if (p1==n1) // NORTH/SOUTH
     if (i==1 && j == 1) return 1./3.;
     else if (i==2 && j==2) return -2./3.;
     else return 0.;
   else if ((n1>p1 && n2>p2) || (n1<p1 && n2<p2)) // NORTHEAST/SOUTHWEST
     if (i==j) return -1./6.;
     else return -1./4.;
   else // NORTHWEST/SOUTHWEST
     if (i==j) return -1./6.;
     else 1./4.;
 }
 
 void
   PFFem::Fapply (MultiFab&       y,
		  int             dst_comp,
		  const MultiFab& x,
		  int             src_comp,
		  int             num_comp,
		  int             level)
 {
   BL_PROFILE("PFFem::Fapply()");
#if BL_SPACEDIM == 2

   const bool tiling = true;
#ifdef _OPENMP
#pragma omp parallel
#endif
   for (MFIter ymfi(y,tiling); ymfi.isValid(); ++ymfi)
     {
       const Box&       tbx  = ymfi.tilebox();
       FArrayBox&       bfab = y[ymfi];
       const FArrayBox& ufab = x[ymfi];

       for (int p1 = tbx.loVect()[0]; p1<=tbx.hiVect()[0]; p1++)
	 for (int p2 = tbx.loVect()[1]; p2<=tbx.hiVect()[1]; p2++)
	   for (int i = 0; i < BL_SPACEDIM; i++)
	     {
	       bfab(amrex::IntVect(p1,p2),i) = 0;
	       for (int n1 = p1-1; n1 <= p1+1; n1++)
		 for (int n2 = p2-1; n2 <= p2+1; n2++)
		   for (int k = 0; k < BL_SPACEDIM; k++)
		     {
		       Real K = (mu + lambda)*Phi(p1,p2,n1,n2,i,k);
		       if (i==k) K += mu*(Phi(p1,p2,n1,n2,0,0)+Phi(p1,p2,n1,n2,1,1));
		       bfab(amrex::IntVect(p1,p2),i) = K*ufab(amrex::IntVect(n1,n2),k);
		     }
	     }
     }
#endif 
 }

}