#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>

#include "GeneralAMRIntegrator.H"

#include <fstream>


using namespace amrex;

void
GeneralAMRIntegrator::CreateCleanDirectory () const
{
  amrex::UtilCreateCleanDirectory(plot_file, false);
}

std::vector<std::string>
GeneralAMRIntegrator::PlotFileName (int lev) const
{
  std::vector<std::string> name;
  name.push_back(plot_file+"/");
  name.push_back(amrex::Concatenate("", lev, 5));
  return name;
}

Array<const MultiFab*>
GeneralAMRIntegrator::PlotFileMF (int fab) const
{
  Array<const MultiFab*> r;
  for (int i = 0; i <= finest_level; ++i) {
    //r.push_back(phi_new[fab][i].get());
  }
  return r;
}

Array<std::string>
GeneralAMRIntegrator::PlotFileVarNames () const
{
  Array<std::string> names;
  // for (int n = 0; n < number_of_grains; n++)
  //   names.push_back(amrex::Concatenate("phi",n));
  names.push_back("phi_sum");
  names.push_back("boundaries");
  return names;
}

void
GeneralAMRIntegrator::WritePlotFile () const
{
  const std::vector<std::string>& plotfilename = PlotFileName(istep[0]);
  const auto& mf = PlotFileMF(0);
  const auto& varnames = PlotFileVarNames();
  amrex::WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1], finest_level+1, mf, varnames,
				 Geom(), t_new[0], istep, refRatio());

  if (ParallelDescriptor::IOProcessor())
    {
      std::ofstream outfile;
      if (istep[0]==0) outfile.open(plot_file+"/output.visit",std::ios_base::out);
      else outfile.open(plot_file+"/output.visit",std::ios_base::app);
      outfile << plotfilename[1] + "/Header" << std::endl;
}
}

void
GeneralAMRIntegrator::InitFromCheckpoint ()
{
  amrex::Abort("GeneralAMRIntegrator::InitFromCheckpoint: todo");
}
