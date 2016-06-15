/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <egull@umich.edu>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/
#ifndef HYB_SIM_MAIN
#define HYB_SIM_MAIN
#endif

#include <iomanip>
#include <alps/params/convenience_params.hpp>
#include"hyb.hpp"

hybridization::hybridization(const alps::params &parms, int crank_)
  : alps::mcbase(parms, crank_),
  crank(crank_),
  local_config(parms),
  hyb_config(parms)
{
  sanity_check(parms); //before doing anything, check whether the input parameters make sense
  show_info(parms,crank);

  //initializing general simulation constants
    nacc.resize(6);
    nprop.resize(6);
    update_type.push_back("change zero state   ");
    update_type.push_back("insert segment      ");
    update_type.push_back("remove segment      ");
    update_type.push_back("insert anti-segment ");
    update_type.push_back("remove anti-segment ");
    update_type.push_back("swap segment        ");
    sweeps=0;                                  //Sweeps currently done
  thermalization_sweeps = parms["THERMALIZATION"].as<unsigned long>();                                 //Sweeps to be done for thermalization
  total_sweeps = parms["SWEEPS"].as<unsigned long>();                                                  //Sweeps to be done in total
  n_orbitals = parms["N_ORBITALS"];                                                //number of orbitals
  sign = 1.;                                                                       //fermionic sign. plus or minus one.

  //initializing physics parameters
  beta = parms["BETA"];                                                            //inverse temperature

  //initializing updates parameters
  N_meas = parms["N_MEAS"];                                                        //number of updates per measurement
  N_hist_orders = parms["N_HISTOGRAM_ORDERS"];                                  //number of orders that are measured for the order histogram

  //initializing measurement parameters
  spin_flip = parms["SPINFLIP"];                                                //whether to perform spin-flip updates
  MEASURE_nnt = parms["MEASURE_nnt"];                                           //measure density-density correlation function in imaginary time
  MEASURE_nnw = parms["MEASURE_nnw"];                                           //measure density-density correlation function in frequency
  MEASURE_nn = parms["MEASURE_nn"];                                              //measure density-density correlation function at equal times
  MEASURE_g2w = parms["MEASURE_g2w"];                                            //measure two-particle Green function
  MEASURE_h2w = parms["MEASURE_h2w"];                                            //measure higher-order two-particle correlator
  MEASURE_freq = parms["MEASURE_freq"];                                          //measure in frequency space
  MEASURE_legendre = parms["MEASURE_legendre"];                                  //measure in legendre polynomials
  MEASURE_sector_statistics = parms["MEASURE_sector_statistics"];                //measure sector statistics
  N_w = parms["N_MATSUBARA"];                                                    //number of Matsubara frequencies for gw
  N_l = parms["N_LEGENDRE"];                                                     //number of Legendre polynomial coefficients
  N_t = parms["N_TAU"];                                                            //number of tau slices for gt
  N_nn = parms["N_nn"];                                                          //number of tau-points on which density density correlator is measured
  N_w2 = parms["N_w2"];                                                          //number of Matsubara frequency points for two-particle measurements
  N_W = parms["N_W"];                                                            //number of bosonic Matsubara frequency points for two-particle measurements
  N_w_aux = (N_w2+N_W>1 ? N_w2+N_W-1 : 0);                                         //number of Matsubara frequency points for the measurment of M(w1,w2)

  //create measurement objects
  create_measurements();

  if(crank==0){
    std::cout<<"Hybridization Expansion Simulation CT-HYB"<<std::endl;
    std::cout<<"Part of the ALPS DMFT Project"<<std::endl;
    std::cout<<"Usage requires citation of the ALPS CT-HYB paper and the ALPS paper"<<std::endl;
    std::cout<<"Refer to the documentation for more information."<<std::endl;
  }
  std::cout<<"process " << crank << " starting simulation"<<std::endl;
}

/*
void hybridization::print_statistics(std::ostream &output) {
    int tot_acc=0,cur_prec = output.precision();
    for (int i=0;i<nacc.size();i++) tot_acc += nacc[i];
    output << std::endl << "|------------- Simulation details after " << sweeps << " sweeps -----------|" << std::endl;
    output << "  Total acceptance rate = " << std::setprecision(2) << std::fixed;
    output << (((double)tot_acc)/sweeps)*100 << "%" << std::endl;
    output << "  Individual acceptance rate for update " << std::endl;
    for (int i=0;i<nacc.size();i++) {
        output << "     " << update_type[i] << " = ";
        output << std::setprecision(2) << std::fixed << (((double)nacc[i])/sweeps)*100 << "%";
        output << " (proposal rate = ";
        output << std::setprecision(2) << std::fixed << (((double)nprop[i])/sweeps)*100 << "%)" << std::endl;
    }
    output << "|-----------------------------------------------------------------|" << std::endl;
    output.unsetf(std::ios_base::fixed);
    output.precision(cur_prec);
}
*/

void hybridization::sanity_check(const alps::params &parms){
//check whether the input parameters make sense before computing
//NOTE: these checks are likely not to be complete, passing all checks does not guarantee all parameters to be meaningful!

//first check that all mandatory parameters exist
  if(!parms.exists("N_TAU")) throw std::invalid_argument("please specify the parameter N_TAU");
  if(!parms.exists("BETA")) throw std::invalid_argument("please specify parameter BETA for inverse temperature");
  if(!parms.exists("N_MEAS")) throw std::invalid_argument("please specify parameter N_MEAS for measurement interval");
  if(!parms.exists("THERMALIZATION") ||
     !parms.exists("SWEEPS") ||
     !parms.exists("N_ORBITALS") ) throw std::invalid_argument("please specify parameters THERMALIZATION, SWEEPS, and N_ORBITALS");

//check paramater that are conditionally required
  if(parms["MEASURE_freq"] && !parms.exists("N_MATSUBARA")) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies to be measured");

  if(parms["MEASURE_legendre"] && !parms.exists("N_LEGENDRE")) throw std::invalid_argument("please specify parameter N_LEGENDRE for # of Legendre coefficients to be measured");
  if(parms["MEASURE_legendre"] && !parms.exists("N_MATSUBARA")) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies");
  if(parms["MEASURE_nnt"] && !parms.exists("N_nn")) throw std::invalid_argument("please specify the parameter N_nn for # of imaginary time points for the density-density correlator");
  if(parms["MEASURE_nnw"] && !parms.exists("N_W")) throw std::invalid_argument("please specify the parameter N_W for # of bosonic frequencies for the density-density correlator");
  if(parms["MEASURE_g2w"] || parms["MEASURE_h2w"] ){
    if(!parms.exists("N_w2") ) throw std::invalid_argument("please specify the parameter N_w2 for # of fermionic Matsubara frequencies for two-particle functions");
    if(!parms.exists("N_W") ) throw std::invalid_argument("please specify the parameter N_W for # of bosonic Matsubara frequencies for two-particle functions");
    if((int)parms["N_w2"]%2!=0) throw std::invalid_argument("parameter N_w2 must be even");
  }
  if(parms["COMPUTE_VERTEX"]){
    if( !(parms["MEASURE_freq"]) ) throw std::invalid_argument("frequency measurement is required for computing the vertex, please set MEASURE_freq=1");

    if(! (parms["MEASURE_g2w"] || parms["MEASURE_h2w"] ) ) throw std::invalid_argument("at least one two-particle quantity is required for computing the vertex, set MEASURE_g2w=1 or MEASURE_h2w=1");
    if((int) parms["N_MATSUBARA"] < ((int)parms["N_w2"]/2 + (int)parms["N_W"] - 1) ) throw std::invalid_argument("for computing the vertex, N_MATSUBARA must be at least N_w2/2+N_W-1");
  }
    VERBOSE = (parms["VERBOSE"]);

return;
}


void hybridization::show_info(const alps::params &parms, int crank){
  if(!(parms["VERBOSE"])) return;

//provide info on what is measured and how long the simulation will run
if(!crank){
  std::cout << "measuring gt" << std::endl;
  if(parms["MEASURE_freq"]) std::cout << "measuring gw" << std::endl << "measuring fw" << std::endl;
  if(parms["MEASURE_legendre"]) std::cout << "measuring gl" << std::endl << "measuring fl" << std::endl;
  if(parms["MEASURE_g2w"]) std::cout << "measuring g2w" << std::endl;
  if(parms["MEASURE_h2w"]) std::cout << "measuring h2w" << std::endl;
  if(parms["MEASURE_nn"]) std::cout << "measuring nn" << std::endl;
  if(parms["MEASURE_nnt"]) std::cout << "measuring nnt" << std::endl;
  if(parms["MEASURE_nnw"]) std::cout << "measuring nnw" << std::endl;
  if(parms["MEASURE_sector_statistics"]) std::cout << "measuring sector statistics" << std::endl;
  if(parms["COMPUTE_VERTEX"]) std::cout << "vertex will be computed" << std::endl;
  if(parms.exists("RET_INT_K")) std::cout << "using retarded interaction" << std::endl;
  if(parms.exists("U_MATRIX")) std::cout << "reading U matrix from file " << parms["U_MATRIX"] << std::endl;
  if(parms.exists("MU_VECTOR")) std::cout << "reading MU vector from file " << parms["MU_VECTOR"] << std::endl;
  std::cout << "Simulation scheduled to run " << parms["timelimit"] << " seconds" << std::endl << std::endl;
}
return;
}


std::ostream &operator<<(std::ostream &os, const hybridization &hyb){
  os<<cred<<"-----------------------------------------------------------------------------------"<<cblack<<std::endl;
  os<<hyb.local_config<<std::endl;
  os<<hyb.hyb_config<<std::endl;
  os<<cred<<"-----------------------------------------------------------------------------------"<<cblack<<std::endl;
  return os;
}
std::ostream &operator<<(std::ostream &os, const segment &s){
  os<<"( "<<s.t_start_<<" , "<<s.t_end_<<" ) ";
  return os;
}

void hybridization::define_parameters(parameters_type & parameters) {
    // If the parameters are restored, they exist
    if (parameters.is_restored()) {
        return;
    }
    // Adds the parameters of the base class
    alps::mcbase::define_parameters(parameters);
    // Adds the convenience parameters (for save/load)
    alps::define_convenience_parameters(parameters);
    // Define ct-hyb related parameters
    parameters
        .description("hybridization expansion simulation")
        .define<bool>("ACCURATE_COVARIANCE", false, "TODO: UNDERSTAND WHAT THIS DOES")
        .define<std::string>("BASEPATH","", "path in hdf5 file to which results are stored")
        .define<double>("BETA", "inverse temperature")
        .define<bool>("COMPUTE_VERTEX", false, "whether to compute the vertex functions or not.")
        .define<std::string>("DELTA","path for hybridization function file")
        .define<bool>("DELTA_IN_HDF5",false,"true if hybridization function file is in hdf5 format")
        .define<bool>("DMFT_FRAMEWORK",false,"true if we need to tie into a dmft framework")
        .define<bool>("GLOBALFLIP", false, "TODO: UNDERSTAND WHAT THIS DOES.")
        .define<double>("J",0,"interaction value for density-density Hund's coupling term J.")
        .define<bool>("K_IN_HDF5",false,"set to true if retarded interaction K is stored in hdf5.")
        .define<bool>("MEASURE_freq",true, "measure in frequency domain")
        .define<bool>("MEASURE_g2w",false, "measure two-particle Green's function in frequency space")
        .define<bool>("MEASURE_h2w",false, "measure two-particle H Green's function in frequency space")
        .define<bool>("MEASURE_legendre",false, "measure legendre Green's function coefficients")
        .define<bool>("MEASURE_nn",false, "measure static density-density correlation functions")
        .define<bool>("MEASURE_nnt",false, "measure density-density correlation functions <n(0) n(t)>")
        .define<bool>("MEASURE_nnw",false, "measure density-density correlation functions in frequency domain")
        .define<bool>("MEASURE_sector_statistics",false, "measure sector statistics")
        .define<bool>("MEASURE_time",false, "measure in the time domain")
        .define<double>("MU", "chemical potential / orbital energy values")
        .define<std::string>("MU_VECTOR", "file name for file with chemical potential / orbital energy values")
        .define<bool>("MU_IN_HDF5", false,"true if the file MU_VECTOR points to a hdf5 file")
        .define<int >("N_HISTOGRAM_ORDERS",200, "orders for the histograms of probability per order")
        .define<int >("N_LEGENDRE",0,"number of legendre coefficients")
        .define<int >("N_MATSUBARA",40,"number of matsubara coefficients")
        .define<int >("N_MEAS","number of updates per measurement")
        .define<int >("N_ORBITALS","number of spin-orbitals (sometimes called flavors)")
        .define<int >("N_TAU","number of imaginary time discretization points")
        .define<int >("N_W",0,"number of bosonic Matsubara frequencies")
        .define<int >("N_nn",0,"number of points for the measurement of the density density correlator")
        .define<int >("N_w2",0,"number of fermionic frequencies for the two-particle measurement")
        .define<std::string>("RET_INT_K","file with the retarted interaction information. See doc for format.")
        .define<bool>("SPINFLIP",false,"TODO: UNDERSTAND THIS PARAMETER")
        .define<unsigned long>("SWEEPS","total number of Monte Carlo sweeps to be done")
        .define<bool>("TEXT_OUTPUT","if this is enabled, we write text files in addition to hdf5 files")
        .define<unsigned long>("THERMALIZATION","thermalization steps")
        .define<double>("U","interaction value. Only specify if you are not reading an U matrix")
        .define<double>("Uprime",0,"interaction value Uprime. Only specify if you are not reading an U matrix")
        .define<std::string>("U_MATRIX","file name for file that contains the interaction matrix")
        .define<bool>("UMATRIX_IN_HDF5",false,"true if we store the U_matrix as /Umatrix in a hdf5 file")
        .define<bool>("VERBOSE",false,"how verbose the code is. true = more output")
        ;
}
