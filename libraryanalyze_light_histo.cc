#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>  

#include "library_access.h"
#include "utility_functions.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "THStack.h"
#include "TColor.h"
#include "TLegend.h"
#include "TMarker.h"

using namespace std;

//SBND active volume 112 Ton --> 0.0028 of DUNE --> ~2.8 SN in SBND -->1 tpc so 1.4?
bool sort_function(std::pair<double, int> pair1, std::pair<double, int> pair2)
  { return (pair1.first < pair2.first); }


int main()
{
  gRandom->SetSeed(0);


  const double MassE = 0.510998910; // mass electron - MeV/c^2
  const double Q = 0.565;//Q value of decay - Ar39
  const double t_singlet = 0.000000006; //6ns
  const double t_triplet = 0.0000015; //1.5 us
  const double scint_time_window = 0.00001; //10 us
  //const double scint_time_window = 0.000005; //5 us
  const double Eav = 20.;
  const int scint_yield = 24000;
  const double quantum_efficiency = 0.2;//Expected value
  
  vector<double> timing;
  vector< std::pair< double, int > > total_timing;

  ///Cosmic Energy Calculation
  const double av_dedx = 2; /// MeV/cm for mip

  const int f_muon = 100; //m^2/s
  const double lariat_drift = 0.47; //m
  const double lariat_beam = 0.9; //m
  const int time_window = 60; //s, the beam fires every 60 seconds
  const double cosmic_pathlen = 40; //cm

  double t_cosmic_eng = f_muon * lariat_drift * lariat_beam * time_window * cosmic_pathlen * av_dedx; ///In MeV per minute
  cout << "Cosmic Energy: " << t_cosmic_eng << " MeV" << endl;

  ///Beam Energy Calculation
  const int ntrig = 30; //roughly /spill
  const double ntrack = 2.7;
  const double trklen = 35; //cm
  double t_beam_eng = ntrig * ntrack * trklen * av_dedx;
  cout << "Beam Energy: " << t_beam_eng << " MeV" << endl; /// The beam comes in a 4 second spill every minute

  ///EVENTS...
  const int scaling = 10000.;
  //const int max_events = ((t_cosmic_eng + t_beam_eng) * scint_yield)/scaling;
  const int max_events = ((t_cosmic_eng + t_beam_eng) * scint_yield)/scaling;
  cout << "Total events: " << max_events << endl;
  cout << "With: " << scaling << " gamma / event" << endl;
  cout << "Total Gammas: " << max_events*scaling << " in: " << time_window << " s" << endl;

  LibraryAccess lar_light;

 //*********** loading the library ****************
  
  //std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160111/opLibArrayPMTswCathd.root";
  //std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160418/timing_files/Lib154PMTs8inch_FullFoilsTPB.root";
  //std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160418/timing_files/Lib154PMTs8inch_NoCathodeNoFoils.root";
  //std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160418/timing_files/Lib154PMTs8inch_OnlyCathodeTPB.root";
  std::string libraryfile = "/pnfs/lariat/scratch/users/pkryczyn/TESTFullLariatPhotonLibrary1cm-10k.root";  


  //fNVoxels = int(gxSteps*gySteps*gzSteps);
  int fNVoxels = (40*90);///Number of voxels on the cathode: #y * #z
  //fNOpChannels = 308;
  bool reflected = true;
  bool reflT = true;
  lar_light.LoadLibraryFromFile(libraryfile, reflected, reflT);

 //*********** end loading the library *************
/*
 // Reading out positions of PMT from txt file *****
  ifstream myfile;
  myfile.open("/home/colton/PhD/Year1/Practice/20160418/timing_files/posPMTs_setup1.txt");
  if(myfile.is_open()){cout << "File opened successfully" << endl;}
  vector<vector<double>> myfile_data;
  while(!myfile.eof()){
    double num_pmt, x_pmt, y_pmt, z_pmt;
    if(myfile >> num_pmt >> x_pmt >> y_pmt >> z_pmt){
      vector<double> line_data({num_pmt, x_pmt, y_pmt, z_pmt});
      //cout << line_data.at(0) << " " << line_data.at(1) << " " << line_data.at(2) << " " << line_data.at(3) << endl;
      myfile_data.push_back(line_data);
    }
    else{break;}
  }
 // End Reading out positions of PMT from txt file
*/

  ///Voxels at the cathode!!!
  std::ifstream voxel_list;
  //voxel_list.open("/home/colton/PhD/Year1/Practice/20160418/timing_files/ion_test/voxel_list.txt");
  voxel_list.open("voxel_list.txt");
  if(voxel_list.is_open()){cout << "File opened successfully" << endl;}
  vector<int> voxel_list_data;
  while(!voxel_list.eof()){
    int voxel_number;
    if(voxel_list >> voxel_number){
      voxel_list_data.push_back(voxel_number);
    }
    else{break;}
  }



//-----------------Ar39 Loop-----------------------------
  //int realisticPMT_IDs[60] = {0, 4, 8, 12, 16, 20, 24, 32, 40, 44, 48, 52, 56, 60, 64, 88, 92, 96, 100, 104, 108, 112, 120, 128, 132, 136, 140, 144, 148, 152, 154, 158, 162, 166, 170, 174, 178, 186, 194, 198, 202, 206, 210, 214, 218, 242, 246, 250, 254, 258, 262, 266, 274, 282, 286, 290, 294, 298, 302, 306};

  const int num_pmts = 5;

  vector<int> vuv_counter;
  vuv_counter.resize(num_pmts);
  for(int k = 0; k < num_pmts; k++){vuv_counter.at(k) = 0;}

  vector<int> vis_counter;
  vis_counter.resize(num_pmts);
  for(int k = 0; k < num_pmts; k++){vis_counter.at(k) = 0;}

  for(int events = 0; events < max_events; events++){

    ///8000 is the number of cathode voxels for SBND
    ///USE fNVoxels
    int rand_num = gRandom->Uniform(fNVoxels);
    int rand_voxel = voxel_list_data.at(rand_num);

    if(events % 50000 == 0){
      cout << "Running... " << events << endl;
    }

    vector<vector<double>> pmt_hits = lar_light.PhotonLibraryAnalyzer(1, scint_yield, quantum_efficiency, rand_voxel);
    
    for(int pmt_loop = 0; pmt_loop < num_pmts; pmt_loop++){
      int num_pmt = pmt_loop;
      //int num_pmt = realisticPMT_IDs[pmt_loop]; ///Needed for SBND
      int events_counter = 0;
        
      int num_VUV = pmt_hits.at(num_pmt).at(2);
      int num_VIS = pmt_hits.at(num_pmt).at(3);

      //if(num_VUV != 0 && num_VIS != 0){cout << num_VUV << ", " << num_VIS << endl;}

      vuv_counter.at(pmt_loop) += num_VUV;
      vis_counter.at(pmt_loop) += num_VIS;

    }//end looping pmts
  }//end loop over events

  int tot_vuv = 0;
  int tot_vis = 0;

  for(int pmts = 0; pmts < num_pmts; pmts++)
  {
    //cout << "PMT: " << realisticPMT_IDs[pmts] << ", VUV: " << vuv_counter.at(pmts) << ", VIS: " << vis_counter.at(pmts) << endl;
    tot_vuv += vuv_counter.at(pmts);
    tot_vis += vis_counter.at(pmts);
  }

  
  cout << "TOTALs - VUV: " <<  tot_vuv << ", VIS: " << tot_vis << " for " << time_window << " s" << endl;
  cout << "VUV: " << tot_vuv/time_window << ", VIS: " << tot_vis/time_window << " Hz" << endl;
/*
  std::sort(total_timing.begin(), total_timing.end(), sort_function);

  energy_ar_list.clear();
  scint_time_list.clear();

  int counter = 0; /// Unique PMTs
  int prev_pmt1;
  int prev_pmt2;
  int flag = 0;
  int trigger_flag = 0;
  int new_gamma1;
  bool triggered = false;
  vector<int> pmt_list;

  for(int gamma1 = 0; gamma1 < total_timing.size()-1; gamma1++)
  {
    ///If previous gamma1 saw a trigger, move up 100ns
    if(triggered == true){gamma1 = new_gamma1;}
    triggered = false;

    double t01 = total_timing.at(gamma1).first;
    for(int gamma2 = gamma1+1; gamma2 < total_timing.size(); gamma2++)
    {
      double t02 = total_timing.at(gamma2).first;
      if(total_timing.at(gamma1).second != total_timing.at(gamma2).second)
      {
	if(t02 - t01 <= 0.01)
	{
	  if(counter == 0)
	  {
	    counter++;
	    pmt_list.push_back(gamma1);
	    pmt_list.push_back(gamma2);
	  }
	  if(counter == 1 && total_timing.at(gamma2).second != prev_pmt1 && total_timing.at(gamma2).second != prev_pmt2)
	  {
	    counter++;
	    pmt_list.push_back(gamma2);
	  }/// If another unique PMT
	}/// if within time

	else
	{
	  /// Number of PMT combinations (1,2,3)
	  if(counter >= 2)
	  {
	    //cout << counter << endl;
	    for(int pmt1 = 0; pmt1 < pmt_list.size(); pmt1++)
	    {
	      int evt = pmt_list.at(pmt1); /// number in total_timing
	      int same_pmt_counter = 0;
	      for(int element = evt+1; element < total_timing.size(); element++)
	      {
		//events within 100 ns
		if(total_timing.at(element).first - total_timing.at(evt).first <= 0.1)
		{
		  if(total_timing.at(element).second == total_timing.at(evt).second)/// now see if two events are on the same PMT
		  {
		    same_pmt_counter++;
		  }
		}/// End loop if events match within 100ns
		else
		{
		  ///***** Number of PE on single PMT*****
		  if(same_pmt_counter >= 6)
		  {
		    flag++;
		    //cout << same_pmt_counter << endl;
		  }
		  break;
		}
	      }/// End loop all events again
	    }/// End loop of 3 PMTs which saw signal within 10ns
	  }/// if 3 unique PMTs
	  counter = 0;
	  pmt_list.clear();
	  break;
	}/// else

	prev_pmt1 = total_timing.at(gamma1).second;
	prev_pmt2 = total_timing.at(gamma2).second;

      }//pm1!=pmt2
    }//gamma2
    /// If all 3 PMTs satisfy the 8 PE timing
    if(flag >= 3)
    {
      trigger_flag++;
	/// Next trigger can only happen +100ns
      for(new_gamma1 = gamma1; new_gamma1 < total_timing.size()-1; new_gamma1++)
      {
	if(total_timing.at(gamma1).first + 0.1 >= time_window * 1000000.){break;}
	if(total_timing.at(new_gamma1).first >= total_timing.at(gamma1).first + 0.1)
	{
	  triggered = true;
	  break;
	} 
      }
    }
    flag = 0;
  }//gamma1
  cout << trigger_flag << endl;


  total_timing.clear();

  TCanvas *can_time = new TCanvas();
  can_time->cd();
  h_pe_time->GetXaxis()->SetTitle("Number of PE within 100ns");
  h_pe_time->Draw();
  can_time->SetLogy();
  can_time->Print("pe_time.pdf");

  TCanvas *can_time2 = new TCanvas();
  can_time2->cd();
  h_single->GetYaxis()->SetTitle("Number of Photoelectrons / 60 PMTs");
  h_single->GetXaxis()->SetTitle("Time [#mu s]");
  h_single->Draw();
  can_time2->Print("pmt_time.pdf");
*/
 //------End Ar39 loop------------------  
  return 0;

}//end main



