/*==========================================================================*/
/* Version: 1.00           date:29/02/08                        	    */
/*                                                              	    */
/* Author: Francesco Montesano, MPE, Garching                    	    */
/*                                                        	            */ 
/* Purpose : Calculation of the linear power spectrum using the functions   */
/* by Heisenstein & Hu [Baryonic Features in the Transfer Function: (also at*/
/* astro-ph/9709112); Power Spectra for Cold Dark Matter and Its Variants:  */
/* (also at astro-ph/9710252)]					            */
/* 									    */
/* http://background.uchicago.edu/~whu/power/power.html	                    */
/*==========================================================================*/

/*==========================================================================*/
/* Input: cosmological parameter, specified in the following  		    */
/* Output: Linear power spectrum without wiggles normalised al large scale  */
/*         with the input one                                               */
/*==========================================================================*/



#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<iomanip>
#include <math.h>
#include <sstream>

using namespace std;

#define PI 3.14159265
#define n_mean 20
#define DEBUGin 0
#define Building 0

/* Convenience from Numerical Recipes in C, 2nd edition */
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static double cubearg;
#define CUBE(a) ((cubearg=(a)) == 0.0 ? 0.0 : cubearg*cubearg*cubearg)
static double pow4arg;
#define POW4(a) ((pow4arg=(a)) == 0.0 ? 0.0 : pow4arg*pow4arg*pow4arg*pow4arg)

typedef struct power_spectrum
{ double k;     // h/Mpc
  double P_k;
} pow_spec;

/*function that read data from file*/
double input_data(string input_file, vector<pow_spec> &indata);

/* Eisenstein & Hu functions*/
double TFnowiggles(double omega0, double f_baryon, double hubble, double Tcmb, double k_hmpc);
double TFsound_horizon_fit(double omega0, double f_baryon, double hubble);

/*==========================================================================*/
/*                               MAIN                                       */ 
/*==========================================================================*/

int main()
{
  int i;

  string input, output;

  vector<pow_spec> in_power;
  vector<pow_spec> noweeggles_power;
  vector<pow_spec> vector_ratio;
  pow_spec load;

  double n_rows_in;  // number of point in the infile

  double Omega_m,      // Relative matter density: CDM+baryons
         Omega_b,      // Relative baryon density
	 hubble_r,     // Hubble constant in unit of 100 km/(s*Mpc) 
	 T_CMB,        // CMB temperature
	 n_spect,      // spectral index
	 baryon_f;     // fraction of baryons with respect to CDM

  /*normalization on Power Spectrum*/
  double sum, A;

  /*Input cosmological parameters*/

  cout << endl;
  cout << "Please enter the cosmological parameter to create the power spectrum"
  << endl;

  cout << "Hubble constant in unit of 100 km/(s*Mpc): [if <=0 default value from Ariel et al, MNRAS 366,189(2006), tab3, col b6] ";
  cin >> hubble_r;
  if(hubble_r<=0.) hubble_r = 0.735;

  cout << "Relative matter density (CDM+barions): [if <=0 default value from Ariel et al, MNRAS 366,189(2006), tab3, col b6] ";
  cin >> Omega_m;
  if(Omega_m<=0.) Omega_m = 0.237;

  cout << "Relative barion density: [if <=0 default value from Ariel et al, MNRAS 366,189(2006), tab3, col b6] ";
  cin >> Omega_b;
  if(Omega_b<=0.) Omega_b = 0.0225/pow(hubble_r,2.);

  cout << "CMB temperature: [if <=0 default value from COBE: T_CBM=2.728] ";
  cin >> T_CMB;
  if(T_CMB <= 0.) T_CMB = 2.728;

  cout << "Spectral index: [if <=0 default value from Ariel et al, MNRAS 366,189(2006), tab3, col b6] ";
  cin >> n_spect;
  if(n_spect<=0.) n_spect = 0.954;

#if DEBUGin
  cout << Omega_m << " " << Omega_b << " " << hubble_r << " " << T_CMB << " " << n_spect << endl;
#endif

  /*input power spectrum*/

#if !Building
  cout << endl;
 try_again_inpow:
  cout << "Insert the name of the file for the input linear power spectrum:\n";
  cin >> input;

  ifstream in(input.c_str());
  if (in.is_open())
    {
      cout << "The file '" << input << "' exists." << endl;
      in.close();
    }
  else {
    cerr << "Invalide file name, try again!\n";
    goto try_again_inpow;
    }
  cout << endl;

  /* output file name and checking*/

 try_again_out1:
  cout << "Please insert the output file name for the linear power spectrum without wiggles:" << endl;
  cin >> output;
  cout << endl;

  ifstream chout_file(output.c_str());
  if(chout_file.is_open()){
    cout << "The file already exists. Choose an other name or delete it" << endl;
    chout_file.close();
    goto  try_again_out1;
  }
#endif

#if Building 
  input = "../Power_spectra/b6_matterpower_z0.dat";

  ostringstream s1;
  s1 << "../Power_spectra/Debug_wiggles/power_nowiggles_" << n_mean << "_z0.dat";
  output = s1.str();
#endif

  ofstream out(output.c_str());

  n_rows_in = input_data(input, in_power);

#if DEBUGin
  ofstream read_data("../Power_spectra/Debug_wiggles/read_data.d");
  read_data.setf(ios_base::scientific);
  read_data.precision(6);
  read_data.width(9);
  for(i=0; i<int(n_rows_in); i++) read_data << in_power[i].k << "\t" << in_power[i].P_k << endl;
#endif

  cout << "input file read" << endl;
  cout << endl;

  /*========================================================================*/
  /* calculation of non wiggles power Spectrum                              */
  /*========================================================================*/

  baryon_f = Omega_b/Omega_m;

  for(i=0; i<int(n_rows_in); i++){
    load.P_k = TFnowiggles(Omega_m, baryon_f, hubble_r, T_CMB, in_power[i].k);

    load.k = in_power[i].k;
    load.P_k = pow(load.k, n_spect) * pow(load.P_k, 2.);

    noweeggles_power.push_back(load);

    load.P_k = in_power[i].P_k / load.P_k;

    vector_ratio.push_back(load);
  }

  /*========================================================================*/
  /* normalization of non wiggles power Spectrum                            */
  /* equal to the input Power spectrum at large scales                      */
  /*========================================================================*/

  sum = 0.;
  for(i=0; i<n_mean; i++) sum = sum + vector_ratio[i].P_k;
  A = sum / n_mean;

  for(i=0; i<int(n_rows_in); i++) out << noweeggles_power[i].k << "\t" << A*noweeggles_power[i].P_k << endl;

  return 0;
}

#include "sub_power_lin_EH.h"
