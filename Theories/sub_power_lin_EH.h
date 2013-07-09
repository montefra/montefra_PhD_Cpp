/*==========================================================================*/
/* Version: 0.00           date:22/01/08                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */ 
/* Purpose : functions called in power_lin_EH.c                             */
/*                                                                          */
/*==========================================================================*/



/*==========================================================================*/
/* Function: read data from file:                                           */ 
/*name file in input                                                        */ 
/*vector made of structures "pow_spec" and dimention of the vector as output*/ 
/*==========================================================================*/


double input_data(string input_file, vector<pow_spec> &indata)
{
  pow_spec read;
  double size_vec;


  ifstream in(input_file.c_str());
  while (true){
    in >>  read.k >> read.P_k;
        
    if (in.eof() == true){
      break;
    }
    else {
      indata.push_back(read);
    }
  }

  in.close();

  size_vec=indata.size();
    
  return(size_vec);

}

/*==========================================================================*/
/* Eisenstein & Hu functions used: 					    */
/* http://background.uchicago.edu/~whu/power/power.html			    */
/*==========================================================================*/

double TFnowiggles(double omega0, double f_baryon, double hubble, 
		double Tcmb, double k_hmpc)
/* Input: omega0 -- CDM density, in units of critical density
	  f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
	  hubble -- Hubble constant, in units of 100 km/s/Mpc
	  Tcmb -- Temperature of the CMB in Kelvin; Tcmb<=0 forces use of
			COBE FIRAS value of 2.728 K
	  k_hmpc -- Wavenumber in units of (h Mpc^-1). */
/* Output: The value of an approximate transfer function that captures the
non-oscillatory part of a partial baryon transfer function.  In other words,
the baryon oscillations are left out, but the suppression of power below
the sound horizon is included. See equations (30) and (31).  */
/* Note: If you prefer to use wavenumbers in units of Mpc^-1, use hubble -> 1
and omega0 -> omega0*hubble^2. */ 
{
    double k, omhh, theta_cmb, k_equality, q, xx, alpha_gamma, gamma_eff;
    double q_eff, T_nowiggles_L0, T_nowiggles_C0;

    k = k_hmpc*hubble;	/* Convert to Mpc^-1 */
    omhh = omega0*hubble*hubble;
    if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
    theta_cmb = Tcmb/2.7;

    k_equality = 0.0746*omhh/SQR(theta_cmb);
    q = k/13.41/k_equality;
    xx = k*TFsound_horizon_fit(omhh, f_baryon, 1.0);

    alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
		SQR(f_baryon);
    gamma_eff = omhh*(alpha_gamma+(1-alpha_gamma)/(1+POW4(0.43*xx)));
    q_eff = q*omhh/gamma_eff;

    T_nowiggles_L0 = log(2.0*2.718282+1.8*q_eff);
    T_nowiggles_C0 = 14.2 + 731.0/(1+62.5*q_eff);
    return T_nowiggles_L0/(T_nowiggles_L0+T_nowiggles_C0*SQR(q_eff));
}


double TFsound_horizon_fit(double omega0, double f_baryon, double hubble)
/* Input: omega0 -- CDM density, in units of critical density
	  f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
	  hubble -- Hubble constant, in units of 100 km/s/Mpc*/
/* Output: The approximate value of the sound horizon, in h^-1 Mpc. */
/* Note: If you prefer to have the answer in  units of Mpc, use hubble -> 1
and omega0 -> omega0*hubble^2. */ 
{
    double omhh, sound_horizon_fit_mpc;
    omhh = omega0*hubble*hubble;
    sound_horizon_fit_mpc = 
	44.5*log(9.83/omhh)/sqrt(1+10.0*pow(omhh*f_baryon,0.75));
    return sound_horizon_fit_mpc*hubble;
}


