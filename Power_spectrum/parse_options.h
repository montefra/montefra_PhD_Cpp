/*=======================================================================*/
/* Author: Francesco Montesano (montefra@mpe.mpg.de)                     */
/*                                                                       */
/* This file contains command line parsing related objects, and defaults */
/* Should contain all that is needed for power spectra codes             */
/* include this file within the main and add to the command line only    */
/* the desired options and arguments                                     */
/*=======================================================================*/

#ifndef PARSE_OPTIONS_H
#define PARSE_OPTIONS_H

/*=======================================================================*/
/* default value for the command line arguments                          */
/*=======================================================================*/

ptrdiff_t dnbins=256, dncells=512;   //default values for the number of bins and cells
double dbsize = 1340.;  //default size of the simulation box
double dbmax = 4000., dbmin=0.;    //default maximum and minimum coordinates of the box
std::string dmas="NGP";  //default value for the mass assigment scheme
int dcorr=0;  //default MAS correction
int dnoise=1;  //default shot noise computation method
double dpw = 0.;  //default value of pw

/*=======================================================================*/
/* command line arguments                                                */
/*=======================================================================*/

//create output file name argument
TCLAP::UnlabeledValueArg<std::string> outfile("outfile", 
    "Output file name. Structure: k, P(k), P(k)+noise", true, "", "outfile (string)"); 
//create input file(s) name argument
TCLAP::UnlabeledValueArg<std::string> infile("infile", "Input catalogue name (ascii)", true, 
    "", "infile (string)"); 

TCLAP::UnlabeledMultiArg<std::string> infiles("infiles", 
    std::string("Input catalogue names. Interpretation of the files:                 ")+
    std::string("1) if number of infiles even: couples of catalogue/random files ")+
    std::string("to load on the same grids and the power spectrum;                          ")+
    std::string("2) if number of infiles odd or 'winonly' given: all files are randoms ")+
    std::string("are loaded in the same grid a one window function is computed;                          ")+
    std::string("3) if 4 infiles and 'crosspk' given: two catalogue/random couples ")+
    std::string("used to computed the cross power spectrum.                                      ")+
    std::string("In cases 1 and 3 the files must ordered as be 'cat1 ran1 [cat2 ran2 ...]'"),
    true, "infiles (string)");

TCLAP::SwitchArg winonly("w", "winonly", "The input file must be used to compute the window function");

TCLAP::SwitchArg crosspk("", "crosspk", "Compute the cross power spectrum using the 4 files given");

//verbose mode
TCLAP::SwitchArg verbose("v", "verbose", "Verbose mode");

//maximim and minimum value to use for the spherical averaged power spectrum
TCLAP::ValueArg<double> kmax("", "kmax", std::string("Maximum k where to evaluate the ")+
    std::string("spherical averaged power spectrum. [Default: Nyquist wavenumber]"), 
    false, -1., "double");
TCLAP::ValueArg<double> kmin("", "kmin", std::string("Minimum k where to evaluate the ")+
    std::string("spherical averaged power spectrum. [Default: 0]"), 
    false, 0., "double");
//number of k-bins
TCLAP::ValueArg<ptrdiff_t> nbins("b", "number-bins", 
    std::string("Number of k-bins of the output power spectrum. Positive: linear binning; ")+
    std::string("negative: logarithic binning [Default: ")+to_string(dnbins)+"]", 
    false, dnbins, "int");

//number of fftw cells
TCLAP::ValueArg<ptrdiff_t> ncells("c", "number-cells", std::string("Number of cells for ")+
    std::string("FFT transform (executed with FFTW). [Default: ")+to_string(dncells)+"]", 
    false, dncells, "unsiged int");

//grid dimensions
TCLAP::ValueArg<double> bsize("s", "box-size", std::string("Side of the cubic box of the ")+
    std::string("of the simulation. [Default ")+to_string(dbsize)+"]", false, dbsize, "double");
TCLAP::ValueArg<double> bmax("", "box-max", std::string("Maximum coordinate of the box used ")+
    std::string("to perfom the FFT. The side of the box, assumed cubic, is 'box-max'-'box-min. '") +
    std::string("[Default ")+to_string(dbmax)+"]", false, dbmax, "double");
TCLAP::ValueArg<double> bmin("", "box-min", std::string("Minimum cooridinate of the box used ")+
    std::string("to perfom the FFT. [Default ")+to_string(dbmin)+"]", false, dbmin, "double");

//choise of the mas
std::vector<std::string> choise;
choise.push_back("NGP");
choise.push_back("CIC");
choise.push_back("TSC");
TCLAP::ValuesConstraint<std::string> choisemas(choise);
choise.clear();
TCLAP::ValueArg<std::string> mas("", "mas", "Mass assigment scheme [Default: "+dmas+"]", false, 
    dmas, &choisemas);
//choise of the correction of the mas
std::vector<int> ichoise;
ichoise.push_back(0);   //no correction
ichoise.push_back(1);   //correct by W(k)
ichoise.push_back(2);   //correct by sum_n W(k+2nk_N)
TCLAP::ValuesConstraint<int> choisecorr(ichoise);
ichoise.clear();
TCLAP::ValueArg<int> corr("", "correction", std::string("MAS correction. 0: no correction; ")+
    std::string("1: correct by dividing each mode by W(k), the Fourier transform of the MAS; ")+
    std::string("2: correct by deviding each mode by sum_n W(k+2*n*k_N), with k_N the Nyquist ")+
    std::string("wavenumber. [Default: ")+to_string(dcorr)+"]", false, dcorr, &choisecorr);

//pw to be used in creating the fkp or pvp weights
TCLAP::ValueArg<double> pw("", "pw", std::string("Value of pw used to compute the FKP or ")+
   std::string("the PVP weight. [Default: ") + to_string(dpw) + "]", false, dpw, "PW (double)");

//compute shot noise method
ichoise.push_back(1);
ichoise.push_back(2);
TCLAP::ValuesConstraint<int> choisenoise(ichoise);
ichoise.clear();
TCLAP::ValueArg<int> whichnoise("", "noise", 
    std::string("choise between two methods to compute the noise [Default: ")+to_string(dnoise)+"]          "+
    std::string("1: noise = (alpha+1)*alpha*sum(w^2_{random})/N^2 (as in FKP 1994)         ")+
    std::string("2: noise = sum(w^2_{cat}) + alpha^2*sum(w^2_{random})/N^2 (as in Cole et al. 2005)"),
    false, dnoise, &choisenoise);

//interpret the input file name as a root of a set of binary files
std::string descr = "N, positive";
TCLAP::FunctionConstraint<int> constrpos(positive, descr);
TCLAP::ValueArg<int> dm("", "dm", std::string("If selected, assumes that the input catalogue ")+
   std::string("is in 'N' Gadget2 binary files: 'infile'+i, with i=0,'N'-1."), 
   false, 0, &constrpos);

//redshift space distortions
ichoise.push_back(0);   //x direction
ichoise.push_back(1);   //y direction
ichoise.push_back(2);   //z direction
TCLAP::ValuesConstraint<int> choisezdist(ichoise);
ichoise.clear();
TCLAP::ValueArg<int> zdist("z", "z-distortion", std::string("Enable redshift ")+
    std::string("along one axis of the simulation box. 0: x axis; 1: y axis; 2: z axis."), 
    false, -1, &choisezdist);

//log transform
TCLAP::ValueArg<double> epsilon("l", "log-transform", std::string("Enable ln(1+delta) transform. ")+
    std::string("If negative, do ln(1+delta), if delta>0, delta, otherwise. ")+
    std::string("If zero, do ln(1+delta). ")+ 
    std::string("If positive, do delta0*ln(1+delta/delta0), with delta0 = fabs( min(delta) )")+
    std::string("+epsilon. Epsilon should be a small number"), false, -1, "EPSILON (double)");

//mass limits
TCLAP::MultiArg<double> mlim("m", "mass-limits", std::string("If not used the power spectrum ")+
    std::string("is compute for the full catalogue. If called only once, the value is used as ")+
    std::string("lower mass limit. If called twice the power spectrum is computed in the given ")+
    std::string("mass bins. If called four times the cross power spectrum between two mass bins ")+
    std::string("is computed."), false, "double" );

#endif
