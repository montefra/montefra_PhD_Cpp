/*==========================================================================
 * version: 1.00, 04/10/2012
 * author: Francesco Montesano, montefra@mpe.mpg.de
 *
 * Given a file containing k and P(k) or r and xi(r)
 * compute the 1D fourier transform via the qawo algorithm from GSL
 *==========================================================================*/

#include "1Dft_corrfunc_pk.h"

#include <fstream>
#include <cmath>
#include <vector>
//command line parsing 
#include <tclap/CmdLine.h> 

//Check the value of 'ncols'. If >=2 returns true,otherwise
//false. Used as constraint in numcols.
bool larger_than_2(int ncols){
  if( ncols >= 2 )
    return( true );
  else
    return( false );
}

/*==========================================================================
 * read the input file and returns the content of the two desired columns
 * Parameters
 * ----------
 * fname: string
 *   name of the file to read
 * columns_number: int
 *   number of columns in the file
 * which_columns: vector of two int
 *   number of the columns to read
 * output
 * ------
 * x, y: vectors
 *   the two columns to read
 * file_exits: bool
 *   true if the file exists, false otherwise
 *==========================================================================*/
bool readfile(std::string fname, int columns_number, 
    std::vector<int> which_columns,
    std::vector<double> &x, std::vector<double> &y);


int main(int argc, char* argv[])
{
  try{    // Catch tclap ecceptions
    // Define the command line object, and insert a message
    // that describes the program. The second argument is the 
    // delimiter (usually space) and the last one is the version number. 
    std::string description= "Given a file containing k and P(k) ";
    description += "or r and xi(r) compute the one dimensional Fourier ";
    description += "transform via the qawo algorithm from GSL";
    TCLAP::CmdLine cmd(description, ' ', "0.01");
   
    //create input and output file name argument and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> infile("infile", "Input file name", 
	true, "", "infile", cmd); 
    TCLAP::UnlabeledValueArg<std::string> outfile("outfile", "Input file name", 
	true, "", "outfile", cmd); 

    //number of columns in the input file. Check on the number of columns
    TCLAP::FunctionConstraint<int> lt2(&larger_than_2, "greater then or equal to 2" );	
    TCLAP::ValueArg<int> numcols("c", "numcols", 
	"Number of columns in the input file [Default: 2]", false, 2, &lt2, cmd);
    //number of the columns to read
    TCLAP::MultiArg<int> columns("n", "columms", 
	"Number of the columms containing k [r] and P(k) [xi(r)]. [Default: 0,1]",
	false, "int", cmd);

    //P(k) to xi(r) or inverse transform
    std::vector<std::string> choises; //vector of choises 
    choises.push_back( "P2xi" );
    choises.push_back( "xi2P" ); 
    TCLAP::ValuesConstraint<std::string> transform_choise( choises );  //convert to tclap object
    TCLAP::ValueArg<std::string> transform( "t", "tranform", 
	"Execute the P(k) to xi(r) or the xi(r) to P(k) Fourier transform [Default: P2xi]", 
	false, "P2xi", &transform_choise, cmd);

    //limits and number of bins
    TCLAP::ValueArg<int> nbins( "b", "nbins", 
	std::string("Number of bins in the output file. Positive: linear; ")+
	std::string(" negative: logarithmic. [Default: -200]"),
	false, -200, "int", cmd);
    TCLAP::ValueArg<double> rkmin( "", "min", 
	"Minimum of r or k. [Default: 0 or 1e-4]", false,
	-999, "double", cmd);
    TCLAP::ValueArg<double> rkmax( "", "max", 
	"Maximum of r or k. [Default: 400 or 1e3]", false,
	-999, "double", cmd);

    // Parse the argv array.
    cmd.parse(argc, argv);

    /*==========================================================================
     * code itself
     *==========================================================================*/

    //check the columns to read in the input file
    if( numcols.getValue() <2 ){ 
      std::cerr << "The input file must have at least 2 columns" << std::endl;
      exit( 9 );
    }
    std::vector<int> vcolums = columns.getValue();  //vector containing the columns to read
    size_t vcolums_size = vcolums.size();   //get the number of elements
    if( vcolums_size < 2 ){
      std::cerr << "Setting the default columns to the first two: 0, 1" << std::endl;
      vcolums.clear();
      vcolums.push_back(0);
      vcolums.push_back(1);
      vcolums_size = 2;
    }
    else if( vcolums_size > 2)
      std::cerr << "Too many columns, keeping the first two" << std::endl;
    
    //read the intput file
    std::vector<double> x, y; //vectors containing x and y
    if( !readfile(infile.getValue(), numcols.getValue(), vcolums, x, y) ){
      std::cerr << "The input file " + infile.getValue() + "does not exist" << std::endl; 
      exit(10);
    }

    //check the output file
    std::fstream out(outfile.getValue().c_str(), std::fstream::in);     
    if( out.is_open() ){
      std::cerr << "The output file " + outfile.getValue() + " already exists. ";
      std::cerr << "Please chose an other file name" << std::endl;
      out.close();
      exit(11);
    }
    out.close();
    //open the output file
    out.open(outfile.getValue().c_str(), std::fstream::out);     
    out.setf(std::ios_base::scientific);
    out.precision(6);
    out.width(9);

    //create the output array
    int number_bins = nbins.getValue(); //number of bins
    double max=rkmax.getValue(), min=rkmin.getValue(); //maximum and minimum values
    //output arrays
    double *xout=new double[abs(number_bins)], *yout=new double[abs(number_bins)]; 
    //set default min and max
    if( transform.getValue().compare( choises[0] ) == 0 ){  //if P(k) to xi(r) required
      if( max < 0 ) max = 400;
      if( min < 0 ) min = 0;
    }
    else{   //if xi(r) to P(k) required
      if( max < 0 ) max = 1e3;
      if( min < 0 ) min = 1e-4;
    }
    //fill xout with the lin or log spaced values
    if( number_bins > 0 ){
      double binsize = (max-min)/(number_bins-1);
      for(int i=0; i<number_bins; ++i) xout[i] = binsize * i + min;
    }
    else{
      number_bins = -number_bins;
      double binsize = log(max/min)/(number_bins-1);
      for(int i=0; i<number_bins; ++i) xout[i] = exp(binsize * i) * min;
    }

    //do the tranformation
    if( transform.getValue().compare( choises[0] ) == 0 ){  //if P(k) to xi(r) required
      pk2corfunc(&x[0], &y[0], x.size(), xout, yout, number_bins );
    }
    else{   //if xi(r) to P(k) required
      corfunc2pk(&x[0], &y[0], x.size(), xout, yout, number_bins );
    }

    //write the result in the file
    for( int i=0; i<number_bins; ++i )
      out << xout[i] << "\t" << yout[i] << std::endl;

    out.close();  //close the output file
    x.clear();   //clear input arrays
    y.clear();
    delete[] xout;  //clear output arrays
    delete[] yout; 


  }
  catch (TCLAP::ArgException &e){  // catch any exceptions
    cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl; 
    exit(1);
  }

  exit(0);
}


/*==========================================================================
 * FUNCTION DEFINITIONS
 *==========================================================================*/
/*==========================================================================
 * read the input file and returns the content of the two desired columns
 * Parameters
 * ----------
 * fname: string
 *   name of the file to read
 * columns_number: int
 *   number of columns in the file
 * which_columns: vector of two int
 *   number of the columns to read
 * output
 * ------
 * x, y: vectors
 *   the two columns to read
 * file_exits: bool
 *   true if the file exists, false otherwise
 *==========================================================================*/
bool readfile(std::string fname, int columns_number, 
    std::vector<int> which_columns,
    std::vector<double> &x, std::vector<double> &y){
  //open the file
  std::ifstream infile(fname.c_str(), std::ifstream::in);
  if( ! infile.is_open() ) return( false );  //if the file do not exists return false

  double *atemp = new double[columns_number];  //array containing each line of the file
  for( ;; ){   //loop until the end of file
    for( int i=0; i<columns_number; ++i ) infile >> atemp[i];  //read the line
    if( infile.eof() ) break;  //break if the end is reached
    x.push_back( atemp[which_columns[0]] );  //push x and y in the output
    y.push_back( atemp[which_columns[1]] );
  }
  infile.close();  //close the file
  delete[] atemp;  //clear the pointer
  return( true );
}

