//classes defining the command line parsing for the power spectrum code

#include <sstream>
#include <string>

//command line parsing 
#include <tclap/CmdLine.h> 

struct defaults{
  const static ptrdiff_t nbins=256, ncells=512;   //default values for the number of bins and cells
};

/* Convert type T to string */
template<class T>
std::string to_string(T arg){
  std::stringstream s;
  s<<arg;
  return(s.str());
}

class parse_base{
  protected:
  void declare();
  public:

  //create input file name argument
  TCLAP::UnlabeledValueArg<std::string> this->infile;
  //create output file name argument
  TCLAP::UnlabeledValueArg<std::string> this->outfile; 
  //number of k-bins
  TCLAP::ValueArg<ptrdiff_t> *nbins;
  //number of fftw cells
  TCLAP::ValueArg<ptrdiff_t> *ncells;

  //others
};
void parse_base::declare(){

  defaults defs;  //declare structure with the default values

  TCLAP::UnlabeledValueArg<std::string> tinfile("infile", "Input catalogue name", true, "", 
    "infile <string>"); 
  this->infile = &tinfile;
  //create output file name argument
  TCLAP::UnlabeledValueArg<std::string> toutfile("outfile", std::string("Output file name.")+
    std::string("Structure: k, P(k), P(k)+noise"), true, "", "outfile <string>"); 
  this->outfile = &toutfile;
  //number of k-bins
  TCLAP::ValueArg<ptrdiff_t> tnbins("b", "number-bins", std::string("Number of k-bins of the ")+
    std::string("output power spectrum. Positive: linear binning; ")+
    std::string("negative: logarithic binning [Default: "+to_string(defs.nbins)+"]"), 
    false, defs.nbins, "int");
  this->nbins = &tnbins;
  //number of fftw cells
  TCLAP::ValueArg<ptrdiff_t> tncells("c", "number-cells", std::string("Number of cells for ")+
    std::string("FFT transform (executed with FFTWi). [Default: "+to_string(defs.nbins)+"]"),
    false, defs.ncells, "unsiged int");
  this->ncells = &tncells;
}

class parse_box: public parse_base{
  /* contains the description, version for the box case 
   * and add the desired arguments to the command line */
  std::string get_description();
  std::string get_version();
  public:
  parse_box(int argc, char* argv[]);
};
std::string parse_box::get_description(){
  //this function returns the description of the code
  std::string description = "Description";
  description += " of the code";
  return(description);
}
std::string parse_box::get_version(){
  //return the version of the code
  return("0.01");
}
parse_box::parse_box(int argc, char* argv[]){
  //constructor: create the command line object, add the options and call parse
  TCLAP::CmdLine cmd(get_description(), ' ', get_version());

  declare();
  cmd.add( infile );
  cmd.add( outfile );
  cmd.add( nbins );
  cmd.add( ncells );

  cmd.parse(argc, argv);
}

//other classes with different descriptions, version and some of the arguments

