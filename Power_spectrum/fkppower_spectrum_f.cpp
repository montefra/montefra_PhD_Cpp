/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/* Definition of the functions                                              */
/*==========================================================================*/

#include "fkppower_spectrum.h"

/*=======================================================================*/
/* Code specific functions definitions                                   */
/*=======================================================================*/
std::string get_description()
/* description of the code   */
{
  std::string description = "This code computes the angular averaged power spectrum ";
  description += "from catalogues with a realistic geometry using the FKP estimator.\n";
  description += "The input file must have the following structure:\n";
  description += "\tx y z w_fkp (w_fc+w_rf-1) w_systematics n(z) z column\n";
  description += "(the last column is ignored)";
  return(description);
}
std::string get_version()
/*version of the code */
{ return("0.01"); }

/*=======================================================================
 * Print error, finalise mpi and exit
 * Parameters
 * ----------
 *  message: error message to print to std err
 *  code: error code
 *  myrank: rank of the processor
 *  root: root processor
 *  com: MPI communicator
 *=======================================================================*/
void on_error(std::string message, int error_code, int myrank, int root,
    MPI_Comm com){
  if(myrank == root)
    std::cerr << message << std::endl;
  MPI_Barrier(com);
  MPI_Finalize();
  exit(error_code);
}

/*==========================================================================
 * Do the checks when the cross power spectrum is required
 * if files.size > 4, files is resized to 4
 * Parameters
 * ----------
 * files: vector of strings with the input files
 * winonly: if the window function required
 * cross: if the cross power spectrum or window function required
 * myrank: rank of the processor
 * root: root processor
 * com: MPI communicator
 *==========================================================================*/
void check_input_files(std::vector<std::string> &files, bool winonly, 
    bool cross, bool two_grids, int myrank, int root, MPI_Comm com){
  size_t ninfiles = files.size(); //number of files

  if(winonly == true){
    if(ninfiles!=1 && ninfiles!=2)
      on_error("Only one or two files are required for the window function",
          1, myrank, root, com);
    if((cross == true || two_grids == true) && ninfiles !=2)
      on_error("To compute the cross window function or to use two grids two files are needed",
          2, myrank, root, com);
  }
  else {
    if(ninfiles!=2 && ninfiles!=4)
      on_error("Only two or four files are required for the power spectrum",
          3, myrank, root, com);
    if((cross == true || two_grids == true) && ninfiles !=4)
      on_error("To compute the cross power spectrum or to use two grids two files are needed",
          4, myrank, root, com);
  }
}

/*==========================================================================
 * Create the header for the output file
 * Parameters
 * ----------
 *  sumscat1: sums of the first catalogue
 *  sumscat2: sums of the second catalogue
 *  sumsran1: sums of the first random
 *  sumsran2: sums of the second random
 *  dimsums: number of elements in the above arrays
 *  nfiles: number of input files
 *  winonly: window function or power spectrum
 * output
 * ------
 *  header: string with the full header
 *==========================================================================*/
std::string create_header(double *sumscat1, double *sumscat2, double *sumsran1,
    double *sumsran2, int dimsums, bool two_grids, bool cross, bool winonly){
  bool twodata=false;  //wether there are two data lines to print
  int random=0;  //add randoms to the header: 0: no randoms; 1: one random; 2: two randoms
  std::string data1, data2, random1, random2; //strings containing the data and randoms strings
  //decide what to print
  if(two_grids==true || cross==true){
    twodata=true;
    data1="#data1";
    data2="#data2";
  }
  else data1="#data";  //data2 not needed
  if(winonly==false){
    if(two_grids==true || cross==true){
      random=2;
      random1="#random1";
      random2="#random2";
    }
    else{
      random=1;
      random1="#random";
    }
  }

  //create the header
  std::string header = "# \t sum(w) \t sum(w^2n(z)) \t sum(w^2)\n";
  header += data1;
  for(int i=0; i<dimsums; ++i) header += "\t"+to_string(sumscat1[i], 4);
  header += "\n";
  if(twodata==true){
    header += data2;
    for(int i=0; i<dimsums; ++i) header += "\t"+to_string(sumscat2[i], 4);
    header += "\n";
  }
  if(random>0){
    header += random1;
    for(int i=0; i<dimsums; ++i) header += "\t"+to_string(sumsran1[i], 4);
    if(random==2){
      header += "\n";
      header += random2;
      for(int i=0; i<dimsums; ++i) header += "\t"+to_string(sumsran2[i], 4);
    }
  }
  return(header);
}

/*==========================================================================
 * class created to read all the files. Provides a simple interface,
 * a constructor (inline) that allocates arrays used while reading the files
 *   and saves all the common information to avoid passing them each time
 * a function that read the files and save the content into the desired grid
 *==========================================================================*/
/*==========================================================================
 * Read the input ascii file and fill the grid
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * zrange: array of two float: only objects within the redshift range considered
 * zrange: range of redshift
 * ignorew: vector<int> (-1,1,2): if !=-1 set to 1 fl.w[1] or fl.w[2]
 * repeatw: bool: if true, assign object fl.w[1] times with fl.w[1]=1
 * output
 * ------
 * sums: array containing sum(w), sum(n*w^2), sum(w^2)
 *==========================================================================*/
double *readfiles::read_file(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
    std::vector<double> zrange, std::vector<int> ignorew, bool repeatw){
  double *sums = init_sums(); //sum(w), sum(n*w^2), sum(w^2)

  std::ifstream in(ifile.c_str());  //open the input file 
  
  if(zrange[1]<=0 && ignorew[0]==-1 && repeatw==false) 
    read_asitis(in, grid, sums);
  else if(zrange[1]>0 && ignorew[0]==-1 && repeatw==false) 
    read_zrange(in, grid, sums, zrange);
  else if(zrange[1]<=0 && ignorew[0]!=-1 && repeatw==false) 
    read_ignorew(in, grid, sums, ignorew);
  else if(zrange[1]<=0 && ignorew[0]==-1 && repeatw==true) 
    read_repeatw(in, grid, sums);
  else if(zrange[1]>0 && ignorew[0]!=-1 && repeatw==false) 
    read_zrange_ignorew(in, grid, sums, zrange, ignorew);
  else if(zrange[1]>0 && ignorew[0]==-1 && repeatw==true) 
    read_zrange_repeatw(in, grid, sums, zrange);
  else if(zrange[1]<=0 && ignorew[0]!=-1 && repeatw==true) 
    read_ignorew_repeatw(in, grid, sums, ignorew);
  else
    read_zrange_ignorew_repeatw(in, grid, sums, zrange, ignorew);

  in.close();  //close it

  return(sums);  //return the number of haloes assigned to the grid
}

/*==========================================================================
 * read in the file with the various options. 
 * For argument refer to the description of 'read_file'
 *==========================================================================*/
void readfiles::read_asitis(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid,
    double* sums){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_zrange(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid,
    double* sums, std::vector<double> zrange){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    if(zrange[0]>fl.z || zrange[1]<fl.z) continue;  //read the next line
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_ignorew(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid,
    double* sums, std::vector<int> ignorew){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    for(size_t i=0; i<ignorew.size(); ++i) fl.w[ignorew[i]] = 1.;
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_repeatw(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid,
    double* sums){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    int nw = fl.w[1];
    fl.w[1] = 1.;
    for(int i=0; i<nw; ++i){
      sums_weights(sums);
      grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
          fl.w[0]*fl.w[1]*fl.w[2]);
    }
  }
}
void readfiles::read_zrange_ignorew(std::ifstream &inif, ps_r2c_c2r_mpi_inplace
    &grid, double* sums, std::vector<double> zrange, std::vector<int> ignorew){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    if(zrange[0]>fl.z || zrange[1]<fl.z) continue;  //read the next line
    for(size_t i=0; i<ignorew.size(); ++i) fl.w[ignorew[i]] = 1.;
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_zrange_repeatw(std::ifstream &inif, ps_r2c_c2r_mpi_inplace
    &grid, double* sums, std::vector<double> zrange){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    if(zrange[0]>fl.z || zrange[1]<fl.z) continue;  //read the next line
    int nw = fl.w[1];
    fl.w[1] = 1.;
    for(int i=0; i<nw; ++i){
      sums_weights(sums);
      grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
          fl.w[0]*fl.w[1]*fl.w[2]);
    }
  }
}
void readfiles::read_ignorew_repeatw(std::ifstream &inif, ps_r2c_c2r_mpi_inplace
    &grid, double* sums, std::vector<int> ignorew){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    for(size_t i=0; i<ignorew.size(); ++i) fl.w[ignorew[i]] = 1.;
    int nw = fl.w[1];
    fl.w[1] = 1.;
    for(int i=0; i<nw; ++i){
      sums_weights(sums);
      grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
          fl.w[0]*fl.w[1]*fl.w[2]);
    }
  }
}
void readfiles::read_zrange_ignorew_repeatw(std::ifstream &inif,
    ps_r2c_c2r_mpi_inplace &grid, double* sums, std::vector<double> zrange, 
    std::vector<int> ignorew){
  for(;;){  //loop through the lines of the file
    read_line(inif);
    if(inif.eof() == true) break;    //break if the end of the file reached
    if(zrange[0]>fl.z || zrange[1]<fl.z) continue;  //read the next line
    for(size_t i=0; i<ignorew.size(); ++i) fl.w[ignorew[i]] = 1.;
    int nw = fl.w[1];
    fl.w[1] = 1.;
    for(int i=0; i<nw; ++i){
      sums_weights(sums);
      grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
          fl.w[0]*fl.w[1]*fl.w[2]);
    }
  }
}
