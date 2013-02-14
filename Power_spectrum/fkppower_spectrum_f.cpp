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

/*==========================================================================
 * Do the checks when the cross power spectrum is required
 * if files.size > 4, files is resized to 4
 * Parameters
 * ----------
 * files: vector of strings with the input files
 * myrank: rank of the processor
 * root: root processor
 * com: MPI communicator
 *==========================================================================*/
void check_crosspk(std::vector<std::string> &files, int myrank, int root,
    MPI_Comm com){
  size_t ninfiles = files.size(); //number of files
  std::cerr << "Not implemented" << std::endl;
  MPI_Barrier(com);
  MPI_Finalize();
  exit(1);
  if(ninfiles<4){
    if(myrank == root)
      std::cerr << "Please give the names of 4 files" << std::endl;
    MPI_Barrier(com);
    MPI_Finalize();
    exit(1);
  }
  else if(ninfiles>4){
    if(myrank == root) 
      std::cerr << "Only the first four files considered" << std::endl;
    files.resize(4);
  }
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
 * ignorew: int (-1,1,2): if !=-1 set to 1 fl.w[1] or fl.w[2]
 * repeatw: bool: if true, assign object fl.w[1] times with fl.w[1]=1
 * output
 * ------
 * sums: array containing sum(w), sum(n*w^2), sum(w^2)
 *==========================================================================*/
double *readfiles::read_file(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
    std::vector<double> zrange, int ignorew, bool repeatw){
  double *sums = init_sums(); //sum(w), sum(n*w^2), sum(w^2)
  
  if(zrange[1]<=0 && ignorew==-1 && repeatw==false) 
    read_asitis(ifile, grid, sums);
  else if(zrange[1]>0 && ignorew==-1 && repeatw==false) 
    read_zrange(ifile, grid, sums, zrange);
  else if(zrange[1]<=0 && ignorew!=-1 && repeatw==false) 
    read_ignorew(ifile, grid, sums, ignorew);
  else if(zrange[1]<=0 && ignorew==-1 && repeatw==true) 
    read_repeatw(ifile, grid, sums);
  else if(zrange[1]>0 && ignorew!=-1 && repeatw==false) 
    read_zrange_ignorew(ifile, grid, sums, zrange, ignorew);
  else if(zrange[1]>0 && ignorew==-1 && repeatw==true) 
    read_zrange_repeatw(ifile, grid, sums, zrange);
  else if(zrange[1]<=0 && ignorew!=-1 && repeatw==true) 
    read_ignorew_repeatw(ifile, grid, sums, ignorew);
  else
    read_zrange_ignorew_repeatw(ifile, grid, sums, zrange, ignorew);

  return(sums);  //return the number of haloes assigned to the grid
}

/*==========================================================================
 * read in the file with the various options. 
 * For argument refer to the description of 'read_file'
 *==========================================================================*/
void readfiles::read_asitis(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
    double* sums){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_zrange(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
    double* sums, std::vector<double> zrange){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
    if(zrange[0]>fl.z || zrange[1]<fl.z) continue;  //read the next line
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_ignorew(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
    double* sums, int ignorew){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
    fl.w[ignorew] = 1.;
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_repeatw(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
    double* sums){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
    int nw = fl.w[1];
    fl.w[1] = 1.;
    for(int i=0; i<nw; ++i){
      sums_weights(sums);
      grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
          fl.w[0]*fl.w[1]*fl.w[2]);
    }
  }
}
void readfiles::read_zrange_ignorew(std::string ifile, ps_r2c_c2r_mpi_inplace
    &grid, double* sums, std::vector<double> zrange, int ignorew){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
    if(zrange[0]>fl.z || zrange[1]<fl.z) continue;  //read the next line
    fl.w[ignorew] = 1.;
    sums_weights(sums);
    grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
        fl.w[0]*fl.w[1]*fl.w[2]);
  }
}
void readfiles::read_zrange_repeatw(std::string ifile, ps_r2c_c2r_mpi_inplace
    &grid, double* sums, std::vector<double> zrange){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
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
void readfiles::read_ignorew_repeatw(std::string ifile, ps_r2c_c2r_mpi_inplace
    &grid, double* sums, int ignorew){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
    fl.w[ignorew] = 1.;
    int nw = fl.w[1];
    fl.w[1] = 1.;
    for(int i=0; i<nw; ++i){
      sums_weights(sums);
      grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
          fl.w[0]*fl.w[1]*fl.w[2]);
    }
  }
}
void readfiles::read_zrange_ignorew_repeatw(std::string ifile,
    ps_r2c_c2r_mpi_inplace &grid, double* sums, std::vector<double> zrange, int
    ignorew){
  std::ifstream in(ifile.c_str());  //open the input file 
  for(;;){  //loop through the lines of the file
    read_line(in);
    if(in.eof() == true) break;    //break if the end of the file reached
    if(zrange[0]>fl.z || zrange[1]<fl.z) continue;  //read the next line
    fl.w[ignorew] = 1.;
    int nw = fl.w[1];
    fl.w[1] = 1.;
    for(int i=0; i<nw; ++i){
      sums_weights(sums);
      grid.assign_particle(fl.pos[0], fl.pos[1], fl.pos[2],
          fl.w[0]*fl.w[1]*fl.w[2]);
    }
  }
}
