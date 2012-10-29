/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/* Definition of the functions                                              */
/*==========================================================================*/

#include "boxpower_spectrum.h"

/*=======================================================================*/
/* Code specific functions definitions                                   */
/*=======================================================================*/
std::string get_description()
/* description of the code   */
{
  std::string description = "This code computes the angular averaged power spectrum ";
  description += "from catalogues with a cubic geometry.";
  return(description);
}
std::string get_version()
/*version of the code */
{ return("1.00"); }

/*==========================================================================
 * read a line from the file, store position and velocity and return the mass
 * Parameters
 * ----------
 * inif: ifstream object with the file to read
 * output:
 * mass: mass of the particle
 * pos[3], vel[3]: position and velocity of the particle
 *==========================================================================*/
double read_line(std::ifstream &inif, double *pos, double *vel){
  double mass;  //mass
  inif >> pos[0] >> pos[1] >> pos[2] >> vel[0] >> vel[1] >> vel[2] >> mass;   //read a line
  return(mass);
}

/*=======================================================================*/
/* read halo catalogues from ascii file                                  */
/*=======================================================================*/
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * no mass limit, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double cell_size){
  std::ifstream in(ifile.c_str());  //open the input file 

  double *pos=new double[3], *vel=new double[3];  //array containing the of position and velocity
  double m;  //mass of the halo
  ptrdiff_t n_haloes=0;  //number of haloes put in the grid

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached
    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    grid.assign_particle(pos[0], pos[1], pos[2], 1.);
    n_haloes +=1 ; //increment the number of objects
  }
  delete [] pos;   //clear velocity and position
  delete [] vel;
  return(n_haloes);  //return the number of haloes assigned to the grid
}
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * lower mass limit, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin: lower limit for the mass
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double cell_size){
  std::ifstream in(ifile.c_str());  //open the input file 

  double pos[3], vel[3];  //array containing the three components of position and velocity
  double m;  //mass of the halo
  ptrdiff_t n_haloes=0;  //number of haloes put in the grid

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached
    if(m<mmin) continue;   //if the mass is smaller than the desired one, go to the next line
    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    grid.assign_particle(pos[0], pos[1], pos[2], 1.);
    n_haloes +=1 ; //increment the number of objects
  }
  return(n_haloes);  //return the number of haloes assigned to the grid
}
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * mass bin, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin, mmax: lower and upper limit for the mass
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double mmax, double cell_size){
  std::ifstream in(ifile.c_str());  //open the input file 

  double pos[3], vel[3];  //array containing the three components of position and velocity
  double m;  //mass of the halo
  ptrdiff_t n_haloes=0;  //number of haloes put in the grid

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached
    if(m<mmin || m>mmax) continue;   //if the mass is not in the desired bin, go to the next line
    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    grid.assign_particle(pos[0], pos[1], pos[2], 1.);
    n_haloes +=1 ; //increment the number of objects
  }
  return(n_haloes);  //return the number of haloes assigned to the grid
}

/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * no mass limit, redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double cell_size,
    int zdist, ptrdiff_t ncells){
  std::ifstream in(ifile.c_str());  //open the input file 

  double pos[3], vel[3];  //array containing the three components of position and velocity
  double m;  //mass of the halo
  ptrdiff_t n_haloes=0;  //number of haloes put in the grid

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached

    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    pos[zdist] += vel[zdist] / cell_size;   // z distortions in z distortions
    if(pos[zdist] < 0.) pos[zdist] += ncells;   //periodic boundary conditions
    else if(pos[zdist] >= ncells) pos[zdist] -= ncells;

    grid.assign_particle(pos[0], pos[1], pos[2], 1.);
    n_haloes +=1 ; //increment the number of objects
  }
  return(n_haloes);  //return the number of haloes assigned to the grid
}
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * lower mass limit, redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin: lower limit for the mass
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double cell_size, int zdist, ptrdiff_t ncells){
  std::ifstream in(ifile.c_str());  //open the input file 

  double pos[3], vel[3];  //array containing the three components of position and velocity
  double m;  //mass of the halo
  ptrdiff_t n_haloes=0;  //number of haloes put in the grid

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached
    if(m<mmin) continue;   //if the mass is smaller than the desired one, go to the next line

    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    pos[zdist] += vel[zdist] / cell_size;   // z distortions in z distortions
    if(pos[zdist] < 0.) pos[zdist] += ncells;   //periodic boundary conditions
    else if(pos[zdist] >= ncells) pos[zdist] -= ncells;

    grid.assign_particle(pos[0], pos[1], pos[2], 1.);
    n_haloes +=1 ; //increment the number of objects
  }
  return(n_haloes);  //return the number of haloes assigned to the grid
}
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * mass bin, redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin, mmax: lower and upper limit for the mass
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double mmax, double cell_size, int zdist, ptrdiff_t ncells){
  std::ifstream in(ifile.c_str());  //open the input file 

  double pos[3], vel[3];  //array containing the three components of position and velocity
  double m;  //mass of the halo
  ptrdiff_t n_haloes=0;  //number of haloes put in the grid

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached
    if(m<mmin || m>mmax) continue;   //if the mass is not in the desired bin, go to the next line

    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    pos[zdist] += vel[zdist] / cell_size;   // z distortions in z distortions
    if(pos[zdist] < 0.) pos[zdist] += ncells;   //periodic boundary conditions
    else if(pos[zdist] >= ncells) pos[zdist] -= ncells;

    grid.assign_particle(pos[0], pos[1], pos[2], 1.);
    n_haloes +=1 ; //increment the number of objects
  }
  return(n_haloes);  //return the number of haloes assigned to the grid
}

/*==========================================================================
 * Read the input ascii file containing haloes and fill the grids
 * cross power, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid1, grid2: 'ps_r2c_c2r_mpi_inplace' objects
 * mass: 4 element vector containing the mass limits
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t *read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid1, ps_r2c_c2r_mpi_inplace &grid2, 
    std::vector<double> mass, double cell_size){
  std::ifstream in(ifile.c_str());  //open the input file 

  double pos[3], vel[3];  //array containing the three components of position and velocity
  double m;  //mass of the halo
  ptrdiff_t *n_haloes = new ptrdiff_t[3];  //number of haloes put in the grid
  n_haloes[0] = n_haloes[1] = n_haloes[2] = 0; //initialise the array 
  bool both = false;   //find how many objects are in both grids

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached
    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    if( m>=mass[0] && m<mass[1] ){
      grid1.assign_particle(pos[0], pos[1], pos[2], 1.);
      n_haloes[0] +=1 ; //increment the number of objects
      both = true;   //set to true if the particle is in the first sample
    }
    if( m>=mass[2] && m<mass[3] ){
      grid2.assign_particle(pos[0], pos[1], pos[2], 1.);
      n_haloes[1] +=1 ; //increment the number of objects
      if( both == true) //if the particle is also in the second sample
	n_haloes[2] +=1;  //add 1 to the third element of the n_haloes
    }
    both = false;  //set back both to false for the next particle
  }
  return(n_haloes);  //return the number of haloes assigned to the grid
}
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grids
 * cross power, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid1, grid2: 'ps_r2c_c2r_mpi_inplace' objects
 * mass: 4 element vector containing the mass limits
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t *read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid1, ps_r2c_c2r_mpi_inplace &grid2, 
    std::vector<double> mass, double cell_size, int zdist, ptrdiff_t ncells){
  std::ifstream in(ifile.c_str());  //open the input file 

  double pos[3], vel[3];  //array containing the three components of position and velocity
  double m;  //mass of the halo
  ptrdiff_t *n_haloes = new ptrdiff_t[3];  //number of haloes put in the grid
  n_haloes[0] = n_haloes[1] = n_haloes[2] = 0; //initialise the array 
  bool both = false;   //find how many objects are in both grids

  for(;;){  //loop through the lines of the file
    m = read_line(in, pos, vel);
    if (in.eof() == true) break;    //breack if the end of the file reached

    for(int i=0; i<3; ++i) pos[i] /= cell_size;   //convert the position in grid units
    pos[zdist] += vel[zdist] / cell_size;   // z distortions in z distortions
    if(pos[zdist] < 0.) pos[zdist] += ncells;   //periodic boundary conditions
    else if(pos[zdist] >= ncells) pos[zdist] -= ncells;

    if( m>=mass[0] && m<mass[1] ){
      grid1.assign_particle(pos[0], pos[1], pos[2], 1.);
      n_haloes[0] +=1 ; //increment the number of objects
      both = true;   //set to true if the particle is in the first sample
    }
    if( m>=mass[2] && m<mass[3] ){
      grid2.assign_particle(pos[0], pos[1], pos[2], 1.);
      n_haloes[1] +=1 ; //increment the number of objects
      if( both == true) //if the particle is also in the second sample
	n_haloes[2] +=1;  //add 1 to the third element of the n_haloes
    }
    both = false;  //set back both to false for the next particle
  }
  return(n_haloes);  //return the number of haloes assigned to the grid
}

/*=======================================================================*/
/* read dark matter catalogues from binary file                          */
/*=======================================================================*/
/*==========================================================================
 * Read the input binary file containing dark matter and fill the grid
 * no redshift space distortions. Gadget 2 binary file
 * Parameters
 * ----------
 * ifile: input file root
 * n_files: number of files
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * output
 * ------
 * n_dm: number of dark matter particles assigned to the grid
 *==========================================================================*/
ptrdiff_t read_dm(std::string ifile, int n_files, ps_r2c_c2r_mpi_inplace &grid, 
    double cell_size){
  std::ifstream in(ifile.c_str());  //open the input file 

  char* fname;   //file name to pass to a pure c routine
  particle_data *P;   //pointer to a structure that contains DM position and velocity
  io_header *header = new io_header;   //allocate the header
  double *pos=new double[3];  //array containing the of position and velocity
  ptrdiff_t n_dm=0;  //number of dark matter particles put in the grid

  for(int i=0; i<n_files; ++i){  //loop through the input files
    std::string temp = ifile + "." + to_string(i);   //create the file name
    fname = new char[temp.size() + 1];   //allocate the char
    strcpy(fname, temp.c_str());   //copy the string into the char
    /* load the header and allocate the memory needed */
    int NumPart = load_header(fname, 1, header);
    P= new particle_data[NumPart];   //allocate the structure

    load_snapshot(fname, 1, P);   //load each snapshot, one each time
    delete [] fname;   //deallocate the char

    for(int j=0; j<NumPart; ++j){   //loop over the particles in each snapshot
      for(int l=0; l<3; ++l) 
        pos[l] = P[j].Pos[l]/cell_size;  //copy the positions

      grid.assign_particle(pos[0], pos[1], pos[2], 1.);
      n_dm += 1;  //increment the number of particles
    }
    delete[] P;   // free the memory before the next file or exiting the loop
    NumPart=0;
  }
  delete header;   //deallocate the header
  delete [] pos;  //deallocate the position array
  return(n_dm);  //return the number of haloes assigned to the grid
}
/*==========================================================================
 * Read the input binary file containing dark matter and fill the grid
 * redshift space distortions. Gadget 2 binary file
 * Parameters
 * ----------
 * ifile: input file root
 * n_files: number of files
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_dm: number of dark matter particles assigned to the grid
 *==========================================================================*/
ptrdiff_t read_dm(std::string ifile, int n_files, ps_r2c_c2r_mpi_inplace &grid, 
    double cell_size, int zdist, ptrdiff_t ncells){
  std::ifstream in(ifile.c_str());  //open the input file 

  char* fname;   //file name to pass to a pure c routine
  particle_data *P;   //pointer to a structure that contains DM position and velocity
  io_header *header = new io_header;   //allocate the header
  double *pos=new double[3];  //array containing the of position
  ptrdiff_t n_dm=0;  //number of dark matter particles put in the grid

  for(int i=0; i<n_files; ++i){  //loop through the input files
    std::string temp = ifile + "." + to_string(i);   //create the file name
    fname = new char[temp.size() + 1];   //allocate the char
    strcpy(fname, temp.c_str());   //copy the string into the char
    /* load the header and allocate the memory needed */
    int NumPart = load_header(fname, 1, header);
    P = new particle_data[NumPart];   //allocate the structure

    load_snapshot(fname, 1, P);   //load each snapshot, one each time
    delete [] fname;   //deallocate the char

    double vel_Mpch = 100. * sqrt(header->time) * sqrt( header->Omega0 * pow(header->time, -3.) + 
        header->OmegaLambda );  // correction to have the velocity in Mpc/h

    for(int j=0; j<NumPart; ++j){   //loop over the particles in each snapshot
      for(int l=0; l<3; ++l) 
        pos[l] = P[j].Pos[l] / cell_size;   //copy the positions in grid units

      pos[zdist] += P[j].Vel[zdist] / (vel_Mpch * cell_size);   // z distortions in cells units
      if( pos[zdist] < 0. ) pos[zdist] += ncells;   //periodic boundary conditions
      else if( pos[zdist] >= ncells ) pos[zdist] -= ncells;

      grid.assign_particle(pos[0], pos[1], pos[2], 1.);
      n_dm += 1;  //increment the number of particles
    }
    delete[] P;   // free the memory before the next file or exiting the loop
    NumPart=0;
  }

  delete header;   //deallocate the header
  delete [] pos;  //deallocate the position array
  return(n_dm);  //return the number of haloes assigned to the grid
}

