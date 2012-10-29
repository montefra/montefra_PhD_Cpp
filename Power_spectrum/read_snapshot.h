/*==========================================================================*/
/* Version: 1.00           date:21/10/08                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: structures and variables to store the content of the binary     */
/* files from the Gadget2 simulations                                       */
/*                                                                          */
/* The original version can be find in the Gadget2 web page:                */
/* http://www.mpa-garching.mpg.de/gadget/                                   */
/*==========================================================================*/

#ifndef READ_SNAPSHOT_H
#define READ_SNAPSHOT_H

#include <cmath>
#include <cstdio>
#include <cstdlib>

/*==========================================================================*/
/* structures containing the header and some of the particle data from the  */
/* snapshot of Gadget2 simulation                                           */
/*==========================================================================*/
struct io_header
{
  int      npart[6];
  double   mass[6];
  double   time, redshift;
  int      flag_sfr, flag_feedback, flag_cooling;
  int      npartTotal[6];
  int      num_files;
  double   BoxSize;
  double   Omega0, OmegaLambda, HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};

struct particle_data
{
  float  Pos[3];
  float  Vel[3];
};

/*==========================================================================*/
/* modified function to read the binary files and to allocate the memory    */
/* for the struct containing the data                                       */
/*==========================================================================*/

int load_header(char *fname, int files, io_header *header);
void load_snapshot(char *fname, int files, particle_data *P);
#endif
