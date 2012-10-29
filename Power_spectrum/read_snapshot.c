/*==========================================================================*/
/* Version: 1.00           date:21/10/08                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: routines to read the binary files from the Gadget2 simulations  */
/*                                                                          */
/* The original version can be find in the Gadget2 web page:                */
/* http://www.mpa-garching.mpg.de/gadget/                                   */
/*==========================================================================*/

#include "read_snapshot.h"

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

int load_header(char *fname, int files, io_header *header)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    pc;
  int NumPart, Ngas;   //number of DM and gas particles

  for(i=0, pc=0; i<files; i++){
    if(files>1) sprintf(buf,"%s.%d",fname,i);
    else sprintf(buf,"%s",fname);

    if(!(fd=fopen(buf,"r"))){
     // if(myrank==root) printf("can't open file `%s`\n",buf);
      exit(0);
    }

    //if(myrank==root) printf("reading the header of `%s' ...\n",buf);

    fread(&dummy, sizeof(dummy), 1, fd);
    fread(header, sizeof(*header), 1, fd);  // read and store the header
    fread(&dummy, sizeof(dummy), 1, fd);

    // counts particles
    if(files==1){
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++) NumPart+= header->npart[k];
      Ngas= header->npart[0];
    }
    else{
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++) NumPart+= header->npartTotal[k];
      Ngas= header->npartTotal[0];
    }
    fclose(fd);
  }
  return( NumPart );  //return the number of particles
}


/* this routine loads particle data from Gadget's default binary file format. (A snapshot may be distributed into multiple files. */
void load_snapshot(char *fname, int files, particle_data *P)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new;
  io_header header1;   //temporary header

  for(i=0, pc=0; i<files; i++, pc=pc_new){
    if(files>1) sprintf(buf,"%s.%d",fname,i);
    else sprintf(buf,"%s",fname);

    fd=fopen(buf,"r");

    //if(myrank==root) printf("reading `%s' ...\n",buf); fflush(stdout);

    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&header1, sizeof(header1), 1, fd);  // read and store the header
    fread(&dummy, sizeof(dummy), 1, fd);

    for(k=0, ntot_withmasses=0; k<5; k++) if(fabs(header1.mass[k])<1e-6) ntot_withmasses+= header1.npart[k];

    SKIP;
    for(k=0,pc_new=pc;k<6;k++){
      for(n=0;n<header1.npart[k];n++){
        fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);   //read the position of all the particles
        pc_new++;
      }
    }
    SKIP;

    SKIP;
    for(k=0,pc_new=pc;k<6;k++){
      for(n=0;n<header1.npart[k];n++){
        fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);   //read the velocity of all the particles
        pc_new++;
      }
    }
    SKIP;

    fclose(fd);
  }
}


