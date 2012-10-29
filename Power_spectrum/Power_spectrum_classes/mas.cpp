/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* MAS for the different cases. Each mas is declared in a class             */
/* function definitions                                                     */
/*==========================================================================*/

#include "mas.h"

/*==========================================================================
 * Compute the weights and the cells to which this weights are to be 
 * applied for the CIC MAS
 * Parameters
 * ----------
 * c: position in cells units
 * icells: number of cells in the dimension considered
 * n: number of cells to used for the MAF (for CIC must be 2)
 * output
 * ------
 * cells: 'n' dimension array containing the cells over which the particle 
 *   is spread
 * w: 'n' dimensional array containing the weights per cell
 *==========================================================================*/
void MAS_common::W_c_CIC(double c, ptrdiff_t ncells, int n, double *w, ptrdiff_t *cells){
  double fc = floor(c);  //take the floor of c
  c -= fc;   //subtract it from c

  for( ptrdiff_t i=0; i<n; ++i) w[i] = 1. - fabs(c-i);  //weights
  cells[0] = fc;  //get the cells to which the weights are applied
  cells[1] = (fc+1 > ncells) ? 0 : fc+1;
}

/*==========================================================================
 * Compute the weights and the cells to which this weights are to be 
 * applied for the TSC MAS
 * Parameters
 * ----------
 * c: position in cells units
 * icells: number of cells in the dimension considered
 * n: number of cells to used for the MAF (for TSC should be 3)
 * output
 * ------
 * cells: 'n' dimension array containing the cells over which the particle 
 *   is spread
 * w: 'n' dimensional array containing the weights per cell
 *==========================================================================*/
void MAS_common::W_c_TSC(double c, ptrdiff_t icells, int n, double *w, ptrdiff_t *cells)
{
  double fc = floor(c);   //take the floor of c
  c = c - fc;   //subtract it from c
  ptrdiff_t disp = (fabs(c) < 0.5) ? 1 : 0;   //and check if the particle is in the first or the second half of the cell
  c = c + disp;  //add the displacement
  ptrdiff_t ncell = fc - disp;  //cell number: just to make less calculations

  for(ptrdiff_t i=0; i<n; ++i){
    w[i] = fabs(c - double(i));
    w[i] = (w[i] < 0.5) ? (0.75 - pow(w[i], 2.)) : (pow(1.5 - w[i], 2.)/2.);

    if(ncell + i < 0) cells[i] = ncell + i + icells;
    else if(ncell + i >= icells) cells[i] = ncell + i - icells;
    else cells[i] = ncell + i;
  }
}


/*==========================================================================
 * CIC MAS for the real grid
 * Parameters
 * ----------
 * rgrid: double FFTW grid
 * x, y, z: particle position in grid units
 * w: weight of the particle
 *==========================================================================*/
void CIC_r2c_c2r_mpi::MAS(double *rgrid, double x, double y, double z, double w)  
{
  int n = 2; //dimention of the MAS
  double wx[n], wy[n], wz[n];  // TSC weights for each dimension
  ptrdiff_t cellx[n], celly[n], cellz[n];  //cells to which the weights are assigned

  W_c_CIC(x, xcells, n, wx, cellx);   //get weights and cells where to apply the weight
  W_c_CIC(y, ycells, n, wy, celly);
  W_c_CIC(z, zcells, n, wz, cellz);

  ptrdiff_t i, j, l;       //loop integral
  for(i=0; i<2; ++i){
    if(cellx[i] >= local_x_start && cellx[i] < local_x_finish)
      for(j=0; j<2; ++j) for(l=0; l<2; ++l)
	rgrid[cellz[i] + 2*zcellso2_1 * (celly[i] + ycells * (cellx[i]-local_x_start))] += 
	  w * wx[i] * wy[j] * wz[l];
  }
}
/*==========================================================================
 * CIC MAS for the complex grid
 * Parameters
 * ----------
 * cgrid: complex FFTW grid
 * x, y, z: particle position in grid units
 * ind: 0 or 1 to assigne the particle to the real or immaginary part of the grid
 * w: weight of the particle
 *==========================================================================*/
void CIC_r2c_c2r_mpi::MAS(fftw_complex *cgrid, double x, double y, double z, int ind, double w)
{
  int n = 2; //dimention of the MAS
  double wx[n], wy[n], wz[n];  // TSC weights for each dimension
  ptrdiff_t cellx[n], celly[n], cellz[n];  //cells to which the weights are assigned

  W_c_CIC(x, xcells, n, wx, cellx);   //get weights and cells where to apply the weight
  W_c_CIC(y, ycells, n, wy, celly);
  W_c_CIC(z, zcells, n, wz, cellz);

  ptrdiff_t i, j, l;       //loop integral
  for(i=0; i<2; ++i){
    if( cellx[i] >= local_x_start && cellx[i] < local_x_finish )
      for(j=0; j<2; ++j) for(l=0; l<2; ++l){
	if( cellz[i]<zcellso2_1 )
	  cgrid[cellz[i] + zcellso2_1 * (celly[i] + ycells * (cellx[i]-local_x_start))][ind] += 
	    w * wx[i] * wy[j] * wz[l];
      }
  }
}


/*==========================================================================
 * TSC MAS for the real grid
 * Parameters
 * ----------
 * rgrid: double FFTW grid
 * x, y, z: particle position in grid units
 * w: weight of the particle
 *==========================================================================*/
void TSC_r2c_c2r_mpi::MAS(double *rgrid, double x, double y, double z, double w)  
{
  int n = 3;  // number of cell to spread the signal
  double wx[n], wy[n], wz[n];  // TSC weights for each dimension
  ptrdiff_t cellx[n], celly[n], cellz[n];  //cells to which the weights are assigned

  W_c_TSC(x, xcells, n, wx, cellx);   //get weights and cells where to apply the weight
  W_c_TSC(y, ycells, n, wy, celly);
  W_c_TSC(z, zcells, n, wz, cellz);

  int i, j, l;       //loop integral
  for(i=0; i<n; ++i){
    if(cellx[i] >= local_x_start && cellx[i] < local_x_finish){
      for(j=0; j<n; ++j) for(l=0; l<n; ++l)
	rgrid[ cellz[l] + 2*zcellso2_1 * (celly[j] + ycells * (cellx[i]-local_x_start)) ] += 
	  w * wx[i] * wy[j] * wz[l];
    }
  }
}
/*==========================================================================
 * TSC MAS for the complex grid
 * Parameters
 * ----------
 * cgrid: complex FFTW grid
 * x, y, z: particle position in grid units
 * ind: 0 or 1 to assigne the particle to the real or immaginary part of the grid
 * w: weight of the particle
 *==========================================================================*/
void TSC_r2c_c2r_mpi::MAS(fftw_complex *cgrid, double x, double y, double z, int ind, double w)
{
  int n = 3;  // number of cell to spread the signal
  double wx[n], wy[n], wz[n];  // TSC weights for each dimension
  ptrdiff_t cellx[n], celly[n], cellz[n];  //cells to which the weights are assigned

  W_c_TSC(x, xcells, n, wx, cellx);   //get weights and cells where to apply the weight
  W_c_TSC(y, ycells, n, wy, celly);
  W_c_TSC(z, zcells, n, wz, cellz);

  int i, j, l;       //loop integral
  for(i=0; i<n; ++i){
    if(cellx[i] >= local_x_start && cellx[i] < local_x_finish){
      for(j=0; j<n; ++j) for(l=0; l<n; ++l)
	if( cellz[l]<zcellso2_1 )
	  cgrid[ cellz[l] + zcellso2_1*(celly[j] + ycells * (cellx[i]-local_x_start)) ][ind] += w * wx[i] * wy[j] * wz[l];
    }
  }
}

