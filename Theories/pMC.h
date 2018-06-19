/*==========================================================================*/
/* Version: xxxx           date:xxxxxxxx                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: base and derived classes for the computation of P_{MC}(k)       */
/*                                                                          */
/*==========================================================================*/

#include <string>
#include <gsl/gsl_spline.h>

/*==========================================================================
 * virtual base class. Cantains a virtual constructor, destructor and two
 * functions: one accept the value of k_star and the other returns 
 * P_{MC}(k) evaluated in k
 *==========================================================================*/
class pmc
{
  protected:
    double kstar;   //value of kstar

  public:
    /*==========================================================================
     * Virtual destructor
     *==========================================================================*/
    //virtual ~pmc()=0;

    /*==========================================================================
     * set the value of kstar
     * Parameters
     * ----------
     * kst: value to be saved
     *==========================================================================*/
    void set_kstar(double kst){ kstar = kst; }

    /*==========================================================================
     * compute or returns P_{MC} evaluated in k
     * Parameters
     * ----------
     * k: wavenumber where to evaluate P_{MC}
     * output
     * ------
     * Pmc: P_{MC}(k)
     *==========================================================================*/
     virtual double pMC(double k) =0;
};

class pmc_interp: public pmc
{
  private:
    gsl_spline *spl1l;   //gsl spline of the 1loop ps
    gsl_interp_accel *acc1l;   //gsl accellation of the 1loop ps

  public:
    /*==========================================================================
     * constructor. Read the file name and allocate the spline and accelerator
     * objects
     * Parameters
     * ----------
     * spl: gsl_spline object: copied to a local pointer
     *==========================================================================*/
    pmc_interp( gsl_spline *spl ) {
      spl1l = spl;    // copy the spline in a local object
      //allocate the accelerator for the P_MC
      acc1l = gsl_interp_accel_alloc();    
    }
    /*==========================================================================
     * destructor. Deallocate the spline and accelerator objects
     *==========================================================================*/
    ~pmc_interp(){
      gsl_interp_accel_free(acc1l);   //gsl accellation of the linear and the 1loop ps
    }
  
    /*==========================================================================
     * returns P_{MC} evaluated in k
     * Parameters
     * ----------
     * k: wavenumber where to evaluate P_{MC}
     * output
     * ------
     * Pmc: P_{MC}(k)
     *==========================================================================*/
    double pMC(double k){
      return( gsl_spline_eval(spl1l, k, acc1l) );
    }
};

class pmc_dump: public pmc
{
  private:
    gsl_spline *spl1l;   //gsl spline of the 1loop ps
    gsl_interp_accel *acc1l;   //gsl accellation of the 1loop ps
  public:
    /*==========================================================================
     * constructor. Read the file name and allocate the spline and accelerator
     * objects
     * Parameters
     * ----------
     * spl: gsl_spline object: copied to a local pointer
     *==========================================================================*/
    pmc_dump( gsl_spline *spl ) {
      spl1l = spl;    // copy the spline in a local object
      //allocate the accelerator for the P_MC
      acc1l = gsl_interp_accel_alloc();    
    }
    /*==========================================================================
     * destructor. Deallocate the spline and accelerator objects
     *==========================================================================*/
    ~pmc_dump(){
      gsl_interp_accel_free(acc1l);   //gsl accellation of the linear and the 1loop ps
    }

    /*==========================================================================
     * returns P_{MC}(k) * exp( k^2/kstar^2 )
     * Parameters
     * ----------
     * k: wavenumber where to evaluate P_{MC}
     * output
     * ------
     * Pmc: P_{MC}(k) * exp( k^2/kstar^2 )
     *==========================================================================*/
    double pMC(double k){
      return( exp( - k*k / (kstar*kstar) ) * gsl_spline_eval(spl1l, k, acc1l) );
    }
};
