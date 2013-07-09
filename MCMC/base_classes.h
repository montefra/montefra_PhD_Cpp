/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: Markov Chain Monte Carlo, abstract base class for reading the
 * data
 * This provide the basic interface that must be implemented to read
 * the input files
 *==========================================================================*/

#include <vector>

#include "parse_ini.h"

/*==========================================================================
 * Base class for reading the input data. 
 * Only the constructor is required. 
 * This class simply provides a template used in the contructor of 
 * 'LikelihoodBase'
 *==========================================================================*/
class ReadDataBase{
  public:
    /*==========================================================================
     * Constructor of the class
     * Parameters
     * ----------
     *  iniparser: instance of ParseIni object
     *==========================================================================*/
    virtual ReadDataBase(ParseIni iniparser) = 0;

    /*==========================================================================
     * Methods that returns the data. 
     * Examples:
     * 
     * template <class T>
     *   virtual std::vector<T> get_data() =0;
     * template <class T>
     *   virtual std::vector<T> get_x() =0;
     * template class<T>
     *   virtual T get_data(x)=0;
     *==========================================================================*/

};

/*==========================================================================
 * Base class that implements the model
 * This class simply provides a template used in the contructor of 
 * 'LikelihoodBase'
 *==========================================================================*/
class ModelBase{
  public:
    /*==========================================================================
     * Methods that returns the model given the parameters and x
     * Examples:
     * 
     * template <class T>
     *   virtual std::vector<T> get_model(std::vector<T> x, ) =0;
     * template <class T>
     *   virtual T get_model(T x, ) =0;
     *==========================================================================*/

};
