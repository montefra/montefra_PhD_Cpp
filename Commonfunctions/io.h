/*==========================================================================*/
/* version: 0.01, 19/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: collection of input-output functions I might need here and there*/
/*==========================================================================*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <gsl/gsl_spline.h>

#include "vector_array.h"

using namespace std;

/*==========================================================================*/
/* Error handling classes used for input and output files                   */
/*==========================================================================*/
/* custom class derived by exception to throw error reading file 'ifile'    */       
class inputexception : public exception
{
  string s;  //output string
  public:
    inputexception(string ifile) : s(ifile) {}  //constructor that takes a single file name
    ~inputexception () throw() {};   //destructor
    virtual const char* what() const throw() { return s.c_str(); }
};
/* custom class derived by exception to throw error regarding 'ofile'      */       
class outputexception : public exception
{
  string s;  //output string
  public:
    outputexception(string ofile) : s(ofile) {}  //constructor that takes a single file name
    ~outputexception () throw() {};   //destructor
    virtual const char* what() const throw() { return s.c_str(); }
};

/* check if file 'fname' exists. returns true if exist, false otherwise     */
bool fileexists(string fname);

/* read the first two columns of the input file name 'fname', store them    */
/* them into vector<double> x, y. The number of columns in the file is 'nc' */
/* returns '0' if no errors occured during the operation                    */
/* returns '1' if the input file does not exist                             */
int readfirst2columns(string fname, int nc, vector<double> &x, vector<double> &y);


/* write arrays 'x' and 'y' of size 'n' in file 'fname'                     */
/* returns '0' if no errors occured during the operation                    */
/* returns '1' if the output file already exist                             */
int writetwoarrays(string fname, double* x, double* y, size_t n);

/* read file 'fname', skipping lines starting with 'comment' or empty and   */
/* returns an vector of strings containing all the other lines              */
vector<string> readlines(string fname, string comment="#");

