/*==========================================================================*/
/* version: 0.01, 19/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: collection of input-output functions I might need here and there*/
/*==========================================================================*/

#include "io.h"

/*==========================================================================*/
/* Check if file 'fname' exists. Returns true if exists, false otherwise     */
/*==========================================================================*/
bool fileexists(string fname)
{
  ifstream check(fname.c_str(), ifstream::in);
  if (check.is_open()){
    check.close();
    return(true);
  }
  else return(false);
}


/*==========================================================================*/
/* read the first two columns of the input file name 'fname', store them    */
/* them into vector<double> x, y. The number of columns in the file is 'nc' */
/* returns '0' if no errors occured during the operation                    */
/* returns '1' if the input file does not exist                             */
/*==========================================================================*/
int readfirst2columns(string fname, int nc, vector<double> &x, vector<double> &y)
{
  int i;   //loop integer
  double xi, yi, dummy;   //first two elements of each line and dummy argument for all the other lines
  if( fileexists(fname) == false ) return(1);
  
  ifstream in(fname.c_str(), ifstream::in);   //open the file object

  for(;;){
    in >> xi >> yi;   //read the first two columns
    for(i=2; i<nc; ++i) in >> dummy;   //skip all the columns from 2 to nc
    if(in.eof() == true) break;   //get out of the loop if the file is finished
    else{
      x.push_back(xi);
      y.push_back(yi);
    }
  }
  in.close();   //close the file object
  return(0);
}

/*==========================================================================*/
/* write arrays 'x' and 'y' of size 'n' in file 'fname'                     */
/* returns '0' if no errors occured during the operation                    */
/* returns '1' if the output file already exist                             */
/*==========================================================================*/
int writetwoarrays(string fname, double* x, double* y, size_t n)
{
  size_t i;  //loop integer
  if( fileexists(fname) == true ) return(1);  //returns 1 if the output file name already exists

  ofstream out(fname.c_str(), ofstream::out);     //open the output file objext
  out.setf(ios_base::scientific);
  out.precision(6);
  out.width(9);
  for(i=0; i<n; ++i) out << x[i] << "\t" << y[i] << endl;
  out.close();
  return(0);
}

/*==========================================================================*/
/* read file 'inifile', skipping lines starting with 'comment' or empty and */
/* returns an vector of strings containing all the other lines              */
/*==========================================================================*/
vector<string> readlines(string fname, string comment)
{
  string str;
  vector<string> lines;
  fstream file;

  ifstream in(fname.c_str());
  while(!in.eof()){
    getline(in, str);
    if(str.compare(0,1,comment) != 0 && str.empty()==false) lines.push_back(str);
  }
  return lines;
}

