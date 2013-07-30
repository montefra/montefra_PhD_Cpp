/*==========================================================================
 * PRIVATE FUNCTIONS
 *==========================================================================*/

#include "parse_ini.h"
/*==========================================================================
 * Construnctor
 * Parameters
 * ----------
 * inifile: string
 *   name of the input file contining the input parameters
 *==========================================================================*/
ParseIni::ParseIni(std::string inifile): ininame(inifile){ 
  this->comment = "#";  //default comment
  this->inis = this->readlines(inifile); 
  this->inisize = this->inis.size();
}

/*==========================================================================
 * read file 'inifile', skipping empty lines,
 * and stripping comments
 * Parameters
 * ----------
 * fname: string
 *   ini file name
 * output
 * ------
 * lines: vector of strings
 *   inifile stripped all all the empty lines and comments
 *==========================================================================*/
std::vector<std::string> ParseIni::readlines(std::string fname){
  std::string str;   //string containing each line
  std::vector<std::string> lines;   //vector containing all the usefull lines

  std::ifstream in(fname.c_str(), std::ifstream::in);   //open the file
  if(!in.is_open()){
    std::cerr << "Could not open file " << fname << std::endl;
    exit(70);
  }
  while(getline(in, str)){   //loop untill no new lines available
    
    if( str.empty() == true ) continue;  //if the line is empty, go to the next line
    else if( str.compare(str.find_first_not_of(" \t"), 1, this->comment) == 0 )  //if the first non empty 
      continue;       //character is a commend, skip the line
    else{
      int theresComment = str.find(this->comment);  //check if there is a commen after something else
      if( theresComment != -1 ) //if a comment is found strip it from the string
	str = str.substr( 0, theresComment );   //extract whatever is before the comment
      lines.push_back(str);  //append the line to the vector
    }
  }
  in.close();   //close the file
  return lines;
}

/*==========================================================================
 * PUBLIC FUNCTIONS
 *==========================================================================*/

/*==========================================================================
 * standard error when retreaving a value from a inifile
 * parameters
 * ----------
 *  param: string
 *    name of the searched parameter
 *  error: int
 *    error code to give to exit
 *==========================================================================*/
void ParseIni::ini_error(std::string param, int error){
  std::cerr << param << " keyword not found or empty in the inifile ";
  std::cerr << get_fname() << std::endl;
  exit(error);
}

/*==========================================================================
 * Read a single parameter following a given string and '='
 * Parameters
 * ----------
 * parname: string
 *   parameter to search
 * output
 * ------
 * value: std::string, int, double, bool
 *   first value following the parameter name and '='
 *   if bool:
 *     if the parameter after '=' is T,t,True,true retruns true
 *     if the parameter after '=' is F,f,False,false retruns false
 * return
 * ------
 * error: int
 *   0: if no error occurs; 1: if the parameter do not exists, 
 *   -1: if the parameter has no argument
 *   if bool:
 *     -2 if the parameter is not among the choises of above
 *==========================================================================*/
int ParseIni::get_param(std::string parname, std::string *value){
  int error=1;   //error code: by default parname not found
  for(size_t i=0; i<this->inisize; ++i){
    if(this->inis[i].compare(0,parname.size(), parname) == 0){  //if a parameter found
      //get everything after the '='
      std::string temp = this->inis[i].substr(inis[i].find("=")+1);   //get the substring
      if( temp.empty() == true )  //if the line is empty
	error = -1;    //set error -1
      else{
	error = 0;    //set error 0
	std::stringstream iniparse(temp);   //initialise stringstream
	iniparse >> *value;  //copy only the first element
	iniparse.str("");  //clear iniparse
	iniparse.clear();  //clear iniparse
      }
      break;  //break
    }
  }
  return(error);
}
int ParseIni::get_param(std::string parname, int *value){
  std::string str;   //string that will contain the int
  int err = get_param(parname, &str); //read the parameter
  int err_conversion;  //erro in the conversion from string to int
  *value = common::to_number<int>(str, &err_conversion);
  if(err_conversion!=0) err = -2;
  return(err);
}
int ParseIni::get_param(std::string parname, size_t *value){
  std::string str;   //string that will contain the int
  int err = get_param(parname, &str); //read the parameter
  int err_conversion;  //erro in the conversion from string to int
  *value = common::to_number<size_t>(str, &err_conversion);
  if(err_conversion!=0) err = -2;
  return(err);
}
int ParseIni::get_param(std::string parname, double *value){
  std::string str;   //string that will contain the int
  int err = get_param(parname, &str); //read the parameter
  int err_conversion;  //erro in the conversion from string to int
  *value = common::to_number<double>(str, &err_conversion);
  if(err_conversion!=0) err = -2;
  return(err);
}
int ParseIni::get_param(std::string parname, bool *value){
  std::string strbool;   //string that will contain the boolean
  int err = get_param(parname, &strbool); //read the parameter
  if(err == 0){   //if the parameter has been succesfully read
    size_t strb_size = strbool.size();  //get the size of the parameter
    //set the parameter to lower case
    if(strbool.compare(0, strb_size, "t") == 0 ||  //check if the string is "t" or "true"
	strbool.compare(0, strb_size, "true") == 0 ) 
      *value = true;    //save true
    else if(strbool.compare(0, strb_size, "f") == 0 || //check if the string is "f" or "false"
	strbool.compare(0, strb_size, "false") == 0 ) 
      *value = false;   //save false
    else err = -2;   //if is not false or true, returns error -2
  }
  return(err);
}
/*==========================================================================
 * Read n parameters and return them in an vector
 * Parameters
 * ----------
 * parname: string
 *   parameter to search
 * n: int
 *   number of parameters to read
 * output
 * ------
 * values: std vector of strint, int, double
 *   first n values following the parameter name and '='
 * return
 * ------
 * error: int
 *   0: if no error occurs; 1: if the parameter do not exists, 
 *   -1: if the parameter has less than n arguments
 *   -2: the values could not be converted to the desired type
 *==========================================================================*/
int ParseIni::get_param_list(std::string parname, size_t n,
    std::vector<std::string> &values){
  int error=1;   //error code: by default parname not found
  for(size_t i=0; i<this->inisize; ++i){
    if(this->inis[i].compare(0,parname.size(), parname) == 0){  //if a parameter found
      //get everything after the '='
      std::string tempstring = this->inis[i].substr(inis[i].find("=")+1);   //get the substring
      if( tempstring.empty() == true )  //if the line is empty
	error = -1;    //set error -1
      else{
	std::stringstream iniparse(tempstring);   //initialise stringstream
	std::vector<std::string> tempvector;  //create a temporary vector
        std::string tempvar;  //create a temporary variable
	while( !iniparse.eof() ){   //copy stringstream into the temporary vector
	  iniparse >> tempvar;
	  tempvector.push_back(tempvar); 
	}
	if(tempvector.size() < n) error = -1;  //if there are less than n arguments error
	else{
	  error = 0;  //no error
	  values.assign( tempvector.begin(), tempvector.begin()+n );  //copy the first n values
	}
	iniparse.str("");  //clear iniparse
	iniparse.clear();  //clear iniparse
	tempvector.clear();  //clear the temporary vector
      }
      break;  //break
    }
  }
  return(error);
}
int ParseIni::get_param_list(std::string parname, size_t n, std::vector<int> &values){
  std::vector<std::string> vstr;   //string that will contain the int
  int err = get_param_list(parname, n, vstr); //read the parameter
  for(size_t i=0; i<vstr.size(); ++i){
    int err_conversion;  //erro in the conversion from string to int
    values.push_back(common::to_number<int>(vstr[i], &err_conversion));
    if(err_conversion!=0){
      err = -2;
      break;
    }
  }
  return(err);
}
int ParseIni::get_param_list(std::string parname, size_t n, std::vector<size_t> &values){
  std::vector<std::string> vstr;   //string that will contain the int
  int err = get_param_list(parname, n, vstr); //read the parameter
  for(size_t i=0; i<vstr.size(); ++i){
    int err_conversion;  //erro in the conversion from string to int
    values.push_back(common::to_number<size_t>(vstr[i], &err_conversion));
    if(err_conversion!=0){
      err = -2;
      break;
    }
  }
  return(err);
}
int ParseIni::get_param_list(std::string parname, size_t n, std::vector<double> &values){
  std::vector<std::string> vstr;   //string that will contain the int
  int err = get_param_list(parname, n, vstr); //read the parameter
  for(size_t i=0; i<vstr.size(); ++i){
    int err_conversion;  //erro in the conversion from string to int
    values.push_back(common::to_number<double>(vstr[i], &err_conversion));
    if(err_conversion!=0){
      err = -2;
      break;
    }
  }
  return(err);
}
/*==========================================================================
 * Read MCMC paramter starting guess, boundaries and sigma
 * Parameters
 * ----------
 * parname: string
 *   name of the parameter in the inifile
 * output
 * ------
 * guess: double
 *   guess of the initial position
 * lower: double
 *   lower limit of the flat prior
 * upper: double
 *   upper limit of the flat prior
 * sigma: double
 *   standard deviation of the gaussian used for the random walk
 * return
 * ------
 * error: int
 *   0: if no error occurs; 1: if the parameter do not exists, 
 *   -1: if less than 4 arguments are found
 *==========================================================================*/
int ParseIni::get_mcmc_params(std::string parname, double *guess, double* lower,
    double *upper, double *sigma){
  //parameter to search in inis
  std::string par = "param["+parname+"]";
  std::vector<double> paras; //vector containing the parameters 
  int error = this->get_param_list(par, 4, paras);
  if(error == 0){  //if there are no problems copy
    *guess = paras[0];   //copy the initial guess
    *lower = paras[1];    //copy the lower bound
    *upper = paras[2];    //copy the upper bound
    *sigma = paras[3];    //copy the step size
  }
  paras.clear();
  return( error );
}

