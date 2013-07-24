/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: read the pk data to use in the mcmc chain
 *==========================================================================*/


/*==========================================================================
 * Constructor of the object
 * Parameter
 * ---------
 *  pk_file: string with the file name containing the name of the file to
 *  read in and the bins to consider
 *==========================================================================*/
Read_pk_files::Read_pk_files(std::string pk_dataset){
  
  //import the dataset file into a inifile object
  dataset = new ParseIni(pk_dataset);
  //get the number of bins
  get_bin_numbers(dataset);

  
}


/*==========================================================================
 * read the number of bins to read from the files and to use for the
 * likelihood
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Read_pk_files::get_bin_numbers(ParseIni ini){
    if(ini.get_param("kitot", &n_bins.kitot) != 0) exit(10);
    if(ini.get_param("kimin", &n_bins.kimin) != 0) exit(11);
    if(ini.get_param("kimax", &n_bins.kimax) != 0) exit(12);
    if(ini.get_param("kjtot", &n_bins.kjtot) != 0) exit(13);
    if(ini.get_param("kjmin", &n_bins.kjmin) != 0) exit(14);
    if(ini.get_param("kjmax", &n_bins.kjmax) != 0) exit(15);

    n_bins.kiuse = n_bins.kimax-n_bins.kimin;
    n_bins.kjuse = n_bins.kjmax-n_bins.kjmin;
}

/*==========================================================================
 * read measured power
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Read_pk_files::read_pk(ParseIni ini){
  std::string pkfile;
  if(ini.get_param("pk_file", &pk_file);

}

