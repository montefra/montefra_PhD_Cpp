/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: container for functions using gsl
 *==========================================================================*/

#include "gsl_funcs.h"

/*==========================================================================
 * Skip commented lines in file
 * Parameters
 * ----------
 *  file: c file pointer
 *    file pointer (c style)
 *  comment: char
 *    comment character
 *==========================================================================*/
void gslf::skip_comment_lines(FILE *file, char comment){
  while(true){
    //get the file position before reading the next char
    fpos_t pos;  
    fgetpos(file, &pos); 
    int c = fgetc(file); //read a char
    if(c == '#'){ //if it is a comment
      do{  //go ahead until the end of the line is hit
        c = fgetc(file);
      }while(c != '\n');
    }
    else{ //if the line does not start with a comment
      //set back the position at the first position in the line
      fsetpos(file, &pos); 
      break;
    }
  }
}

/*==========================================================================
 * read a vector of size dim from file 'file_name'
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the vector
 *  dim: int
 *    dimention of the input vector
 * output
 * ------
 *  vector: gsl vector
 *==========================================================================*/
gsl_vector *gslf::read_gsl_vector(std::string file_name, int dim){
  // create a gsl vector and read in the whole file
  gsl_vector *vector = gsl_vector_alloc(dim);
  FILE *f = fopen(file_name.c_str(), "r");   //file unit
  if(f == NULL){
    std::cerr << "Error opening file " << file_name << std::endl;
    exit(80);
  }
  skip_comment_lines(f, '#');

  if(gsl_vector_fscanf(f, vector) !=0){
    std::cerr << "The vector in file  " <<file_name << " can't be read." << std::endl;
    exit(81);
  }
  fclose(f);

  return vector;
}
/*==========================================================================
 * read a vector of size dim from file 'file_name', cut 'cut_dim' using imin as
 * the minimum index
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the vector
 *  dim: int
 *    dimention of the input vector
 *  cut_dim: int
 *    dimension of the cut vector
 *  imin: int
 *    minimum index in the cut vector
 * output
 * ------
 *  vector: gsl vector
 *==========================================================================*/
gsl_vector *gslf::read_cut_gsl_vector(std::string file_name, int dim, int cut_dim,
    int imin){
  //read the full vector
  gsl_vector *temp_input = read_gsl_vector(file_name, dim);
  //allocate the output vector
  gsl_vector *vector = gsl_vector_alloc(cut_dim);
  //take a view of the subvector
  gsl_vector_view vview = gsl_vector_subvector(temp_input, imin, cut_dim);
  //copy the interesting part of the vector
  gsl_vector_memcpy(vector, &(vview.vector));
  gsl_vector_free(temp_input);
  return vector;
}

/*==========================================================================
 * read a matrix of size xdim*ydim from file 'file_name'
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the matrix
 *  xdim, ydim: int
 *    dimentions of the input matrix
 * output
 * ------
 *  matrix: gsl matrix
 *    output gsl matrix
 *==========================================================================*/
gsl_matrix *gslf::read_gsl_matrix(std::string file_name, int xdim, int
    ydim){
  // create a gsl matrix and read in the whole file
  gsl_matrix *matrix = gsl_matrix_alloc(xdim, ydim); 
  FILE *f = fopen(file_name.c_str(), "r");   //file unit
  if(f == NULL){
    std::cerr << "Error opening file " << file_name << std::endl;
    exit(82);
  }
  skip_comment_lines(f, '#');

  if(gsl_matrix_fscanf(f, matrix) !=0){
    std::cerr << "The matrix in file  " <<file_name << " can't be read." << std::endl;
    exit(83);
  }
  fclose(f);

  return(matrix);
}


/*==========================================================================
 * read a matrix of size xdim*ydim from file 'file_name', cut to 
 * 'cut_xdim*cut_ydim' using xmin and ymin as the mimimum x and y indeces
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the matrix
 *  xdim, ydim: int
 *    dimentions of the input matrix
 *  cut_xdim, cut_ydim: int
 *    dimension of the cut matrix
 *  xmin, ymin: int
 *    minimum x and y index in the cut matrix
 * output
 * ------
 *  matrix: gsl matrix
 *    output cut gsl matrix
 *==========================================================================*/
gsl_matrix *gslf::read_cut_gsl_matrix(std::string file_name, int xdim, int ydim, int
    cut_xdim, int cut_ydim, int xmin, int ymin){
  //read the full matrix
  gsl_matrix *temp_input = read_gsl_matrix(file_name, xdim, ydim);
  //allocate the output matrix
  gsl_matrix *matrix = gsl_matrix_alloc(cut_xdim, cut_ydim);
  //take a view of the submatrix
  gsl_matrix_view mview = gsl_matrix_submatrix(temp_input, xmin, ymin,
      cut_xdim, cut_ydim);
  //copy the interesting part of the matrix
  gsl_matrix_memcpy(matrix, &(mview.matrix));

  gsl_matrix_free(temp_input);
  return(matrix);
}


/*==========================================================================
 * read file fname (made of two columns) and save it into a spline object
 * Parameters
 * ----------
 *  fname: string
 *    name of the file to read
 * output
 * ------
 *  spl: gsl spline object
 *    containing the read file
 *==========================================================================*/
gsl_spline *gslf::read_2_spline(std::string fname){
  std::vector<double> x,y;  //temporary vector to read the file
  
  x.push_back(0);   //add point P(k=0)=0 to the input file: avoid errors when the spline is evaluated
  y.push_back(0);   //and shouldn't influence the output 
  
  std::ifstream in(fname.c_str(), std::ifstream::in);   //open the file object
  if(!in.is_open()){
    std::cerr << "Could not open file " << fname << std::endl;
    exit(84);
  }
  double xi, yi; //temporary columns
  while(in >> xi >> yi){
    x.push_back(xi);
    y.push_back(yi);
  }
  in.close();

  x.push_back(x[x.size()-1] *100. );   //add point P(k=very large)=0 to the input file: avoid errors when the spline is evaluated
  y.push_back(0);   //and shouldn't influence the output 

  gsl_spline *spl = gsl_spline_alloc(gsl_interp_linear, x.size());
  if(gsl_spline_init(spl, &x[0], &y[0], x.size()) !=0){
    std::cerr << "Initialization of the spline " << gsl_spline_name(spl);
    std::cerr << " for the data in file " << fname << " failed" << std::endl;
    exit(85);
  }
  x.clear();
  y.clear();
  return spl;  //return the spline
}
