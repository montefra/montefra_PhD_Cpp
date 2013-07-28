/*==========================================================================
 * Convert 'value' to string. Templatized, so accept all standard types     
 * Parameters
 * ----------
 *  value: value to convert to string
 * output
 * ------
 *  str: string 
 *==========================================================================*/
template <typename T>  //set template name
string to_string(T const& value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}
