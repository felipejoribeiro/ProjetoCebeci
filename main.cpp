// Foi utilizada a biblioteca boost, para se ter types de grande precisão. Possível achar em https://www.boost.org

#include <boost/lambda/lambda.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include "functions.h"
using boost::multiprecision::cpp_dec_float_50;


int main()
{
    cpp_dec_float_50 variavel = cpp_dec_float_50(1) / 7;

    variavel = variavel ;
    std::cout.precision(std::numeric_limits<cpp_dec_float_50>::digits10);
    std::cout << variavel << std::endl;

    return 0;

}

