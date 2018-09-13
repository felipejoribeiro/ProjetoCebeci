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
    cpp_dec_float_50 variavel1; // declarou-se uma variável com alta resolução
    cpp_dec_float_50 variavel2; // declarou-se outra variável com alta precisão.

    variavel1 = cpp_dec_float_50(1) * 1/7;

    variavel2 = cpp_dec_float_50(1) * 1/3;

    std::cout.precision(std::numeric_limits<cpp_dec_float_50>::digits10);
    std::cout << variavel1 << std::endl;

    variavel1 = cpp_dec_float_50(4) + variavel2;

    std::cout << variavel1 << std::endl;

    return 0;

}

