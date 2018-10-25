
program calculadora
    implicit none
    double precision:: prt , Ret , pr

pr = 0.025d0
Ret = 395.d0



prt = ((4.52901632 * 10.d0 ** (-12.d0) ) * Ret**3 - &
(5.73952059d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + &
(9.397008473d0 * 10.d0 ** (-5.d0) )* Ret + 0.873117480)* (pr/0.71)**(-0.008d0)




print*, prt

end program