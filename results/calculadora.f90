
program calculadora
    implicit none
    double precision:: prt , Ret , pr , vc

pr = 0.71d0
Ret = 395.d0



prt = ((4.52901632 * 10.d0 ** (-12.d0) ) * Ret**3 - &
(5.73952059d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + &
(9.397008473d0 * 10.d0 ** (-5.d0) )* Ret + 0.873117480)* (pr/0.71)**(-0.008d0)

vc = exp( 0.164405721012d0 * log(Ret)**3.d0 - 2.87424334318d0 * log(Ret)**2.d0 +  16.3562873171d0 * log(Ret) - &                      ! genetic cebeci com ajuste molecular
            26.6310370449d0 )


print*, prt , vc

end program