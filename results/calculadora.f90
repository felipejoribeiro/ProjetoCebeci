
program calculadora
    implicit none
    double precision:: prt , Ret , pr , vc
    character*200:: metodo , dirname
    integer :: N
Call getcwd( dirname )
N = 400
pr = 10.d0
Ret = 1020.d0
metodo = '_Prt0905_A26'


prt = ((4.52901632 * 10.d0 ** (-12.d0) ) * Ret**3 - &
(5.73952059d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + &
(9.397008473d0 * 10.d0 ** (-5.d0) )* Ret + 0.873117480)* (pr/0.71)**(-0.008d0)

vc = exp( 0.164405721012d0 * log(Ret)**3.d0 - 2.87424334318d0 * log(Ret)**2.d0 +  16.3562873171d0 * log(Ret) - &                      ! genetic cebeci com ajuste molecular
            26.6310370449d0 )


Prt = - 4.56041707672d0 * 10.d0 ** (-10.d0) * Ret**3 +  9.56902551372d0 * 10.d0 **(-7.d0) * Ret**2 &                               ! Otimizado sem a otimização de cebeci
- 0.000617158206068d0 * Ret +  1.01789506426


print*, prt , vc

    write(*,FMT=201) trim(dirname) , Ret, Pr , N , trim(metodo)
    201     format( A ,'/results/graficos/image',F5.0,'_', F5.2,'_', I3 , A ,".txt")
end program