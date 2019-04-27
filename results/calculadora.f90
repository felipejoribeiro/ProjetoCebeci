
program calculadora
    implicit none
    double precision:: prt , Ret , pr , vc
    character*200:: metodo , dirname
    integer :: N
Call getcwd( dirname )
N = 400
pr = 0.71d0
Ret = 150.d0
! metodo = '_Prt0905_A26'


prt = ((4.52901632 * 10.d0 ** (-12.d0) ) * Ret**3 - &
(5.73952059d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + &
(9.397008473d0 * 10.d0 ** (-5.d0) )* Ret + 0.873117480)* (pr/0.71)**(-0.008d0)

vc = exp( 0.164405721012d0 * log(Ret)**3.d0 - 2.87424334318d0 * log(Ret)**2.d0 +  16.3562873171d0 * log(Ret) - &                      ! genetic cebeci com ajuste molecular
            26.6310370449d0 )


! Prt = - 4.56041707672d0 * 10.d0 ** (-10.d0) * Ret**3 +  9.56902551372d0 * 10.d0 **(-7.d0) * Ret**2 &                               ! Otimizado sem a otimização de cebeci
! - 0.000617158206068d0 * Ret +  1.01789506426

prt = (Ret**(0.164405721* (log(Ret))**2 - 2.874243343* log(Ret) + 16.356287317 ) )/(exp(26.631037045))






Prt = -2.48916601371e-10 * Ret**3 +  3.60362337151e-07 * Ret**2 +  3.79213671785e-05 * Ret +  0.71234674305                                 ! genetic with 2 temperature



vc = exp( 0.0395059904287 * log(Ret)**3 -0.758759596012 * log(Ret)**2  + 4.66369525666 * log(Ret) -5.6703426304 )                    ! genetic with 2 temperatures cebeci








print*, prt , vc

    ! write(*,FMT=201) trim(dirname) , Ret, Pr , N , trim(metodo)
    ! 201     format( A ,'/results/graficos/image',F5.0,'_', F5.2,'_', I3 , A ,".txt")
end program