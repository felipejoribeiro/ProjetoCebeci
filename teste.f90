!variaveis globais definidas (para que as funções as tenham)

module prandtll
    double precision , dimension(:) , allocatable :: vPrt
    double precision:: Ret , Re, Pr , Prt , vs
end module


! Programa principal:

program teste

    use prandtll
    implicit none

    !! Declaram-se as variáveis do programa:
    
    integer :: N , i , ii , iii, iiii, temperatura , prandtlanal, prandtlvalor, acree, iiiii
    double precision :: R , dy , Nu, k1 , k2 , k3, k4 , ly , ff , L , um , f , prrm, vprrm, acre, ul
    double precision , dimension(:) , allocatable :: T , u , e 
    double precision , dimension(:, :) , allocatable :: s , us, Prr
    character*200 :: dirname , m, NN , pp, RR, filename, filenamet , aa

    ! Parameters: 
    
    Ret = 395.d0   ! 150.d0  !180.d0   !395.d0   ! 640.d0    !1020.d0         
    pr = 10.d0     !0.71d0   !10.d0
    Re =  14062.d0 !4560.d0  !5683.d0  !14062.d0 ! 24428.d0  !41441.d0
    Nu = 83.00d0 
    
    ! Controles numéricos:

    N = 100
    acre = 1.d-9
    R = 1.d0
    dy = (R/(dble(N) - 0.5d0)) * Ret/R;    ! i_1 = dy/2 ... i_n = R 
    acree = int(- log10(acre))

    ! Méta modelos:

    prt = ((1.3d0 * 10.d0 ** (-11.d0) ) * Ret**3 - & 
    (7.1d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + 0.0001d0 * Ret + 0.87d0)* (pr/0.71)**(0.000006)

    vs = (Ret**(log(Ret) * 0.045d0) * exp(5.3) ) / (Ret ** 0.61d0)

    ! Criando vetor espaço:

    print*, prt , vs 





end program
