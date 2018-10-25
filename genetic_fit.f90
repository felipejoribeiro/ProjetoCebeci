! Esse programa executa um fit a partir de otimização.
module data_HDE                                                                                              ! módulo com parametros base.
    implicit none
    integer :: NP, itermax , strategy, refresh, iwrite , Dim_XC
    integer, dimension(3) :: method
    double precision :: VTR , CR_XC
    double precision :: F_XC , F_CR
    double precision, dimension(:) , allocatable :: XCmin, XCmax
    double precision, dimension(:) , allocatable :: bestmem_XC
    integer :: nfeval
    double precision :: bestval
end module data_HDE

module prandtll
    double precision , dimension(:) , allocatable :: vPrt , T , u , e
    double precision , dimension(:,:) , allocatable :: Tdns
    integer , dimension(3) :: p
    double precision:: Ret , Re, Pr , Prt , vc , dy , incre , um, L2 , L1 , Li
    integer:: N
    character*100:: dirname , filename
end module



program Rosen
use data_HDE
use prandtll
implicit none
integer :: i
integer , dimension(8) :: time
intrinsic date_and_time
external FTN
!...........................................................................................................
! Setagem de parametros do algorítmo evolutivo.                                                            .
!...........................................................................................................
! Parametros
NP=10                                                                        ! Indivíduos na população
itermax=100                                                                  ! Número máximo de iterações
refresh=4                                                                    ! Intervalos de atualização do console
VTR=-1.0d-3                                                                  ! Acurácia esperada do argumento de saída
iwrite=7                                                                     ! Unidade de escrita de arquivo.
Dim_XC=1                                                                     ! Número de variáveis de entrada
filename = "/results/FITPRANDTL.txt"                                         ! Nome do arquivo.
N = 400                                                                      ! Número de células
allocate(XCmin(Dim_XC) , XCmax(Dim_XC), bestmem_XC(Dim_XC))
! Parametros metodológicos
method=(/0, 1, 0/)                                                           ! Metodologia de mutação
strategy=6                                                                   ! Stratégias de mutação
CR_XC=0.5d0                                                                  ! Crossover factor for real decision parameters.
F_XC=0.8d0                                                                   ! Mutation scaling factor for real decision parameters.
F_CR=0.8d0                                                                   ! Random combined factor
XCmin(1)= -10.0d0                                                            ! Limite inferior para domínio de algorítmos de entrada
XCmax(1)=10.0d0                                                              ! Limite superior para domínio de algorítmos de entrada
!...........................................................................................................
! Rotina de escrita em arquivo resultante                                                                  .
!...........................................................................................................

open(iwrite,file='Resultados.txt')
call date_and_time(values=time)                                                                               ! Chamou tempo

write(unit=iwrite, FMT=11) time(1:3), time(5:7)
11  format(1x, 'Beginning of Program. Date:', I4, '/', I2,'/', I2, ', Time:', 1x , I2,':',I2,':',I2)          ! Escrita tempo e hora de início.

print*, filename

call DE_Fortran90(FTN, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC,&
                CR_XC, strategy, refresh, iwrite, bestmem_XC, &                                               ! Chamou o algorítmo evolutivo.
                bestval, nfeval, F_CR, method)

write(iwrite,205) NP, nfeval, method(1:3)                                                                     ! Escrita Número da população, número de chamadas de função e metodologia mutagênica
205     format(2x, 'NP=', I4, 4x, 'No. function call =', I9, /2x, 'mehtod(1:3) =',3I3)

write(iwrite,FMT=201) F_XC, CR_XC, F_CR                                                                       ! Escrita Fatores de escala.
201     format(2x, 'F_XC =',F6.3, 2x, 'CR_XC =', F6.3, 2x, 'F_CR =', F6.3)

write(iwrite,FMT=200) bestval                                                                                 ! Escrita Melhor valor de função (em termos de variável de retorno)
200     format(/2x, 'Bestval=', ES14.7)

do i=1,Dim_XC
    write(iwrite,FMT=202) i,bestmem_XC(i)                                                                     ! Escrita Melhor parametro (variaveis de controle)
end do
202     format(2x, 'best_XC(',I3,') =',ES14.7)

call date_and_time(values=time)                                                                               ! Chamou tempo

write(unit=iwrite, FMT=10)time(1:3), time(5:7)                                                                ! Escrita tempo e hora de final
10  format(/1x, 'End of Program. Date:', I4, '/', I2,'/', I2, ', Time: ', I2,':',I2,':',I2)

deallocate(XCmin, XCmax, bestmem_XC)


end program Rosen








! Funçao a ser otimizada! Aqui entra o código!
subroutine FTN(X, objval)
  use data_HDE
  use prandtll
  implicit none
  double precision, dimension(Dim_XC), intent(in) :: X
  double precision, intent(out) :: objval
  double precision :: re1 , re2 , re3
  integer :: i


  ! Declaram-se as variáveis do programa

    Ret = 395.d0                                                                 ! Reynolds turbulento
    Pr = 10.d0                                                                   ! Prandtl molecular

        ! Méta modelos a partir da referência

        prt = (3.19791882062d-10 * Ret**3 - 1.08216023658d-06 * Ret**2 +0.00116281300928*Ret+0.449206978959)* &
        ((Pr/0.71)**(-0.008d0) + x(1) * (Pr - 0.71))

        ! Amostrando resultados

        print*, "Prt (10) = " , prt
        re1 = prt



        Pr = 0.025d0                                                                        ! Prandtl molecular

        ! Méta modelos a partir da referência

        prt = (3.19791882062d-10 * Ret**3 - 1.08216023658d-06 * Ret**2 +0.00116281300928*Ret+0.449206978959)* &
        ((Pr/0.71)**(-0.008d0) + x(1) * (Pr - 0.71))

        ! Amostrando resultados

        print*, "Prt (0.025) = " , prt
        re2 = prt



        Pr = 0.71d0                                                                       ! Prandtl molecular

        ! Méta modelos a partir da referência

        prt = (3.19791882062d-10 * Ret**3 - 1.08216023658d-06 * Ret**2 +0.00116281300928*Ret+0.449206978959)* &
        ((Pr/0.71)**(-0.008d0) + x(1) * (Pr - 0.71))


        ! Amostrando resultados

        print*, "Prt (0.71) = " , prt
        re3 = prt



        ! Encerramento

        ! print*, " "

        write(13, *) "resultados:" , re1 ,  re2 , re3
        print*, "Fim"

  objval = (re1 - 0.88268260782709396d0)**2 + (re2 -   0.92602145583845075d0 )**2 + (re3 - 0.75938280042931627d0)**2

end subroutine FTN






! Desenvolvimento computacional
subroutine Program()

    use prandtll
    ! Identificando diretório atual
    Call getcwd( dirname )
    ! Criação do vetor espaco discretizado
    call SpaceVector()
    ! Simulação do vetor velocidade e valor médio
    call VelocitySimu()
    call VelocityMedia()
    ! Setagem do vetor Prandtl turbulento
    call Prandtlvector()
    ! Simulação do vetor temperatura
    call TemperatureSimu()
    ! Importando o DNS
    Call DNSinput()
    ! Tirando norma L2
    call L2norm()
    ! Tirando norma L1
    call L1norm()
    ! Tirando norma Li
    call Linorm()
    return

    end subroutine Program






! Tem como objetivo setar os parametros de forma devida, de acordo com o DNS
! p(1) = Número de células no DNS da temperatura.
! p(2) = Número de células no DNS do prandtl.
! p(3) = Número de células no DNS da velocidade.
! Adequar de acordo com os DNS's...
subroutine AdequaParametro()

    use prandtll
    if(Ret == 1020.d0 .and. Pr == 0.71d0)then
        p(1) = 224
        Re = 41441.d0
    elseif(Ret == 150.d0 .and. Pr == 0.71d0)then
        p(1) = 73
        Re = 4560.d0
    elseif(Ret == 150.d0 .and. Pr == 0.025d0)then
        p(1) = 73
        Re = 4560.d0
    elseif(Ret == 640.d0 .and. Pr == 0.71d0)then
        p(1) = 128
        Re = 24428.d0
    elseif(Ret == 640.d0 .and. Pr == 0.025d0)then
        p(1) = 128
        Re = 24428.d0
    elseif(Ret == 395.d0 .and. Pr == 10.d0)then
        p(1) = 240
        Re = 14062.d0
    elseif(Ret == 395.d0 .and. Pr == 0.71d0)then
        p(1) = 96
        Re = 14062.d0
    elseif(Ret == 395.d0 .and. Pr == 7.d0)then
        p(1) = 240
        Re = 14062.d0
    elseif(Ret == 395.d0 .and. Pr == 5.d0)then
        p(1) = 240
        Re = 14062.d0
    elseif(Ret == 395.d0 .and. Pr == 2.d0)then
        p(1) = 240
        Re = 14062.d0
    elseif(Ret == 395.d0 .and. Pr == 1.d0)then
        p(1) = 240
        Re = 14062.d0
    elseif(Ret == 395.d0 .and. Pr == 0.025d0)then
        p(1) = 96
        Re = 14062.d0
    end if
    return

    end subroutine AdequaParametro






! Desenvolvimento do vetor espaco
subroutine SpaceVector()

    use prandtll
    implicit none
    integer :: i
    do i = 1 , N
        e(i) = (dble(i) - 0.5d0)*dy
    end do
    return

    end Subroutine SpaceVector






! Desenvolve o vetor velocidade
subroutine VelocitySimu()

    use prandtll
    implicit none
    double precision :: k1 , k2 , k3 , ff
    integer :: i
    u(N) = 0.d0
    do i = N , 2 , -1

        k1 = ff(e(i) - dy)
        k2 = ff(e(i) - 0.5d0*dy)
        k3 = ff(e(i))

        u(i-1) = u(i) - dy * (k1 + 2.d0 * k2 + k3)/4.d0  ! $u_{i-1} = u_{i} * \frac{y ( \rho_1 + 2 \rho_2 + \rho_3 )}{4} $
    end do
    return

    end subroutine VelocitySimu





! Calcula velocidade média da velocidade
subroutine VelocityMedia()

    use prandtll
    implicit none
    integer :: i
    um = 0.d0
    do i = 1 , N
        um = um + u(i) !* dy
    end do
    um = um / N  !(Ret)
    return

    end subroutine VelocityMedia






! du/dy(y)
function ff(position)

    use prandtll
    implicit none
    double precision :: ff , L
    double precision, intent(in) :: position
    ff = ( -2.d0 * position * (1.d0/Ret)/( 1.d0 + sqrt(1.d0 + 4.d0 * (L(position))**2.d0 *Ret * position)))
    return

    end function ff






! L(Y)
function L(position)

    use prandtll
    implicit none
    double precision :: L
    double precision, intent(in) :: position
    L = ((0.14d0 - 0.08d0 * (position/Ret)**2.d0 - 0.06d0*(position/Ret)**4.d0 )*(1.d0 - exp((position/Ret - 1.d0)*Ret/vc)))
    return

    end function L






! Setagem do vetor Prandtl turbulento
subroutine Prandtlvector()

    use prandtll
    implicit none
    integer :: i
    do i= 1 , N
        vPrt(i) = prt
    end do
    return

    end subroutine Prandtlvector






! Desenvolve o vetor temperatura
subroutine TemperatureSimu()

    use prandtll
    implicit none
    integer :: i
    double precision :: k2 , f
    do i = 1 , N
        T(i) = 0.d0
    end do
    k2 = 10.d0
    do while ( abs( T(1) - k2  ) > incre )
    k2 = T(1)
    T(1) = T(2) + ( dy**2.d0 *(u(1)/um ))/f(e(1) + dy/2.d0 , 1)
        do i = 2 , N - 1
            T(i) =( (dy**2.d0)*(u(i)/um) &
            + T(i-1)*f(e(i) - dy/2.d0 , i - 1) + T(i+1)*f(e(i) + dy/2.0d0, i))/(f(e(i)-dy/2.d0,i-1) &
            + f(e(i) + dy/2.d0, i) )
        end do
        ! print*, T(1)
    end do
    return

    end subroutine TemperatureSimu






function somatoria(i)

    use prandtll
    implicit none
    integer :: i , ii
    double precision :: soma , somatoria
    soma = 0.d0
    do ii = 1 , i
        soma = soma + u(i) * dy
    end do
    somatoria = soma
    return

    end function






! f(Y)
function f(s, i)

    use prandtll
    implicit none
    double precision, intent(in) :: s
    integer, intent(in) :: i
    double precision :: f , ff , L
    f = ( Ret/Pr - ((((L(s))**2.d0 )*Ret**3.d0)/vPrt(i) ) * ff(s)  )
    return

end function f






! importanto o DNS
subroutine DNSinput()

    use prandtll
    implicit none
    integer :: i
    double precision :: ly
    character(len=:), allocatable :: m
    ! Abertura de arquivo
    if (Ret == 1020.d0 .and. Pr == 0.71d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_1000_071.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_1000_071.txt'
        open(unit=10,file= m)
    elseif(Ret == 150.d0 .and. Pr == 0.71d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_150_071.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_150_071.txt'
        open(unit=10,file=m)
    elseif(Ret == 150.d0 .and. Pr == 0.025d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_150_0025.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_150_0025.txt'
        open(unit=10,file=m)
    elseif(Ret == 640.d0 .and. Pr == 0.71d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_640_071.txt')):: m)
        m = trim(dirname) // '/DNS/DNS_RE_640_071.txt'
        open(unit=10,file=m)
    elseif(Ret == 640.d0 .and. Pr == 0.025d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_640_0025.txt')):: m)
        m = trim(dirname) // '/DNS/DNS_RE_640_0025.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0 .and. Pr == 0.71d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395_071.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395_071.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0 .and. Pr == 10.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395_10.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395_10.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0 .and. Pr == 7.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395_7.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395_7.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0 .and. Pr == 5.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395_5.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395_5.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0 .and. Pr == 2.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395_2.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395_2.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0 .and. Pr == 1.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395_1.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395_1.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0 .and. Pr == 0.025d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395_0025.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395_0025.txt'
        open(unit=10,file=m)
    end if
    ! Leitura
    Do i = 1 ,  p(1)
            read(10,*) ly , Tdns(1, i), Tdns(2, i)
    End Do
    close (10)
    deallocate(m)
    return

    end subroutine DNSinput






! Tira espaços do meio de string
subroutine SOUT(string)

    character(len=*) :: string
    integer :: stringLen
    integer :: last, actual
    stringLen = len (string)
    last = 1
    actual = 1
    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do
    return

    end subroutine SOUT






! Tira a norma L2 da temperatura
subroutine L2norm()

    use prandtll
    implicit none
    double precision :: ly , k1
    integer :: i , ii , iii
    k1 = 0.d0
    iii = 0
    do i = 2 , N

        do ii = 1 , p(1)
            if((Ret - Tdns(1 , ii)) < e(i) .and. (Ret - Tdns(1 , ii)) > e(i-1) )then
                ly = (T(i) - T(i-1))*((Ret - Tdns(1,ii)) - e(i-1))/(e(i) - e(i-1))
                ly = T(i-1) + ly
                k1 = k1 + (ly - Tdns(2,ii))**2.d0
                iii = iii + 1
            end if
        end do

    end do
    L2 = sqrt(k1 / (iii))
    return

    end subroutine L2norm






subroutine L1norm()

    use prandtll
    implicit none
    double precision :: ly , k1
    integer :: i , ii , iii
    k1 = 0.d0
    iii = 0
    do i = 2 , N
        do ii = 1 , p(1)
            if((Ret - Tdns(1 , ii)) < e(i) .and. (Ret - Tdns(1 , ii)) > e(i-1) )then
                ly = (T(i) - T(i-1))*((Ret - Tdns(1,ii)) - e(i-1))/(e(i) - e(i-1))
                ly = T(i-1) + ly
                k1 = k1 + abs(ly - Tdns(2,ii))
                iii = iii + 1
            end if
        end do
    end do
    L1 = k1 / iii
    return

    end subroutine L1norm






! Tirando a L infinito da temperatura
subroutine Linorm()

    use prandtll
    implicit none
    double precision :: ly , k1
    integer :: i , ii , iii
    k1 = 0.d0
    iii = 0
    do i = 2 , N
        do ii = 1 , p(1)
            if((Ret - Tdns(1 , ii)) < e(i) .and. (Ret - Tdns(1 , ii)) > e(i-1) )then
                ly = (T(i) - T(i-1))*((Ret - Tdns(1,ii)) - e(i-1))/(e(i) - e(i-1))
                ly = T(i-1) + ly
                if(abs(ly - Tdns(2,ii)) > k1)then
                    k1 = abs(ly - Tdns(2,ii))
                end if
            end if
        end do
    end do
    Li = k1
    return

    end subroutine Linorm







! Algorítmo de otimização
subroutine DE_Fortran90(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC, &
           CR_XC, strategy, refresh, iwrite, bestmem_XC, bestval, nfeval, &
           F_CR, method)
!.......................................................................
!
! Differential Evolution for Optimal Control Problems
!
!.......................................................................
!  This Fortran 90 program translates from the original MATLAB
!  version of differential evolution (DE). This FORTRAN 90 code
!  has been tested on Compaq Visual Fortran v6.1.
!  Any users new to the DE are encouraged to read the article of Storn and Price.
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of the
!    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary
!    Comutation, 842-844.
!
!  This Fortran 90 program written by Dr. Feng-Sheng Wang
!  Department of Chemical Engineering, National Chung Cheng University,
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!.........................................................................
!                obj : The user provided file for evlauting the objective function.
!                      subroutine obj(xc,fitness)
!                      where "xc" is the real decision parameter vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the real decision parameters.
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
!                VTR : The expected fitness value to reach.
!                 NP : Population size.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for real decision parameters.
!              CR_XC : Crossover factor for real decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced if
!                      "refresh < 1".
!             iwrite : The unit specfier for writing to an external data file.
! bestmen_XC(Dim_XC) : The best real decision parameters.
!              bestval : The best objective function.
!             nfeval : The number of function call.
!         method(1) = 0, Fixed mutation scaling factors (F_XC)
!                   = 1, Random mutation scaling factors F_XC=[0, 1]
!                   = 2, Random mutation scaling factors F_XC=[-1, 1]
!         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
!                        in the mutation operation
!                   = other, fixed combined factor provided by the user
!         method(3) = 1, Saving results in a data file.
!                   = other, displaying results only.
use prandtll
implicit none
integer, intent(in) :: NP, Dim_XC, itermax, strategy, iwrite, refresh
double precision, intent(in) :: VTR, CR_XC
double precision :: F_XC, F_CR
double precision, dimension(Dim_XC), intent(in) :: XCmin, XCmax
double precision, dimension(Dim_XC), intent(inout) :: bestmem_XC
double precision, intent(out) :: bestval
integer, intent(out) :: nfeval
double precision, dimension(NP,Dim_XC) :: pop_XC, bm_XC, mui_XC, mpo_XC, popold_XC, rand_XC, ui_XC
integer :: i, ibest, iter
integer, dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
integer, dimension(4) :: ind
double precision :: tempval
double precision, dimension(NP) :: val
double precision, dimension(Dim_XC) :: bestmemit_XC
double precision, dimension(Dim_XC) :: rand_C1
integer, dimension(3), intent(in) :: method
external  obj
intrinsic max, min, random_number, mod, abs, any, all, maxloc
interface
function randperm(num)
implicit none
integer, intent(in) :: num
integer, dimension(num) :: randperm
end function randperm
end interface
call getcwd(dirname)
open(unit=13 , file = trim(dirname)//trim(filename) )
!!-----Initialize a population --------------------------------------------!!

pop_XC=0.0d0
do i=1,NP
  call random_number(rand_C1)
  pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
  end do

!!--------------------------------------------------------------------------!!

!!------Evaluate fitness functions and find the best member-----------------!!
val=0.0d0
nfeval=0
ibest=1
call obj(pop_XC(1,:), val(1))
bestval=val(1)
nfeval=nfeval+1
do i=2,NP
  call obj(pop_XC(i,:), val(i))
  nfeval=nfeval+1
  if (val(i) < bestval) then
      ibest=i
      bestval=val(i)
      end if
  end do

bestmemit_XC=pop_XC(ibest,:)
bestmem_XC=bestmemit_XC
!!--------------------------------------------------------------------------!!

bm_XC=0.0d0
rot=(/(i,i=0,NP-1)/)
iter=1
!!------Perform evolutionary computation------------------------------------!!
do while (iter <= itermax)
  popold_XC=pop_XC

!!------Mutation operation--------------------------------------------------!!
  ind=randperm(4)
  a1=randperm(NP)
  rt=mod(rot+ind(1),NP)
  a2=a1(rt+1)
  rt=mod(rot+ind(2),NP)
  a3=a2(rt+1)
  rt=mod(rot+ind(3),NP)
  a4=a3(rt+1)
  rt=mod(rot+ind(4),NP)
  a5=a4(rt+1)
  bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)
!----- Generating a random sacling factor--------------------------------!
  select case (method(1))
  case (1)
  call random_number(F_XC)
  case(2)
  call random_number(F_XC)
  F_XC=2.0d0*F_XC-1.0d0
  end select

!---- select a mutation strategy-----------------------------------------!
  select case (strategy)
  case (1)
  ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

  case default
  ui_XC=popold_XC(a3,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

  case (3)
  ui_XC=popold_XC+F_XC*(bm_XC-popold_XC+popold_XC(a1,:)-popold_XC(a2,:))

  case (4)
  ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:)-popold_XC(a4,:))

  case (5)
  ui_XC=popold_XC(a5,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:) &
  -popold_XC(a4,:))
  case (6) ! A linear crossover combination of bm_XC and popold_XC
  if (method(2) == 1) call random_number(F_CR)
  ui_XC=popold_XC+F_CR*(bm_XC-popold_XC)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

  end select
!!--------------------------------------------------------------------------!!
!!------Crossover operation-------------------------------------------------!!
  call random_number(rand_XC)
  mui_XC=0.0d0
  mpo_XC=0.0d0
  where (rand_XC < CR_XC)
  mui_XC=1.0d0
!           mpo_XC=0.0d0
  elsewhere
!           mui_XC=0.0d0
  mpo_XC=1.0d0
  end where

  ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC
!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!
  do i=1,NP
!!------Confine each of feasible individuals in the lower-upper bound-------!!
    ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
    call obj(ui_XC(i,:), tempval)
    nfeval=nfeval+1
    if (tempval < val(i)) then
    pop_XC(i,:)=ui_XC(i,:)
    val(i)=tempval
    if (tempval < bestval) then
    bestval=tempval
    bestmem_XC=ui_XC(i,:)
    end if
    end if
    end do
  bestmemit_XC=bestmem_XC
  if( (refresh > 0) .and. (mod(iter,refresh)==0)) then
  if (method(3)==1) write(unit=iwrite,FMT=203) iter
  write(unit=*, FMT=203) iter
  do i=1,Dim_XC
    if (method(3)==1) write(unit=iwrite, FMT=202) i, bestmem_XC(i)
    write(*,FMT=202) i,bestmem_XC(i)
    end do
  if (method(3)==1) write(unit=iwrite, FMT=201) bestval
  write(unit=*, FMT=201) bestval
  end if
  iter=iter+1
  if ( bestval <= VTR .and. refresh > 0) then
  write(unit=iwrite, FMT=*) ' The best fitness is smaller than VTR'
  write(unit=*, FMT=*) 'The best fitness is smaller than VTR'
  exit
  endif
  end do
!!------end the evolutionary computation------------------------------!!
201 format(2x, 'bestval =', ES14.7, /)
202 format(5x, 'bestmem_XC(', I3, ') =', ES12.5)
203 format(2x, 'No. of iteration  =', I8)
write(13, *) bestmem_XC(1) , bestmem_XC(2) , bestval
close(13)
end subroutine DE_Fortran90


function randperm(num)
  implicit none
  integer, intent(in) :: num
  integer :: number, i, j, k
  integer, dimension(num) :: randperm
  double precision, dimension(num) :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
           number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
           number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
end function randperm