!variaveis globais definidas (para que as funções que as tenham)
module prandtll
    double precision , dimension(:) , allocatable :: vPrt , T , u , e
    integer , dimension(3) :: p
    double precision:: Ret , Re, Pr , Prt , vc , dy , incre , um
    integer:: N , increlog
    character*200 :: dirname
end module


! Programa principal

program teste

    use prandtll
    implicit none

    ! Declaram-se as variáveis do programa

    double precision :: R

    ! Parameters

    Ret = 395.d0                                                                       ! 150.d0  !180.d0   !395.d0   ! 640.d0    !1020.d0
    Pr = 10.d0                                                                         !0.71d0   !10.d0

    ! Controles numéricos

    N = 100                                                                            ! Número de células
    incre = 1.d-9                                                                      ! incremento para convergência do método implícito
    increlog = int(- log10(incre))                                                     ! Número de casas após a vírgula.
    R = 1.d0                                                                           ! Raio do canal
    dy = (R/(dble(N) - 0.5d0)) * Ret/R;                                                ! i_1 = dy/2 ... i_n = R

    ! Méta modelos a partir da referência

    prt = ((1.3d0 * 10.d0 ** (-11.d0) ) * Ret**3 - &
    (7.1d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + 0.0001d0 * Ret + 0.87d0)* (pr/0.71)**(0.000006)

    vc = (Ret**(log(Ret) * 0.045d0) * exp(5.3) ) / (Ret ** 0.61d0)

    ! Adequação aos parâmetros padrão

    call AdequaParametro()

    ! Adequação numérica final (usuário)

    !!...

    ! Alocando-se os alocáveis

    allocate(e(N))
    allocate(u(N))

    ! Desenvolvimento do método

    call Program()

    ! Desalocando-se os desalocáveis

    deallocate(e)
    deallocate(u)

end program




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
    ! Simulação do vetor temperatura
    call TemperatureSimu()






    end subroutine Program




! Tem como objetivo setar os parametros de forma devida, de acordo com o DNS
! p(1) = Número de células no DNS da temperatura.
! p(2) = Número de células no DNS do prandtl.
! p(3) = Número de células no DNS da velocidade.
! Adequar de acordo com os DNS's...
subroutine AdequaParametro()

    use prandtll
    if(Ret == 1020.d0)then
        p(1) = 224
        p(2) = p(1)
        p(3) = P(1)
        Re = 41441.d0
    elseif(Ret == 150.d0)then
        p(1) = 73
        p(2) = p(1)
        p(3) = 64
        Re = 4560.d0
    elseif(Ret == 640.d0)then
        p(1) = 128
        p(2) = p(1)
        p(3) = p(1)
        Re = 24428.d0
    elseif(Ret == 180.d0)then
        p(1) = 64
        p(2) = p(1)
        p(3) = p(1)
        Re = 5683.d0
    elseif(Ret == 395.d0)then
        p(1) = 240
        p(2) = 96
        p(3) = p(2)
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


! Desenvolve o vetor temperatura
subroutine TemperatureSimu()

    use prandtll
    implicit none
    return

    end subroutine TemperatureSimu

! Do geito