!variaveis globais definidas (para que as funções que as tenham)
module prandtll
    double precision , dimension(:) , allocatable :: vPrt , T , u , e
    double precision , dimension(:,:) , allocatable :: Tdns
    integer , dimension(3) :: p
    double precision:: Ret , Re, Pr , Prt , vc , dy , incre , um
    integer:: N
    character(len=:), allocatable :: dirname
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
    R = 1.d0                                                                           ! Raio do canal
    dy = (R/(dble(N) - 0.5d0)) * Ret/R;                                                ! i_1 = dy/2 ... i_n = R

    ! Méta modelos a partir da referência

    prt = ((1.3d0 * 10.d0 ** (-11.d0) ) * Ret**3 - &
    (7.1d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + 0.0001d0 * Ret + 0.87d0)* (pr/0.71)**(-0.04d0)

    vc = (Ret**(log(Ret) * 0.045d0) * exp(5.3) ) / (Ret ** 0.61d0)

    ! Adequação aos parâmetros padrão
    call AdequaParametro()
    ! Adequação numérica final (usuário)
    vc = 25.677d0
    ! Alocando-se os alocáveis
    allocate(e(N))
    allocate(u(N))
    allocate(vPrt(N))
    allocate(T(N))
    allocate (Tdns(2 , p(1)))
    dirname = '                                                 '
    ! Desenvolvimento do método
    call Program()
    ! Desalocando-se os desalocáveis
    deallocate(Tdns)
    deallocate(T)
    deallocate(vPrt)
    deallocate(e)
    deallocate(u)

end program






! Desenvolvimento computacional
subroutine Program()

    use prandtll
    ! Identificando diretório atual
    Call getcwd( dirname )
    dirname = trim(dirname)
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
    call DNSinput()
    ! Tirando norma L2
    call L2norm()



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
        print*, t(1)
    end do
    return

    end subroutine TemperatureSimu




! f(Y)
function f(s, i)

    use prandtll
    implicit none
    double precision, intent(in) :: s
    integer, intent(in) :: i
    double precision :: f , ff , L
    f = ( Ret/Pr - ((((L(s))**2 )*Ret**3)/vPrt(i) ) * ff(s)  )
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
    if (Ret == 1020.d0)then
        m = trim(dirname) // '/DNS/DNS_RE_1000.txt'
        open(unit=10,file= m)
    elseif(Ret == 150.d0)then
        m = trim(dirname) // '/DNS/DNS_RE_150.txt'
        open(unit=10,file=m)
    elseif(Ret == 640.d0)then
        m = trim(dirname) // '/DNS/DNS_RE_640.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0)then
        m = dirname // '/DNS/DNS_RE_395_10.txt'
        open(unit=10,file=m)
    elseif(Ret == 180.d0)then
        m = trim(dirname) // '/DNS/DNS_RE_180.txt'
        open(unit=10,file=m)
    end if
    ! Leitura
    if(Ret == 180.d0)then
        Do i = 1 ,  p(1)
            read(10,*) ly , Tdns(1, i), Tdns(2, i)
        End Do
    else
        Do i = 1 ,  p(1)
            read(10,*) ly , Tdns(1, i), Tdns(2, i)
        End Do
    end if
    close (10)
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




subroutine L2norm()

    use prandtll
    implicit none
    k1 = 0
    k2 = 0
    do i = 2 , N
        do ii = 1 , iii
          if((Ret - s(1 , ii)) < e(i) .and. (Ret - s(1 , ii)) > e(i-1) )then
               ly = (T(i) - T(i-1))*((Ret - s(1,ii)) - e(i-1))/(e(i) - e(i-1))
               ly = T(i-1) + ly
               k1 = k1 + (ly - s(2,ii))**2
               k2 = k2 + 1
          end if
        end do
    end do
    k1 = sqrt(k1 / (k2 - 1))
    return

    end subroutine L2norm


