!variaveis globais definidas (para que as funções que as tenham)
module prandtll
    double precision , dimension(:) , allocatable :: vPrt , T , u , e
    double precision , dimension(:,:) , allocatable :: Tdns
    integer , dimension(3) :: p
    double precision:: Ret , Re, Pr , Prt , vc , dy , incre , um, L2 , L1 , Li, acuracialvo, valorInicial, incremento
    integer:: N , SIM
    character*100:: dirname
end module


! Programa principal

program teste

    use prandtll
    implicit none

    ! Declaram-se as variáveis do programa

    double precision :: R

    ! Controle dos parametros

        Ret = 150.d0
        Pr = 0.71d0


        ! Controles numéricos

        incremento = 0.1d0                                                                 ! incremento inicial
        acuracialvo = 0.001d0                                                              ! Acurácia para os números de Prandtl turbulento ideais.
        valorInicial = 0.5d0                                                               ! Valor do qual o Prandtl turbulento ira iniciar.
        N = 400                                                                            ! Número de células.
        incre = 1.d-9                                                                      ! incremento para convergência do método implícito.
        R = 1.d0                                                                           ! Raio do canal.
        dy = (R/(dble(N) - 0.5d0)) * Ret/R;                                                ! i_1 = dy/2 ... i_n = R

        ! Méta modelos a partir da referência

        vc = (Ret**(log(Ret) * 0.04510621d0) * exp(5.27528132d0) ) / (Ret ** 0.60941173d0)

        vc = 26

        ! Adequação aos parâmetros padrão
        call AdequaParametro()



        !inicia o algorítimo evolutivo
        call evolutivo()


        ! Amostrando resultados
        print*, "------------------------------------------------------------------------------"
        Print*, "Fim da simulação !"
        Print*, "Ret =" , Ret
        Print*, "Pr =" , Pr
        Print*, "Valor de Cebeci =" , vc
        Print*, "N =", N
        print*,"Prandtl turbulento ideal encontrado : " , Prt

end program




! Algorítimo evolutivo
subroutine evolutivo()

    use prandtll
    implicit none
    character(Len = 200) :: nome1 , nome2 , nome3 , nome4, nome5
    double precision :: const1 , const2 , const3
    write(nome1 , "(f8.2)") Ret
    write(nome2 , "(f4.2)") Pr
    write(nome3 , "(I5)") N
    print*, "Algoritmo de otimização iniciado!"
    print*, "#######################################"
    print*, "# Reynolds tau = " , trim(nome1) , "             #"
    print*, "# Prandtl molecular = " , trim(nome2), "            #"
    print*, "# número de células = " , trim(nome3) , "           #"
    print*, "#######################################"
    SIM = 0
    print*, "---------------------------------------------------------"
    Prt = valorInicial
    ! Primeira simulação
    call simulacao()
    const1 = L2
    Prt = Prt + incremento
    call simulacao()
    const2 = L2
    Prt = Prt - 2.d0 * incremento
    call simulacao()
    const3 = L2
    Prt = Prt + incremento
    print*, "----------"
        do while(incremento > acuracialvo)
        if(const1 < const2 .and. const1 < const3)then
        incremento = incremento/2.d0
        Prt = Prt + incremento
        call simulacao()
        const2 = L2
        Prt = Prt - 2.d0 * incremento
        call simulacao()
        const3 = L2
        write(nome5 , "(f4.2)") Prt
        Prt = Prt + incremento
        print*, "Prt = " ,  trim(nome5) , " incremento =" , incremento
        print*, "----------"
        elseif(const1 < const2 .and. const1 > const3)then
        Prt = Prt - incremento
        const2 = const1
        const1 = const3
        Prt = Prt - incremento
        call simulacao()
        const3 = L2
        Prt = Prt + incremento
        print*, "<--"
        print*, "----------"
        elseif(const1 > const2 .and. const1 < const3)then
        Prt = Prt + incremento
        const3 = const1
        const1 = const2
        Prt = Prt + incremento
        call simulacao()
        const2 = L2
        Prt = Prt - incremento
        print*, "-->"
        print*, "----------"
        end if
        end do
    end subroutine evolutivo




!Simulação numérica
subroutine simulacao()

    use prandtll
    implicit none
    character(Len = 200) :: nome1 , nome2 , nome3
    SIM = SIM + 1
    if(SIM < 10)then
    Write(nome1, "(I1)") SIM
    elseif(SIM >=10 .and. SIM <100)then
    Write(nome1, "(I2)") SIM
    elseif(SIM >=100 .and. SIM <1000)then
    Write(nome1, "(I3)") SIM
    end if
    print*, "simulação: " , trim(nome1)
    ! Alocando-se os alocáveis
    allocate(e(N))
    allocate(u(N))
    allocate(vPrt(N))
    allocate(T(N))
    allocate (Tdns(2 , p(1)))
    ! Desenvolvimento do método
    call Program()
    ! Desalocando-se os desalocáveis
    deallocate(Tdns)
    deallocate(T)
    deallocate(vPrt)
    deallocate(e)
    deallocate(u)
    Write(nome2, "(f4.2)") L2
    print*, " L2: " , trim(nome2)
end subroutine simulacao





! Desenvolvimento computacional
subroutine Program()

    use prandtll
    implicit none
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






subroutine PrandtlvectorDNS640()

    use prandtll
    implicit none
    integer :: i , ii , iii
    double precision :: k1 , k2 , k3 , ly
    double precision , dimension(:,:) , allocatable :: TPrt
    character(len=:), allocatable :: m
    allocate(TPrt(2 , 128))
    allocate(character(len=len(dirname // '/DNS/Prt_RE_640_071.txt')) :: m)
    m = trim(dirname) // '/DNS/Prt_RE_640_071.txt'
    open(unit=10,file= m)
    Do i = 1 ,  128
            read(10,*) ii , TPrt(1, i), k1 , k2 , k3 , TPrt(2, i)
    End Do
    close(10)
    k1 = 0
    k2 = 0
    do i = 2 , 128
        do ii = 1 , N
        if( e(ii) > ( TPrt(1, i-1)) .and. e(ii) < ( TPrt(1, i)) )then
            k1 = (TPrt(2,i) - TPrt(2,i-1))/(TPrt(1,i) - TPrt(1, i-1))
            vPrt(N + 1 - ii) = TPrt(2, i-1) + k1 * (e(ii) - TPrt(1, i-1))
        end if
     end do
    end do
    vPrt(1) = vPrt(2)
    deallocate(m)
    return

    end subroutine PrandtlvectorDNS640






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
    end do
    return

    end subroutine TemperatureSimu




    ! Desenvolve o vetor temperatura
subroutine TemperatureSimu2()

    use prandtll
    implicit none
    integer :: i
    double precision :: k1 , k2 , k3 , f , somatoria
    do i = 1 , N
        T(i) = 0.d0
    end do
    do i = N , 2 , -1

        k1 = (somatoria(i))/(f(e(i)-dy , i-2)*um)
        k2 = (somatoria(i))/(f(e(i)-0.5d0*dy , i - 1)*um)
        k3 = (somatoria(i))/(f(e(i) , i )*um)

        T(i-1) = T(i) - dy * (k1 + 2.d0 * k2 + k3)/4.d0
    end do

    print*, T
    return

    end subroutine TemperatureSimu2



    !Desenvolvimento da integral velocidade
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




! Tirando norma L1 da temperatura
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