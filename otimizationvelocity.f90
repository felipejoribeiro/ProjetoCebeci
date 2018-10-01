!variaveis globais definidas (para que as funções que as tenham)
module prandtll
    double precision , dimension(:) , allocatable :: vPrt , u , e
    double precision , dimension(:,:) , allocatable :: Udns
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

        Ret = 1020.d0


        ! Controles numéricos

        incremento = 1.d0                                                                 ! incremento inicial
        acuracialvo = 0.00001d0                                                            ! Acurácia para os números de Prandtl turbulento ideais.
        valorInicial = 20.d0                                                               ! Valor do qual o Prandtl turbulento ira iniciar.
        N = 400000                                                                          ! Número de células.
        incre = 1.d-9                                                                      ! incremento para convergência do método implícito.
        R = 1.d0                                                                           ! Raio do canal.
        dy = (R/(dble(N) - 0.5d0)) * Ret/R;                                                ! i_1 = dy/2 ... i_n = R

        ! Adequação aos parâmetros padrão
        call AdequaParametro()

        !inicia o algorítimo evolutivo
        call evolutivo()


        ! Amostrando resultados
        print*, "------------------------------------------------------------------------------"
        Print*, "Fim da simulação !"
        Print*, "Ret =" , Ret
        Print*, "N =", N
        Print*, "Prt =" , prt
        print*,"Valor de cebeci ideal: " , vc
        print*, "L2 =" , L2

end program




! Algorítimo evolutivo
subroutine evolutivo()

    use prandtll
    implicit none
    character(Len = 200) :: nome1 , nome2 , nome3 , nome4, nome5
    double precision :: const1 , const2 , const3
    write(nome1 , "(f8.2)") Ret
    write(nome3 , "(I8)") N
    print*, "Algoritmo de otimização iniciado!"
    print*, "#######################################"
    print*, "# Reynolds tau = " , trim(nome1) , "             #"
    print*, "# número de células = " , trim(nome3) , "        #"
    print*, "#######################################"
    SIM = 0
    print*, "---------------------------------------------------------"
    vc = valorInicial
    ! Primeira simulação
    call simulacao()
    const1 = L2
    vc = vc + incremento
    call simulacao()
    const2 = L2
    vc= vc - 2.d0 * incremento
    call simulacao()
    const3 = L2
    vc = vc + incremento
    print*, "----------"
        do while(incremento > acuracialvo)
        if(const1 < const2 .and. const1 < const3)then
        incremento = incremento/2.d0
        vc = vc + incremento
        call simulacao()
        const2 = L2
        vc = vc - 2.d0 * incremento
        call simulacao()
        const3 = L2
        write(nome5 , "(f8.4)") vc
        vc = vc + incremento
        print*, "Constante de cebeci = " ,  trim(nome5) , " incremento =" , incremento
        print*, "----------"
        elseif(const1 < const2 .and. const1 > const3)then
        vc = vc - incremento
        const2 = const1
        const1 = const3
        vc = vc - incremento
        call simulacao()
        const3 = L2
        vc = vc + incremento
        print*, "<--"
        print*, "----------"
        elseif(const1 > const2 .and. const1 < const3)then
        vc = vc + incremento
        const3 = const1
        const1 = const2
        vc = vc + incremento
        call simulacao()
        const2 = L2
        vc = vc - incremento
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
    allocate (Udns(2 , p(1)))
    ! Desenvolvimento do método
    call Program()
    ! Desalocando-se os desalocáveis
    deallocate(Udns)
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
    ! Simulação do vetor velocidade
    call VelocitySimu()
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
    if(Ret == 1020.d0)then
        p(1) = 224
        Re = 41441.d0
    elseif(Ret == 150.d0)then
        p(1) = 64
        Re = 4560.d0
    elseif(Ret == 640.d0)then
        p(1) = 128
        Re = 24428.d0
    elseif(Ret == 395.d0)then
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






! importanto o DNS
subroutine DNSinput()

    use prandtll
    implicit none
    integer :: i
    double precision :: ly
    character(len=:), allocatable :: m
    ! Abertura de arquivo
    if (Ret == 1020.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_1000.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_1000.txt'
        open(unit=10,file= m)
    elseif(Ret == 150.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_150.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_150.txt'
        open(unit=10,file=m)
    elseif(Ret == 640.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_640.txt')):: m)
        m = trim(dirname) // '/DNS/DNS_RE_640.txt'
        open(unit=10,file=m)
    elseif(Ret == 395.d0)then
        allocate(character(len=len(dirname // '/DNS/DNS_RE_395.txt')) :: m)
        m = trim(dirname) // '/DNS/DNS_RE_395.txt'
        open(unit=10,file=m)
    end if
    ! Leitura
    Do i = 1 ,  p(1)
            read(10,*) ly , Udns(1, i), Udns(2, i)
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
            if((Ret - Udns(1 , ii)) < e(i) .and. (Ret - Udns(1 , ii)) > e(i-1) )then
                ly = (u(i) - u(i-1))*((Ret - Udns(1,ii)) - e(i-1))/(e(i) - e(i-1))
                ly = u(i-1) + ly
                k1 = k1 + (ly - Udns(2,ii))**2.d0
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
            if((Ret - Udns(1 , ii)) < e(i) .and. (Ret - Udns(1 , ii)) > e(i-1) )then
                ly = (u(i) - u(i-1))*((Ret - Udns(1,ii)) - e(i-1))/(e(i) - e(i-1))
                ly = u(i-1) + ly
                k1 = k1 + abs(ly - Udns(2,ii))
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
            if((Ret - Udns(1 , ii)) < e(i) .and. (Ret - Udns(1 , ii)) > e(i-1) )then
                ly = (u(i) - u(i-1))*((Ret - Udns(1,ii)) - e(i-1))/(e(i) - e(i-1))
                ly = u(i-1) + ly
                if(abs(ly - Udns(2,ii)) > k1)then
                    k1 = abs(ly - Udns(2,ii))
                end if
            end if
        end do
    end do
    Li = k1
    return

    end subroutine Linorm