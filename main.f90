!variaveis globais definidas (para que as funções que as tenham)
module prandtll
    double precision , dimension(:) , allocatable :: vPrt , T , u , e
    double precision , dimension(:,:) , allocatable :: Tdns
    integer , dimension(3) :: p
    double precision:: Ret , Re, Pr , Prt , vc , vct , dy , incre , um, L2 , L1 , Li
    integer:: N , y , Parameters
    character*100:: dirname , filename , metodo
    Logical:: criar_imagem , salvar_resultados
end module


! Programa principal

program teste

    use prandtll
    implicit none

    ! Declaram-se as variáveis do programa

    double precision :: R
    character*100 :: titulo

    ! Controle dos parametros
    y = 12
    do Parameters = 1 , y

        !   1 - Ret = 150  , Pr = 0.025
        !   2 - Ret = 150  , Pr = 0.71
        !   3 - Ret = 395  , Pr = 0.025
        !   4 - Ret = 395  , Pr = 0.71
        !   5 - Ret = 395  , Pr = 1
        !   6 - Ret = 395  , Pr = 2
        !   7 - Ret = 395  , Pr = 5
        !   8 - Ret = 395  , Pr = 7
        !   9 - Ret = 395  , Pr = 10
        !   10- Ret = 640  , Pr = 0.025
        !   11- Ret = 640  , Pr = 0.71
        !   12- Ret = 1020 , Pr = 0.71

        call DeterminaParametros()

        ! Controles numéricos

        N = 400                                                                                                    ! Número de células
        incre = 1.d-9                                                                                              ! incremento para convergência do método implícito
        R = 1.d0                                                                                                   ! Raio do canal
        dy = (R/(dble(N) - 0.5d0)) * Ret/R;                                                                        ! i_1 = dy/2 ... i_n = R
        criar_imagem = .true.                                                                                      ! cria-se arquivo de imagem?
        salvar_resultados = .true.                                                                                 ! Salva-se resultados?
        metodo = '_Genetic2temperature' !'_Prt(Ret)_Avelocity'!'_Prt(Ret)_A26' !'_classico'  !'_Prt0905_A26'       ! metodo de execução do programa.
        filename = '/results/ResultadosGeraisGenetic2temperature.txt'                                              ! nome do arquivo resultados salvo.
        titulo = 'simulação com dois vc, um diamico e um genetic'                                                  ! nome da janela.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Méta modelos a partir da referência



        ! Prt = 0.905d0                                                                                                                       ! Otimização unica para 1020


        ! Prt = 0.71                                                                                                                          !valor classico


        ! vct = 26                                                                                                                             ! Valor classico


        ! Prt = - 4.56041707672d0 * 10.d0 ** (-10.d0) * Ret**3 +  9.56902551372d0 * 10.d0 **(-7.d0) * Ret**2 &                               ! Otimizado sem a otimização de cebeci
        ! - 0.000617158206068d0 * Ret +  1.01789506426



        ! prt = ((4.52901632 * 10.d0 ** (-12.d0) ) * Ret**3 - &
        !     (5.73952059d0 * 10.d0 **(-8.d0)) * Ret**2.d0 + &
        !     (9.397008473d0 * 10.d0 ** (-5.d0) )* Ret + 0.873117480)* (pr/0.71)**(-0.008d0)                                                    ! Otimizado com a otimização de cebeci



        vc = (Ret**(log(Ret) * 0.04510621d0) * exp(5.27528132d0) ) / (Ret ** 0.60941173d0)                                                    ! Otimizado para o menor erro quanto a velocidade.


        ! vct = vc


        ! prt = (3.19791882062d-10 * Ret**3 - 1.08216023658d-06 * Ret**2 +0.00116281300928*Ret+0.449206978959)*(pr/0.71)**(-0.008d0)            ! genetic prandtl



        ! vc = exp( 0.164405721012d0 * log(Ret)**3.d0 - 2.87424334318d0 * log(Ret)**2.d0 +  16.3562873171d0 * log(Ret) - &                      ! genetic cebeci
        !     26.6310370449d0 )




        ! prt = (3.19791882062d-10 * Ret**3 - 1.08216023658d-06 * Ret**2 +0.00116281300928d0*Ret+0.449206978959d0)                          ! genetic prandtl com ajuste molecular





        ! vc = exp( 0.164405721012d0 * log(Ret)**3.d0 - 2.87424334318d0 * log(Ret)**2.d0 +  16.3562873171d0 * log(Ret) - &                      ! genetic cebeci com ajuste molecular
        !     26.6310370449d0 )



        Prt = -2.48916601371e-10 * Ret**3 +  3.60362337151e-07 * Ret**2 +  3.79213671785e-05 * Ret +  0.71234674305                                 ! genetic with 2 temperature



        vct = exp( 0.0395059904287 * log(Ret)**3 -0.758759596012 * log(Ret)**2  + 4.66369525666 * log(Ret) -5.6703426304 )                    ! genetic with 2 temperatures cebeci


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        ! Adequação aos parâmetros padrão
        call AdequaParametro()
        ! Adequação numérica final (usuário)



        print*, titulo
        print*, "Initiated algorithm"



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

        ! Amostrando resultados

        print*, "------------------------------------------------"
        print*, "Ret = " , Ret
        print*, "Pr = " , Pr
        print*, "Prt = " , prt
        print*, "A = ", vc
        print*, "A_termico = ", vct
        print*, "L1 = " , L1
        print*, "L2 = " , L2
        print*, "Li = " , Li

        end do

        ! Encerramento

        print*, " "
        print*, "End of simulations!"
        read(*,*)

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
    ! Registrando resultado
    if(salvar_resultados)then
        call EscreverArquivo()
        end if
    ! Salvar foto
    if(criar_imagem)then
        call CriarImagem()
        end if
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




! f(Y)
function f(s, i)

    use prandtll
    implicit none
    double precision, intent(in) :: s
    integer, intent(in) :: i
    double precision :: f , ff2 , Lt
    f = ( Ret/Pr - ((((Lt(s))**2.d0 )*Ret**3.d0)/prt ) * ff2(s)  )
    return

end function f



! L(Y)
function Lt(position)

    use prandtll
    implicit none
    double precision :: Lt
    double precision, intent(in) :: position
    Lt = ((0.14d0 - 0.08d0 * (position/Ret)**2.d0 - 0.06d0*(position/Ret)**4.d0 )*(1.d0 - exp((position/Ret - 1.d0)*Ret/vct)))
    return

    end function Lt





! du/dy(y)
function ff2(position)

    use prandtll
    implicit none
    double precision :: ff2 , Lt
    double precision, intent(in) :: position
    ff2 = ( -2.d0 * position * (1.d0/Ret)/( 1.d0 + sqrt(1.d0 + 4.d0 * (Lt(position))**2.d0 *Ret * position)))
    return

    end function ff2





! importando o DNS
subroutine DNSinput()

    use prandtll
    implicit none
    integer :: i, ii
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
                k1 = k1 + (ly - Tdns(2,ii))**2
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


! Determinando parametros
subroutine DeterminaParametros()

    use prandtll
    implicit none
    allocate(Tdns(12 , 2))
    Tdns(1, 1) = 150.d0
    Tdns(1, 2) = 0.025d0
    Tdns(2, 1) = 150.d0
    Tdns(2, 2) = 0.71d0
    Tdns(3, 1) = 395.d0
    Tdns(3, 2) = 0.025d0
    Tdns(4, 1) = 395.d0
    Tdns(4, 2) = 0.71d0
    Tdns(5, 1) = 395.d0
    Tdns(5, 2) = 1.d0
    Tdns(6, 1) = 395.d0
    Tdns(6, 2) = 2.d0
    Tdns(7, 1) = 395.d0
    Tdns(7, 2) = 5.d0
    Tdns(8, 1) = 395.d0
    Tdns(8, 2) = 7.d0
    Tdns(9, 1) = 395.d0
    Tdns(9, 2) = 10.d0
    Tdns(10, 1) = 640.d0
    Tdns(10, 2) = 0.025d0
    Tdns(11, 1) = 640.d0
    Tdns(11, 2) = 0.71d0
    Tdns(12, 1) = 1020.d0
    Tdns(12, 2) = 0.71d0
    Ret = Tdns(Parameters , 1)
    Pr = Tdns(Parameters , 2)
    deallocate(Tdns)
    return

    end subroutine DeterminaParametros




! Escreve Arquivo em pasta resultados
subroutine EscreverArquivo()

    use prandtll
    implicit none
    integer :: i
    double precision :: k1 , k2 , k3 , k4 , k5
    character(len = :), allocatable :: m

    if(parameters == 1)then
        allocate(character(len=len(dirname // '/results/temp.txt')) :: m)
        m = trim(dirname) // '/results/temp.txt'
        open(unit=10,file=m)
        deallocate(m)
        Write(10,*) Ret , Pr , L1 , L2 , Li
        close(10)
    elseif(parameters /= y)then
        allocate(character(len=len(dirname // filename)) :: m)
        m = trim(dirname) // trim(filename)
        open(unit=10,file=m)
        deallocate(m)
            allocate(character(len=len(dirname // '/results/temp.txt')) :: m)
            m = trim(dirname) // '/results/temp.txt'
            open(unit=11,file=m)
            deallocate(m)
            do i = 1 , parameters - 1
                read(11,*) k1 , k2 , k3 , k4 , k5
                write(10,*) k1 , k2 , k3 , k4 , k5
            end do
            close(11)
        Write(10,*) Ret , Pr , L1 , L2 , Li
        close(10)
        allocate(character(len=len(dirname // filename)) :: m)
        m = trim(dirname) // trim(filename)
        open(unit=10,file=m)
        deallocate(m)
        allocate(character(len=len(dirname // '/results/temp.txt')) :: m)
        m = trim(dirname) // '/results/temp.txt'
        open(unit=11,file=m)
        deallocate(m)
        do i = 1 , parameters
        read(10,*) k1 , k2 , k3 , k4 , k5
        write(11,*) k1 , k2 , k3 , k4 , k5
        end do
        close(11)
        close(10)
    else
        allocate(character(len=len(dirname // filename)) :: m)
        m = trim(dirname) // trim(filename)
        open(unit=10,file=m)
        deallocate(m)
            allocate(character(len=len(dirname // '/results/temp.txt')) :: m)
            m = trim(dirname) // '/results/temp.txt'
            open(unit=11,file=m)
            deallocate(m)
            do i = 1 , parameters - 1
                read(11,*) k1 , k2 , k3 , k4 , k5
                write(10,*) k1 , k2 , k3 , k4 , k5
            end do
            close(11)
        Write(10,*) Ret , Pr , L1 , L2 , Li
        close(10)
        allocate(character(len=len(dirname // filename)) :: m)
        m = trim(dirname) // trim(filename)
        open(unit=10,file=m)
        deallocate(m)
        allocate(character(len=len(dirname // '/results/temp.txt')) :: m)
        m = trim(dirname) // '/results/temp.txt'
        open(unit=11,file=m)
        deallocate(m)
        do i = 1 , parameters
        read(10,*) k1 , k2 , k3 , k4 , k5
        write(11,*) k1 , k2 , k3 , k4 , k5
        end do
        close(11)
        close(10)
        CALL SYSTEM("cd results ; rm -f temp.txt")
    end if
    return

    end subroutine EscreverArquivo



! Cria imagem temperatura.
subroutine CriarImagem()

    use prandtll
    implicit none
    character*200 :: imagename
    integer:: i

    write(imagename,FMT=187) trim(dirname) , Ret, Pr , N , trim(metodo)
    187     format( A ,'/results/graficos/image',F5.0,'_', F4.2,'_', I3 , A ,".txt")
    open(unit=10,file=imagename)
    do i = 1 , N
        write(10,*) T(i)
        end do
    close(10)
    return
    end subroutine CriarImagem