! ---------------------------------- Pràctica 6 ------------------------------------- !
! Autor: Javier Rozalén Sarmiento
! Grup: B1B
! Data: 19/11/2019
!
! Funcionalitat: es programa el mètode d'integració de Montecarlo i es posa a prova amb 
! diverses funcions a integrar.

program practica6 
    implicit none
    double precision u,d,integral,error,nu,sigmau,nd,sigmad
    double precision, allocatable :: posi(:)
    double precision pi,L,p,g,prod,cosinus,sinus,f
    integer N,i,j,N1
    external u,d,p,g,cosinus
    common/cts/pi,L 
    pi=dacos(-1.d0)
    L=pi ! Es defineix per fer més clares les fórmules

    ! -------------------------------- Execici 1 --------------------------------------- !
    ! Càlcul del nombre de quarks de valència dins del protó (I1)
    open(11,file="P6-1920-res.dat")
    write(11,*) "#N,N_u,Sigma_u,N_d,Sigma_d"
    write(11,*) ""

    do i=150,45000,150
        call montecarlocru(1.d-14,1.d0,i,u,integral,error)
        nu=integral
        sigmau=error
        call montecarlocru(1.d-14,1.d0,i,d,integral,error)
        nd=integral
        sigmad=error
        write(11,*) i,nu,sigmau,nd,sigmad
    enddo
    print*,"El valor de nu és: ", nu, "amb error",sigmau
    print*,"El valor de nd és: ", nd, "amb error",sigmad

    ! Generació d'un milió de nombres aleatoris
    N1=1000000
    allocate(posi(N1))
    call acceptrebuig(N1,posi,-L,L,0.35d0,p)

    ! Càlcul de la integral I2 a partir del vector posi
    call write(11)
    write(11,*) "#N,I_2(N),sigma_I2(N)"

    do j=10000,N1,10000
        integral=0.d0
        error=0.d0
        do i=1,j
            integral=integral+g(posi(i))
            error=error+(g(posi(i)))**2.d0
        enddo
        integral=integral/dble(j)
        error=((dble(j))**(-0.5d0))*dsqrt((error/dble(j))-integral**2.d0)
        write(11,*) j,integral,error        
    enddo
    print*, "El valor de I2 és: ",integral, "amb error ",error

    ! -------------------------------- Execici 2 --------------------------------------- !
    ! Càlcul de la integral I3 
    call write(11)
    write(11,*) "#N,I_3(N),sigma_I3(N)"

    do j=10000,250000,10000
        integral=0.d0
        error=0.d0
        do i=1,j,3
            integral=integral+f(posi(i),posi(i+1),posi(i+2))
            error=error+(f(posi(i),posi(i+1),posi(i+2)))**2.d0
        enddo
        integral=integral/dble(j)
        error=((dble(j))**(-0.5d0))*dsqrt((error/dble(j))-integral**2.d0)
        write(11,*) j,integral,error
    enddo
    print*, "El valor de I3 és: ",integral, "amb error ",error

    close(11)
end program practica6

! Funció u --> densitat del quark up
double precision function u(x)
    ! x --> Variable independent
    implicit none
    double precision x

    u=5.109d0*((1.d0-x)**3.d0)*x**(-0.1998d0)

    return
end function u

! Funció d --> densitat del quark down
double precision function d(x)
    ! x --> Variable independent
    implicit none
    double precision x

    d=3.058d0*((1.d0-x)**4.d0)*x**(-0.197d0)

    return
end function d

! Funció p --> densitat de l'àtom ultrafred en una caixa 1-D
double precision function p(x)
    ! x --> Variable independent
    implicit none
    double precision x
    double precision pi,L
    common/cts/pi,L

    p=(1.d0/L)*(dsin((pi*(x-L))/(2.d0*L)))**2.d0

    return
end function p

! Funció g --> Funció donada a l'apartat 1) c)
double precision function g(x)
    ! x --> Variable independent
    implicit none
    double precision x
    double precision pi,L
    common/cts/pi,L

    g=(dsin((8.d0*pi*(x-L))/(2.d0*L)))**2.d0

    return
end function g

! Funció prod --> Productori estrany
double precision function f(x1,x2,x3)
    ! x --> Variable independent
    implicit none
    double precision x1,x2,x3
    double precision pi,L,cosinus,sinus,sins,prod,p
    common/cts/pi,L

    prod=abs(cosinus(x1)-cosinus(x2))*abs(cosinus(x1)-cosinus(x3))*abs(cosinus(x2)-cosinus(x3))
    sins=abs(sinus(x1)*sinus(x2)*sinus(x3))
    f=(prod*sins)**2.d0/(p(x1)*p(x2)*p(x3))

    return
end function f

! Funció cosinus --> Funció donada a l'apartat 2
double precision function cosinus(x)
    ! x --> Variable independent
    implicit none
    double precision x
    double precision pi,L
    common/cts/pi,L

    cosinus=dcos((pi*(x-L))/(2.d0*L))

    return
end function cosinus

! Funció sinus --> Funció donada a l'apartat 2
double precision function sinus(x)
    ! x --> Variable independent
    implicit none
    double precision x
    double precision pi,L
    common/cts/pi,L

    sinus=dsin((pi*(x-L))/(2.d0*L))

    return
end function sinus

! Subrutina montecarlocru --> Calcula una integral 1-D amb Montecarlo "cru"
subroutine montecarlocru(a,b,N,funci,integral,error)
    ! a,b --> Extrems de l'interval d'integració (input)
    ! N --> Nombre de números aleatoris que es volen generar (input)
    ! funci --> Funció a integrar (input)
    ! integral --> Valor de la integral (a calcular) (output)
    ! error --> Error comès en el càlcul de integral (output)
    implicit none
    double precision a,b,funci,integral,error
    integer N
    double precision x
    integer i
    integral=0.d0
    error=0.d0

    do i=1,N
        integral=integral+funci((b-a)*rand()+a)
        error=error+(funci((b-a)*rand()+a))**2.d0
    enddo

    integral=integral*(b-a)/dble(N)
    error=(dble(N))**(-0.5d0)*((error*((b-a)**2.d0)/dble(N))-integral**2.d0)**0.5d0

    return
end subroutine montecarlocru

! Subrutina acceptrebuig --> genera un vector amb ndat números distribuïts segons funci
subroutine acceptrebuig(ndat,posi,a,b,M,funci)
    ! ndat --> nombre de números aleatoris que es vol generar (input)
    ! posi --> vector que contindrà els ndat números aleatoris (output)
    ! a,b --> extrems de l'interval on està definida la nova variable aleatòria (input)
    ! M --> cota superior (input)
    ! funci --> densitat de probabilitat segons la qual està distribuida la nova var. aleat. (output)
    ! Nota --> la subrutina també calcula la variància i la desviació estàndard de la nova
    ! distribució. Es pot esborrar tranquilament aquesta part sense afectar la generació de nombres. 
    implicit none
    double precision posi(ndat),a,b,M,funci,x1,x2,p,x
    double precision valmitj,var,desvest
    integer ndat,iseed,counter,i
    counter=0 ! Variable que porta el compte dels nombres aleatoris generats

    ! Generació dels nombres aleatoris en distribució uniforme [0,1]
 1  x1=rand()
    x2=rand()
    ! Canvi de variable x1-->x, x2-->p 
    x=(b-a)*x1+a 
    p=M*x2 
    ! Ara, x és U(a,b) i p és U(0,M). Iniciem la comprovació
    if (funci(x).ge.p) then 
        posi(counter)=x
        counter=counter+1
        if (counter.le.ndat) then 
            goto 1
        endif
    else 
        goto 1
    endif

    return
end subroutine acceptrebuig

! Subrutina write --> Escriu dues línies en blanc en un arxiu
subroutine write(arxiu)
    ! arxiu --> número de l'arxiu
    implicit none
    integer arxiu

    write(arxiu,*) ""
    write(arxiu,*) ""

    return
end subroutine