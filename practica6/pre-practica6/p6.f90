! ---------------------------------- Pre-pràctica 6 ------------------------------------- !
! Autor: Javier Rozalén Sarmiento
! Grup: B1B
! Data: 19/11/2019
!
! Funcionalitat: es programa el mètode d'integració de Montecarlo i es posa a prova amb 
! diverses funcions a integrar.

program pre_practica6
    implicit none
    double precision e,pi,a,b,fun1,fun2,I1,I2,err1,err2
    double precision integral,error,I1_exacte,I2_exacte
    double precision x1,x2,phi,p,int1,int2,int3,err3
    double precision sum1,sum2,sum3,g
    double precision normals(1050000),p_x(1050000)
    integer N,k,i,s
    common/cts/e,pi
    external fun1,fun2,p
    e=exp(1.d0)
    pi=acos(-1.d0)
    I1_exacte=e*(e**2.d0+pi**2.d0)**0.5d0+pi**2.d0*log(e/pi+(1.d0+(e/pi)**2.d0)**0.5d0)
    I2_exacte=pi*(1061.d0/288.d0)-(5.d0/12.d0)*pi**3.d0

    ! -------------------------------- Execici 1 --------------------------------------- !
    ! Generació d'un fitxer i càlcul de les dues primeres integrals amb Montecarlo cru
    open(11,file="P6-1920-res.dat")
    open(12,file="extra.dat") ! Fitxer extra per la figura P6-1920-fig1.png
    do k=2500,150000+1,2500
        call montecarlo(fun1,-e,e,k,integral,error)
        I1=integral
        err1=error
        call montecarlo(fun2,-pi,pi,k,integral,error)
        I2=integral
        err2=error
        write(11,*) k,I1,err1,I2,err2
        write(12,*) err1,abs(I1-I1_exacte),err2,abs(I2-I2_exacte)
    enddo
    close(12)

    ! Generació de 1050000 nombres gaussians centrats en 0 amb sigma=1
    call boxmuller(1050000,normals)

    ! Generació de 1050000 nombres distribuits segons p(x)
    call acceptrebuig(1050000,p_x,-pi,pi,0.4d0,p)

    ! Càlcul de les tres darreres integrals amb sampleig d'importància
    ! S'aprofitarà el mateix bucle per fer les 3 integrals
    int1=0.d0
    err1=0.d0
    int2=0.d0
    err2=0.d0
    int3=0.d0
    err3=0.d0
    do k=1,1050000
        ! Es calcula el terme general de cada sumatori 1 vegada i s'utilitza tant per
        ! la integral com per l'error
        sum1=(p_x(k))**2.d0
        sum2=exp(abs(p_x(k))-(p_x(k)**2.d0)/2.d0)*(1.d0+p_x(k)**2.d0)/(tan(p_x(k)))**2.d0
        sum3=exp(-normals(k)**2.d0/2.d0)*(sin(normals(k)))**4.d0*normals(k)**2.d0
        int1=int1+sum1
        err1=err1+sum1**2.d0
        int2=int2+sum2
        err2=err2+sum2**2.d0
        int3=int3+sum3
        err3=int3+sum3**2.d0
        if (mod(k,5000).eq.0) then
            ! Es multipliquen les integrals pels factors adients en cada cas; el mateix pels errors
            int1=int1*(4.d0/5.d0)*(1.d0-exp(-pi))/dble(k)
            err1=(dble(k))**(-0.5d0)*((err1/dble(k))-int1**2.d0)**0.5d0
            int2=int2/(dble(k)*((5.d0/4.d0)*1.d0/(1.d0-exp(-pi))))
            err2=(dble(k))**(-0.5d0)*((err2/dble(k))-int2**2.d0)**0.5d0
            int3=int3*sqrt(2.d0*pi)/dble(k)
            err3=(dble(k))**(-0.5d0)*((err3*2.d0*pi/dble(k))-int3**2.d0)**0.5d0
            write(11,*) ""
            write(11,*) k,int1,err1,int2,err2,int3,err3
        endif
    enddo

    ! -------------------------------- Execici 2 --------------------------------------- !
    ! Per tal d'obtenir una distribució gaussiana de 5 variables simplement s'agafaran els
    ! nombres gaussians ja generats de 5 en 5
    integral=0.d0
    error=0.d0
    sum1=0.d0
    write(11,*) ""
    write(11,*) ""
    do k=1,1050000,5
        sum1=g(normals(k),normals(k+1),normals(k+2),normals(k+3),normals(k+4))
        integral=integral+sum1
        error=error+sum1**2.d0
        if (mod(k,1500).eq.0) then
            write(11,*) k/5,integral/dble(k/5)
        endif
    enddo
    close(11)
end program pre_practica6

! Subrutina montecarlo --> Calcula una integral definida pel mètode de Montecarlo cru
subroutine montecarlo(funci,a,b,N,integral,error)
    implicit none
    double precision a,b,funci,integral,x_k,error,sumand
    integer N,k
    integral=0.d0
    error=0.d0
    do k=1,N
        x_k=rand()
        sumand=funci((b-a)*x_k+a)
        integral=integral+sumand
        error=error+sumand**2.d0
    enddo
    integral=integral*(b-a)/dble(N)
    error=(dble(N))**(-0.5d0)*((error*((b-a)**2.d0)/dble(N))-integral**2.d0)**0.5d0
    return
end subroutine montecarlo

! Subrutina acceptrebuig --> Genera nombres aleatoris distribuits segons funci
subroutine acceptrebuig(ndat,xnums,a,b,M,funci)
    implicit none
    double precision xnums(ndat),a,b,M,funci,x1,x2,p,x
    double precision valmitj,var,desvest
    integer ndat,iseed,counter,i
    counter=0

    ! Generació dels nombres aleatoris
 1   x1=rand()
    x2=rand()
    x=(b-a)*x1+a
    p=M*x2
    if (funci(x).ge.p) then 
        xnums(counter)=x
        counter=counter+1
        if (counter.le.ndat) then 
            goto 1
        endif
    else 
        goto 1
    endif
    return
end subroutine acceptrebuig

! Subrutina boxmuller --> genera un vector amb ndat números amb distribució normal centrada a 0
subroutine boxmuller(ndat,xnormal)
    ! ndat --> nombre de números aleatoris que es vol generar (input)
    ! xnormal --> vector que contindrà els ndat números aleatoris (output)
    implicit none
    double precision xnormal(ndat)
    integer ndat
    double precision pi,r,phi
    integer i
    pi=dacos(-1.d0)

    do i=1,ndat-1,2
        r=dsqrt(-2.d0*log(rand()))
        phi=2.d0*pi*rand()
        xnormal(i)=r*dcos(phi)
        xnormal(i+1)=r*dsin(phi)
    enddo

    return
end subroutine boxmuller

! Funció fun1 --> Correspon al primer integrand de l'apartat a)
double precision function fun1(x)
    implicit none
    double precision x,e,pi
    common/cts/e,pi 
    fun1=(pi**2.d0+x**2.d0)**0.5d0
    return
end function fun1

! Funció fun2 --> Correspon al segon integrand de l'apartat a)
double precision function fun2(x)
    implicit none
    double precision x
    fun2=sin(x)*(x+3.d0*(x**2.d0)*sin(x)-x**3.d0)*(cos(x))**2.d0
    return
end function fun2

! Funció p --> Distribució donada a l'apartat b)
double precision function p(x)
    implicit none
    double precision x,e,pi
    common/cts/e,pi
    p=(5.d0/4.d0)*exp(-abs(x))*((sin(x))**2.d0)/(1.d0-exp(-pi))
    return
end function p

! Funció g --> Correspon a la part de l'integrand de I6 que no és la 
! distribució normal multivariable
double precision function g(x1,x2,x3,x4,x5)
    implicit none
    double precision x1,x2,x3,x4,x5,pi,e
    common/cts/e,pi
    ! Funció g que correspon a la que es dona a l'enunciat
    g=dexp(x1*dcos(x2+x3))*((x3*x4*x5)**2.d0+x5*dsin(x5)*(dcos(x3+x4))**2.d0)
    ! Ara es divideix per la gaussiana multivariable
    g=g*dexp(-0.5d0*(x1**2.d0+x2**2.d0+x3**2.d0+x4**2.d0+x5**2.d0))*((2.d0*pi)**2.5d0)
    return
end function g