program main
    implicit none
    integer,parameter::N=5
    integer method
    real*4,dimension(N)::Fe,Fw,De,Dw,AE,AW,AP,p,y,x,f,u,a
    real*4,parameter::r=1.0,g=0.1,lt=1.0,pA=1.0,pB=0.0
    real*4 V
    real*4 dx,L(2:N),d(1:N-1),c(1:N-1),e(2:N),F0,D0
    integer::i

	write(*,*) '输入选择:1. 中心差分；2. 上风格式；3. 混合格式'
    read(*,*) method
    write(*,*) '速度 V='
    read(*,*) V
    dx=lt/N
    F0=r*V
    D0=g/dx
    
    do i=1,N
        Fe(i)=F0
        De(i)=D0
        Fw(i)=F0
        Dw(i)=D0
    end do
    De(N)=2*D0
    Dw(1)=2*D0
    !write(*,*) (DW(i),i=1,N)


    select case(method)
    case(1)
        do i=1,N
            AW(i)=Dw(i)+0.5*Fw(i)
            AE(i)=De(i)-0.5*Fe(i)
        end do
        AW(1)=Dw(1)+Fw(1)
        AE(N)=De(N)-Fe(N)
        !write(*,*) (AW(i),i=1,N)
    case(2)
        do i=1,N
            AW(i)=Dw(i)+max(Fw(i),0.0)
            AE(i)=De(i)+max(-Fe(i),0.0)
        end do
        AW(1)=Dw(1)+Fw(1)
        AE(N)=De(N)
    case(3)
        do i=1,N
            AW(i)=max(Fw(i),Dw(i)+0.5*Fw(i),0.0)
            AE(i)=max(-Fe(i),De(i)-0.5*Fe(i),0.0)
        end do
        AW(1)=Dw(1)+Fw(1)
        AE(N)=De(N)
    end select
    do i=1,N
        AP(i)=AE(i)+AW(i)+(Fe(i)-Fw(i))
    end do
    write(*,200) AW,AE,AP
    200 format('AW=',5(/F10.5),/,'AE=',5(/F10.5),/,'AP=',5(/F10.5))
    do i=1,N
        p(i)=1-(exp(r*V*(i-0.5)*dx/g)-1)/(exp(r*V/g)-1)
    end do

    do i=1,N
        a(i)=AP(i)
        f(i)=0
    end do
    f(1)=AW(1)*pA
    f(N)=AE(N)*pB
    write(*,210) f
    210 format('f=',5(/F10.5))
    do i=1,N-1
        c(i)=-AE(i)
    end do
    do i=2,N
        e(i)=-AW(i)
    end do

    do i=1,N-1
        d(i)=c(i)
    end do
    u(1)=a(1)
    do i=2,N
        L(i)=e(i)/u(i-1)
        u(i)=a(i)-L(i)*c(i-1)
    end do
  
    y(1)=f(1)
    do i=2,N
        y(i)=f(i)-L(i)*y(i-1)
    end do
  
    x(n)=y(n)/u(n)
    do i=n-1,1,-1
        x(i)=(y(i)-c(i)*x(i+1))/u(i)
    end do

    write(*,100) x,p
    100 format(/,'数值解：',/,5(/F10.6),/,'解析解：',/,5(/F10.6))

end program


