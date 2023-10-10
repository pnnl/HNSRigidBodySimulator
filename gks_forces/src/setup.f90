subroutine setup
include "setup.inc"
integer NX
integer ia,i,j,k,k1,k2
real*8 rv(3)
real*8 U(3,3),UIU(3)
open(10,FILE="setup.in")
pi=4.0d0*atan(1.0d0)
i=0
100   continue
read(10,*,END=110) NX
j=i+1
read(10,*) Nc(j),tolz
do ia=1,Nc(j)
  read(10,*) rcB(:,ia,j),ZQ(ia,j),alpha(ia,j)
enddo
do ia=1,Nc(j)
  MU(:,ia,j)=0.0d0
enddo
read(10,*) Nl(j)
do ia=1,Nl(j)
  read(10,*) rlB(:,ia,j),sig(ia,j),eps(ia,j)
enddo
read(10,*) Nm(j)
do ia=1,Nm(j)
  read(10,*) ch(ia,j),rmB(:,ia,j),m(ia,j)
enddo
do k1=1,NX
  i=i+1
  Nc(i)=Nc(j)
  do ia=1,Nc(j)
    rcB(:,ia,i)=rcB(:,ia,j)
    ZQ(ia,i)=ZQ(ia,j)
    MU(:,ia,i)=MU(:,ia,j)
    alpha(ia,i)=alpha(ia,j)
  enddo
  Nl(i)=Nl(j)
  do ia=1,Nl(j)
    rlB(:,ia,i)=rlB(:,ia,j)
    sig(ia,i)=sig(ia,j)
    eps(ia,i)=eps(ia,j)
  enddo
  Nm(i)=Nm(j)
  do ia=1,Nm(j)
    rmB(:,ia,i)=rmB(:,ia,j)
    m(ia,i)=m(ia,j)
    ch(ia,i)=ch(ia,j)
  enddo
enddo
goto 100
110   continue
N=i
do i=1,N
  mTOT(i)=0.0d0
  rv(:)=0.0d0
  do ia=1,Nm(i)
    mTOT(i)=mTOT(i)+m(ia,i)
    rv(:)=rv(:)+m(ia,i)*rmB(:,ia,i)
  enddo
  do ia=1,Nm(i)
    rmB(:,ia,i)=rmB(:,ia,i)-rv(:)/mTOT(i)
  enddo
  do ia=1,Nc(i)
    rcB(:,ia,i)=rcB(:,ia,i)-rv(:)/mTOT(i)
  enddo
  do ia=1,Nl(i)
    rlB(:,ia,i)=rlB(:,ia,i)-rv(:)/mTOT(i)
  enddo
  U(:,:)=0.0d0
  do ia=1,Nm(i)
    U(1,1)=U(1,1)+m(ia,i)*(rmB(2,ia,i)**2+rmB(3,ia,i)**2)
    U(1,2)=U(1,2)-m(ia,i)*rmB(1,ia,i)*rmB(2,ia,i)
    U(1,3)=U(1,3)-m(ia,i)*rmB(1,ia,i)*rmB(3,ia,i)
    U(2,1)=U(2,1)-m(ia,i)*rmB(2,ia,i)*rmB(1,ia,i)
    U(2,2)=U(2,2)+m(ia,i)*(rmB(1,ia,i)**2+rmB(3,ia,i)**2)
    U(2,3)=U(2,3)-m(ia,i)*rmB(2,ia,i)*rmB(3,ia,i)
    U(3,1)=U(3,1)-m(ia,i)*rmB(3,ia,i)*rmB(1,ia,i)
    U(3,2)=U(3,2)-m(ia,i)*rmB(3,ia,i)*rmB(2,ia,i)
    U(3,3)=U(3,3)+m(ia,i)*(rmB(1,ia,i)**2+rmB(2,ia,i)**2)
  enddo
  call diag(3,U,U,UIU)
  In(:,i)=UIU(:)
  do ia=1,Nm(i)
    rv(:)=matmul(transpose(U),rmB(:,ia,i))
    rmB(:,ia,i)=rv(:)
  enddo
  do ia=1,Nc(i)
    rv(:)=matmul(transpose(U),rcB(:,ia,i))
    rcB(:,ia,i)=rv(:)
  enddo
  do ia=1,Nl(i)
    rv(:)=matmul(transpose(U),rlB(:,ia,i))
    rlB(:,ia,i)=rv(:)
  enddo
enddo
close(10)
NAT=0
do i=1,N
  NAT=NAT+Nm(i)
enddo
return
end
