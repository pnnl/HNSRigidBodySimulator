subroutine force(id,rCM,Q,FTOT,tauL,PE)
integer id
include "setup.inc"
real*8 rCM(3,NP,0:DEG)
real*8 Q(0:3,NP,0:DEG)
real*8 FTOT(3,NP),tauL(3,NP)
real*8 PE
real*8 A(3,3)
real*8 Xc(3,NcP,NP),Xm(3,NmP,NP),Xl(3,NlP,NP)
real*8 FCL(3,NcP,NP),FLJ(3,NlP,NP)
real*8 phi(NcP,NP),E(3,NcP,NP)
real*8 z,norm,norm2
integer iter
integer ia,i,k,k1,k2
do i=1,N
  A(1,1)=Q(3,i,id)**2+Q(0,i,id)**2-Q(1,i,id)**2-Q(2,i,id)**2
  A(1,2)=2.0d0*(Q(0,i,id)*Q(1,i,id)+Q(3,i,id)*Q(2,i,id))
  A(1,3)=2.0d0*(Q(0,i,id)*Q(2,i,id)-Q(3,i,id)*Q(1,i,id))
  A(2,1)=2.0d0*(Q(0,i,id)*Q(1,i,id)-Q(3,i,id)*Q(2,i,id))
  A(2,2)=Q(3,i,id)**2-Q(0,i,id)**2+Q(1,i,id)**2-Q(2,i,id)**2
  A(2,3)=2.0d0*(Q(1,i,id)*Q(2,i,id)+Q(3,i,id)*Q(0,i,id))
  A(3,1)=2.0d0*(Q(0,i,id)*Q(2,i,id)+Q(3,i,id)*Q(1,i,id))
  A(3,2)=2.0d0*(Q(1,i,id)*Q(2,i,id)-Q(3,i,id)*Q(0,i,id))
  A(3,3)=Q(3,i,id)**2-Q(0,i,id)**2-Q(1,i,id)**2+Q(2,i,id)**2
  do ia=1,Nc(i)
    Xc(:,ia,i)=rCM(:,i,id)+matmul(transpose(A),rcB(:,ia,i))
  enddo
  do ia=1,Nl(i)
    Xl(:,ia,i)=rCM(:,i,id)+matmul(transpose(A),rlB(:,ia,i))
  enddo
  do ia=1,Nm(i)
    Xm(:,ia,i)=rCM(:,i,id)+matmul(transpose(A),rmB(:,ia,i))
  enddo
enddo
PE=0.0d0
FTOT(:,:)=0.0d0
tauL(:,:)=0.0d0
iter=0
100   continue
iter=iter+1
call field(Xc,phi,E)
norm=  0.0d0
norm2= 0.0d0
do i=1,N
  do ia=1,Nc(i)
    do k=1,3
      z= alpha(ia,i)*E(k,ia,i)
      norm=norm+(z-MU(k,ia,i))**2
      norm2=norm2+z**2
      MU(k,ia,i)=z
    enddo
  enddo
enddo
norm=sqrt(norm/real(N))
if(norm.gt.tolz)goto 100
call coul(Xc,PE,FCL)
do i=1,N
  do ia=1,Nc(i)
    do k=1,3
      PE=PE+0.5d0*MU(k,ia,i)**2/alpha(ia,i)
    enddo
  enddo
enddo
do i=1,N
  do ia=1,Nc(i)
    FTOT(:,i)=FTOT(:,i)+FCL(:,ia,i)
    tauL(1,i)=tauL(1,i)+(Xc(2,ia,i)-rCM(2,i,id))*FCL(3,ia,i)
    tauL(1,i)=tauL(1,i)-(Xc(3,ia,i)-rCM(3,i,id))*FCL(2,ia,i)
    tauL(2,i)=tauL(2,i)+(Xc(3,ia,i)-rCM(3,i,id))*FCL(1,ia,i)
    tauL(2,i)=tauL(2,i)-(Xc(1,ia,i)-rCM(1,i,id))*FCL(3,ia,i)
    tauL(3,i)=tauL(3,i)+(Xc(1,ia,i)-rCM(1,i,id))*FCL(2,ia,i)
    tauL(3,i)=tauL(3,i)-(Xc(2,ia,i)-rCM(2,i,id))*FCL(1,ia,i)
  enddo
enddo
call LJ(Xl,PE,FLJ)
do i=1,N
  do ia=1,Nl(i)
    FTOT(:,i)=FTOT(:,i)+FLJ(:,ia,i)
    tauL(1,i)=tauL(1,i)+(Xl(2,ia,i)-rCM(2,i,id))*FLJ(3,ia,i)
    tauL(1,i)=tauL(1,i)-(Xl(3,ia,i)-rCM(3,i,id))*FLJ(2,ia,i)
    tauL(2,i)=tauL(2,i)+(Xl(3,ia,i)-rCM(3,i,id))*FLJ(1,ia,i)
    tauL(2,i)=tauL(2,i)-(Xl(1,ia,i)-rCM(1,i,id))*FLJ(3,ia,i)
    tauL(3,i)=tauL(3,i)+(Xl(1,ia,i)-rCM(1,i,id))*FLJ(2,ia,i)
    tauL(3,i)=tauL(3,i)-(Xl(2,ia,i)-rCM(2,i,id))*FLJ(1,ia,i)
  enddo
enddo
return
end
