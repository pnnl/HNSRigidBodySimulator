include "setup.inc"
real*8 rCM(3,NP,0:DEG),VCM(3,NP,0:DEG)
real*8 Q(0:3,NP,0:DEG),LL(3,NP,0:DEG)
integer id
integer nt,iwrite,irestart,ithermal,ianneal,ifA
real*8  dt,Temp
real*8 WB(3)
real*8 A(3,3)
real*8 FTOT(3,NP),tauL(3,NP)
real*8 Xm(3,NmP,NP)
integer i,k,k1,k2,it
real*8 rv(3)
real*8 t0,t,KE,PE,norm
real*8 KEold
real*8 tauf(3)
real*8 dFTOT
real*8 dtauB(3)
real*8 dtauL(3)
integer idum
real*8  gasdev
real*8  kB,ANG,PS,KCAL
data kB/3.166830E-6/
data ANG/0.529177258/
data PS/2.41888439E-5/
data KCAL/627.5095/
read(*,*) t0,nt,dt,iwrite,irestart,Temp,ithermal,ianneal
read(*,*) etaCM
read(*,*) etaB(:)
t0=t0/PS
dt=dt/PS
dFTOT=sqrt(2.0d0*kB*abs(Temp)*etaCM/dt)
dtauB(:)=sqrt(2.0d0*kB*abs(Temp)*etaB(:)/dt)
call setup
id=0
open(20,FILE='state')
read(20,*) ifA
do i=1,N
  read(20,*) rCM(:,i,id)
enddo
if(ifA.eq.1)then
  do i=1,N
    do k1=1,3
      read(20,*) A(k1,:)
    enddo
    call setQ(Q,i,id,A)
  enddo
else
  do i=1,N
    read(20,*) Q(:,i,id)
  enddo
endif
if(Temp.ge.0.0d0)then
  idum=88
  rv(:)=0.0d0
  norm=0.0d0
  do i=1,N
    do k=1,3
      VCM(k,i,id)=sqrt(kB*Temp/mTOT(i))*gasdev(idum)
    enddo
    rv(:)=rv(:)+mTOT(i)*VCM(:,i,id)
    norm=norm+mTOT(i)
  enddo
  do i=1,N
    VCM(:,i,id)=VCM(:,i,id)-rv(:)/norm
  enddo
  do i=1,N
    if(Nm(i).ge.3)then
    do k=1,3
      WB(k)=sqrt(kB*Temp/In(k,i))*gasdev(idum)
    enddo
    A(1,1)=Q(0,i,id)**2+Q(1,i,id)**2-Q(2,i,id)**2-Q(3,i,id)**2
    A(1,2)=2.0d0*(Q(1,i,id)*Q(2,i,id)+Q(0,i,id)*Q(3,i,id))
    A(1,3)=2.0d0*(Q(1,i,id)*Q(3,i,id)-Q(0,i,id)*Q(2,i,id))
    A(2,1)=2.0d0*(Q(1,i,id)*Q(2,i,id)-Q(0,i,id)*Q(3,i,id))
    A(2,2)=Q(0,i,id)**2-Q(1,i,id)**2+Q(2,i,id)**2-Q(3,i,id)**2
    A(2,3)=2.0d0*(Q(2,i,id)*Q(3,i,id)+Q(0,i,id)*Q(1,i,id))
    A(3,1)=2.0d0*(Q(1,i,id)*Q(3,i,id)+Q(0,i,id)*Q(2,i,id))
    A(3,2)=2.0d0*(Q(2,i,id)*Q(3,i,id)-Q(0,i,id)*Q(1,i,id))
    A(3,3)=Q(0,i,id)**2-Q(1,i,id)**2-Q(2,i,id)**2+Q(3,i,id)**2
    LL(:,i,id)=matmul(transpose(A(:,:)),In(:,i)*WB(:))
    else
    LL(:,i,id)=0.0d0
    endif
  enddo
else
  do i=1,N
    read(20,*) VCM(:,i,id)
  enddo
  do i=1,N
    read(20,*) LL(:,i,id)
  enddo
endif
close(20)
open(10,FILE="E.out")
open(11,FILE="traj.xyz")
KE=0.0d0
do it=0,nt-1
  t=t0+real(it)*dt

  call dstep(0,0,1, 0.5d0*dt,rCM,VCM,Q,LL)
  call dstep(1,0,2, 0.5d0*dt,rCM,VCM,Q,LL)
  call dstep(2,0,3, 1.0d0*dt,rCM,VCM,Q,LL)
  call dstep(3,0,4,-0.5d0*dt,rCM,VCM,Q,LL)

  do i=1,N
    rCM(:,i,0)=(rCM(:,i,1)+2.0d0*rCM(:,i,2)+rCM(:,i,3)-rCM(:,i,4))/3.0d0
    VCM(:,i,0)=(VCM(:,i,1)+2.0d0*VCM(:,i,2)+VCM(:,i,3)-VCM(:,i,4))/3.0d0
    LL(:,i,0)= ( LL(:,i,1)+2.0d0* LL(:,i,2)+ LL(:,i,3)- LL(:,i,4))/3.0d0
    Q(:,i,0)=  (  Q(:,i,1)+2.0d0*  Q(:,i,2)+  Q(:,i,3)-  Q(:,i,4))/3.0d0
  enddo

  do i=1,N
    do k=1,3
      rv(k)=dtauB(k)*gasdev(idum)
    enddo
    dtauL(:)=matmul(transpose(A(:,:)),rv(:))
    LL(:,i,0)=LL(:,i,0)+dt*dtauL(:)
    do k=1,3
      VCM(k,i,0)=VCM(k,i,0)+dt*dFTOT*gasdev(idum)/mTOT(i)
    enddo
  enddo
  id=0
  do i=1,N
    if(Nm(i).ge.3)then
    norm=sqrt(dot_product(Q(:,i,id),Q(:,i,id)))
    Q(:,i,id)=Q(:,i,id)/norm
    endif
  enddo
  KEold=KE
  call force(id,rCM,Q,FTOT,tauL,PE)
  KE=0.0d0
  do i=1,N
    KE=KE+0.5*mTOT(i)*dot_product(VCM(:,i,id),VCM(:,i,id))
    if(Nm(i).ge.3)then
    A(1,1)=Q(0,i,id)**2+Q(1,i,id)**2-Q(2,i,id)**2-Q(3,i,id)**2
    A(1,2)=2.0d0*(Q(1,i,id)*Q(2,i,id)+Q(0,i,id)*Q(3,i,id))
    A(1,3)=2.0d0*(Q(1,i,id)*Q(3,i,id)-Q(0,i,id)*Q(2,i,id))
    A(2,1)=2.0d0*(Q(1,i,id)*Q(2,i,id)-Q(0,i,id)*Q(3,i,id))
    A(2,2)=Q(0,i,id)**2-Q(1,i,id)**2+Q(2,i,id)**2-Q(3,i,id)**2
    A(2,3)=2.0d0*(Q(2,i,id)*Q(3,i,id)+Q(0,i,id)*Q(1,i,id))
    A(3,1)=2.0d0*(Q(1,i,id)*Q(3,i,id)+Q(0,i,id)*Q(2,i,id))
    A(3,2)=2.0d0*(Q(2,i,id)*Q(3,i,id)-Q(0,i,id)*Q(1,i,id))
    A(3,3)=Q(0,i,id)**2-Q(1,i,id)**2-Q(2,i,id)**2+Q(3,i,id)**2
    WB(:)=matmul(A,LL(:,i,id))/In(:,i)
    KE=KE+0.5d0*dot_product(WB(:),In(:,i)*WB(:))
    endif
  enddo
  if((ianneal.eq.1).and.(KE.lt.KEold))then
    KE=0.0d0
    VCM(:,:,id)=0.0d0
    LL(:,:,id)=0.0d0
  endif
  if(ithermal.gt.0)then
  if(mod(it,ithermal).eq.0)then
  rv(:)=0.0d0
  norm=0.0d0
  do i=1,N
    do k=1,3
      VCM(k,i,id)=sqrt(kB*Temp/mTOT(i))*gasdev(idum)
    enddo
    rv(:)=rv(:)+mTOT(i)*VCM(:,i,id)
    norm=norm+mTOT(i)
  enddo
  do i=1,N
    VCM(:,i,id)=VCM(:,i,id)-rv(:)/norm
  enddo
  do i=1,N
    if(Nm(i).ge.3)then
    do k=1,3
      WB(k)=sqrt(kB*Temp/In(k,i))*gasdev(idum)
    enddo
    A(1,1)=Q(0,i,id)**2+Q(1,i,id)**2-Q(2,i,id)**2-Q(3,i,id)**2
    A(1,2)=2.0d0*(Q(1,i,id)*Q(2,i,id)+Q(0,i,id)*Q(3,i,id))
    A(1,3)=2.0d0*(Q(1,i,id)*Q(3,i,id)-Q(0,i,id)*Q(2,i,id))
    A(2,1)=2.0d0*(Q(1,i,id)*Q(2,i,id)-Q(0,i,id)*Q(3,i,id))
    A(2,2)=Q(0,i,id)**2-Q(1,i,id)**2+Q(2,i,id)**2-Q(3,i,id)**2
    A(2,3)=2.0d0*(Q(2,i,id)*Q(3,i,id)+Q(0,i,id)*Q(1,i,id))
    A(3,1)=2.0d0*(Q(1,i,id)*Q(3,i,id)+Q(0,i,id)*Q(2,i,id))
    A(3,2)=2.0d0*(Q(2,i,id)*Q(3,i,id)-Q(0,i,id)*Q(1,i,id))
    A(3,3)=Q(0,i,id)**2-Q(1,i,id)**2-Q(2,i,id)**2+Q(3,i,id)**2
    LL(:,i,id)=matmul(transpose(A(:,:)),In(:,i)*WB(:))
    else
    LL(:,i,id)=0.0d0
    endif
  enddo
  endif
  endif
  if(mod(it,iwrite).eq.0)then
  do i=1,N
    A(1,1)=Q(0,i,id)**2+Q(1,i,id)**2-Q(2,i,id)**2-Q(3,i,id)**2
    A(1,2)=2.0d0*(Q(1,i,id)*Q(2,i,id)+Q(0,i,id)*Q(3,i,id))
    A(1,3)=2.0d0*(Q(1,i,id)*Q(3,i,id)-Q(0,i,id)*Q(2,i,id))
    A(2,1)=2.0d0*(Q(1,i,id)*Q(2,i,id)-Q(0,i,id)*Q(3,i,id))
    A(2,2)=Q(0,i,id)**2-Q(1,i,id)**2+Q(2,i,id)**2-Q(3,i,id)**2
    A(2,3)=2.0d0*(Q(2,i,id)*Q(3,i,id)+Q(0,i,id)*Q(1,i,id))
    A(3,1)=2.0d0*(Q(1,i,id)*Q(3,i,id)+Q(0,i,id)*Q(2,i,id))
    A(3,2)=2.0d0*(Q(2,i,id)*Q(3,i,id)-Q(0,i,id)*Q(1,i,id))
    A(3,3)=Q(0,i,id)**2-Q(1,i,id)**2-Q(2,i,id)**2+Q(3,i,id)**2
    do ia=1,Nm(i)
      Xm(:,ia,i)=rCM(:,i,id)+matmul(transpose(A(:,:)),rmB(:,ia,i))
    enddo
  enddo
  write(10,'(6e16.8)') t*PS,KE*KCAL,PE*KCAL,(KE+PE)*KCAL
  write(11,'(i5)') NAT
  write(11,'(i5,4e14.5)') it,t*PS,KE*KCAL,PE*KCAL,(KE+PE)*KCAL
  do i=1,N
    do ia=1,Nm(i)
      write(11,'(a2,3f14.5)') ch(ia,i),Xm(:,ia,i)
    enddo
  enddo
  endif
enddo
close(10)
close(11)
open(10,FILE='state0')
id=0
write(10,'(i5)') 0
do i=1,N
  write(10,'(3f20.10)') rCM(:,i,id)
enddo
do i=1,N
  write(10,'(4f20.10)') Q(:,i,id)
enddo
close(10)
open(10,FILE='state1')
id=0
write(10,'(i5)') 1
do i=1,N
  write(10,'(3f20.10)') rCM(:,i,id)
enddo
do i=1,N
  do k1=1,3
    write(10,'(3f20.10)') A(k1,:)
  enddo
enddo
close(10)
end
