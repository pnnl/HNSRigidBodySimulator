subroutine dstep(id0,id1,id2,dt,rCM,VCM,Q,LL)
include "setup.inc"
real*8 rCM(3,NP,0:DEG),VCM(3,NP,0:DEG)
real*8 Q(0:3,NP,0:DEG),LL(3,NP,0:DEG)
integer id0,id1,id2
real*8 dt  
real*8 WB(3),A(3,3),FTOT(3,NP),tauL(3,NP),tauf(3),PE
call force(id0,rCM,Q,FTOT,tauL,PE)
do i=1,N
  rCM(:,i,id2)=rCM(:,i,id1)+dt*VCM(:,i,id0)
  VCM(:,i,id2)=VCM(:,i,id1)+dt*(FTOT(:,i)-etaCM*VCM(:,i,id0))/mTOT(i)
  if(Nm(i).ge.3)then
  A(1,1)=Q(3,i,id0)**2+Q(0,i,id0)**2-Q(1,i,id0)**2-Q(2,i,id0)**2
  A(1,2)=2.0d0*(Q(0,i,id0)*Q(1,i,id0)+Q(3,i,id0)*Q(2,i,id0))
  A(1,3)=2.0d0*(Q(0,i,id0)*Q(2,i,id0)-Q(3,i,id0)*Q(1,i,id0))
  A(2,1)=2.0d0*(Q(0,i,id0)*Q(1,i,id0)-Q(3,i,id0)*Q(2,i,id0))
  A(2,2)=Q(3,i,id0)**2-Q(0,i,id0)**2+Q(1,i,id0)**2-Q(2,i,id0)**2
  A(2,3)=2.0d0*(Q(1,i,id0)*Q(2,i,id0)+Q(3,i,id0)*Q(0,i,id0))
  A(3,1)=2.0d0*(Q(0,i,id0)*Q(2,i,id0)+Q(3,i,id0)*Q(1,i,id0))
  A(3,2)=2.0d0*(Q(1,i,id0)*Q(2,i,id0)-Q(3,i,id0)*Q(0,i,id0))
  A(3,3)=Q(3,i,id0)**2-Q(0,i,id0)**2-Q(1,i,id0)**2+Q(2,i,id0)**2
  WB(:)=matmul(A,LL(:,i,id0))/In(:,i)
  Q(3,i,id2)=Q(3,i,id1)+0.5d0*dt*(-Q(0,i,id0)*WB(1)-Q(1,i,id0)*WB(2)-Q(2,i,id0)*WB(3))
  Q(0,i,id2)=Q(0,i,id1)+0.5d0*dt*(+Q(3,i,id0)*WB(1)-Q(2,i,id0)*WB(2)+Q(1,i,id0)*WB(3))
  Q(1,i,id2)=Q(1,i,id1)+0.5d0*dt*(+Q(2,i,id0)*WB(1)+Q(3,i,id0)*WB(2)-Q(0,i,id0)*WB(3))
  Q(2,i,id2)=Q(2,i,id1)+0.5d0*dt*(-Q(1,i,id0)*WB(1)+Q(0,i,id0)*WB(2)+Q(3,i,id0)*WB(3))
  tauf(:)=matmul(transpose(A(:,:)),etaB(:)*WB(:))
  LL(:,i,id2)=LL(:,i,id1)+dt*(tauL(:,i)-tauf(:))
  endif
enddo
return
end
