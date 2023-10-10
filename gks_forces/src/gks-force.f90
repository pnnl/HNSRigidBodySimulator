include "setup.inc"
real*8 rCM(3,NP,0:DEG)
real*8 Q(0:3,NP,0:DEG)
real*8 A(3,3)
real*8 FTOT(3,NP),tauL(3,NP)
integer id,i,k
real*8 PE
call setup
id=0
open(10,FILE='positions.out')
do i=1,N
  read(10,*) (rCM(k,i,id),k=1,3),(Q(k,i,id),k=0,2),Q(3,i,id)
enddo
close(10)
call force(id,rCM,Q,FTOT,tauL,PE)
open(10,FILE='gks-force.txt')
write(10,'(e20.10)') PE
do i=1,N
  write(10,'(3e20.10)') (FTOT(k,i),k=1,3)
enddo
do i=1,N
  write(10,'(3e20.10)') (tauL(k,i),k=1,3)
enddo
close(10)
end
