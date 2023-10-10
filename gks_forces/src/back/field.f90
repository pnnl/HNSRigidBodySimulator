subroutine field(Xc,phi,E)
include "setup.inc"
real*8  Xc(3,NcP,NP),phi(NcP,NP),E(3,NcP,NP)
real*8  r2,r3
integer ia,ib
real*8  D0,D1,D2
real*8  rm
integer i,j,k
real*8  r,ri
real*8  rv(3)
phi(:,:)=0.0d0
E(:,:,:)=0.0d0
do i=1,N
  do ia=1,Nc(i)
    do j=1,N
      do ib=1,Nc(j)
        rv(:)=Xc(:,ia,i)-Xc(:,ib,j)
        if(i.ne.j)then
          r2=dot_product(rv,rv)
          r=sqrt(r2)
          ri=1.0d0/r
          r2=1.0d0/r2
          r3=r2*ri
          D0=ri
          D1=-ri*r2
          D2=3.0d0*r2*r3
          rm=dot_product(rv(:),MU(:,ib,j))
          phi(ia,i)=phi(ia,i)+ZQ(ib,j)*D0-rm*D1
          E(:,ia,i)=E(:,ia,i)-ZQ(ib,j)*rv(:)*D1
          E(:,ia,i)=E(:,ia,i)+rv(:)*rm*D2+MU(:,ib,j)*D1
        endif
      enddo
    enddo
  enddo
enddo
return
end
