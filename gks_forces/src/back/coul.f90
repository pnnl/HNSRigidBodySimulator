subroutine coul(Xc,PE,FCL)
include "setup.inc"
real*8  Xc(3,NcP,NP),PE,FCL(3,NcP,NP)
real*8  r2,r3,r4
integer ia,ib
real*8  D0,D1,D2,D3
real*8  q2,mr,rm,mm,mq,qm
integer i,j,k
real*8  r
real*8  ri
real*8  rv(3)
PE=0.0d0
FCL(:,:,:)=0.0d0
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
          r4=r2*r2
          D0=ri
          D1=(-ri)*r2
          D2=(3.0d0*r2)*r3
          D3=(-15.0d0*r3)*r4
          q2=ZQ(ia,i)*ZQ(ib,j)
          mr=dot_product(MU(:,ia,i),rv(:))
          rm=dot_product(rv(:),MU(:,ib,j))
          mm=dot_product(MU(:,ia,i),MU(:,ib,j))
          mq=mr*ZQ(ib,j)
          qm=ZQ(ia,i)*rm
          PE=PE+0.5d0*(q2*D0+mq*D1-qm*D1-mr*rm*D2-mm*D1)
          FCL(:,ia,i)=FCL(:,ia,i)-q2*rv(:)*D1
          FCL(:,ia,i)=FCL(:,ia,i)-rv(:)*mq*D2-MU(:,ia,i)*ZQ(ib,j)*D1
          FCL(:,ia,i)=FCL(:,ia,i)+rv(:)*qm*D2+ZQ(ia,i)*MU(:,ib,j)*D1
          FCL(:,ia,i)=FCL(:,ia,i)+rv(:)*mr*rm*D3+rv(:)*mm*D2
          FCL(:,ia,i)=FCL(:,ia,i)+MU(:,ia,i)*rm*D2+mr*MU(:,ib,j)*D2
        endif
      enddo
    enddo
  enddo
enddo
return
end
