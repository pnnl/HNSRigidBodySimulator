subroutine LJ(Xl,PE,FLJ)
include "setup.inc"
real*8 Xl(3,NlP,NP)
real*8 PE
real*8 FLJ(3,NlP,NP)
real*8 r,r6
real*8 sigij,epsij
real*8 rv(3)
integer ia,ib,i,j,k
FLJ(:,:,:)=0.0d0
do i=1,N
  do ia=1,Nl(i)
    do j=1,N
      do ib=1,Nl(j)
        sigij=(sig(ia,i)+sig(ib,j))/2.0d0
        epsij=sqrt(eps(ia,i)*eps(ib,j))
        rv(:)=Xl(:,ia,i)-Xl(:,ib,j)
        if(i.ne.j)then
          r=sqrt(dot_product(rv,rv))
          r6=(sigij/r)**6
          FLJ(:,ia,i)=FLJ(:,ia,i)+24.0d0*epsij*r6*(2.0d0*r6-1.0d0)*rv(:)/r**2
          PE=PE+2.0d0*epsij*r6*(r6-1.0d0)
        endif
      enddo
    enddo
  enddo
enddo
return
end
