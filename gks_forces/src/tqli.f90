subroutine tqli(d,e,n,z)
real*8 d(n),e(n),z(n,n)
real*8 dd,g,r,s,c,p,f,b
integer n,i,l,m,iter,k
if (n.gt.1) then
  do i=2,n
    e(i-1)=e(i)
  enddo
  e(n)=0.0d0
  do l=1,n
    iter=0
1   do m=l,n-1
      dd=abs(d(m))+abs(d(m+1))
      if (abs(e(m))+dd.eq.dd) goto 2
    enddo
    m=n
2   if(m.ne.l)then
      if(iter.eq.30)then
        print*, 'tqli: too many iterations'
        stop
      endif
      iter=iter+1
      g=(d(l+1)-d(l))/(2.0d0*e(l))
      r=sqrt(g**2+1.0d0)
      g=d(m)-d(l)+e(l)/(g+sign(r,g))
      s=1.0d0
      c=1.0d0
      p=0.0d0
      do i=m-1,l,-1
        f=s*e(i)
        b=c*e(i)
        if(abs(f).ge.abs(g))then
          c=g/f
          r=sqrt(c**2+1.0d0)
          e(i+1)=f*r
          s=1.0d0/r
          c=c*s
        else
          s=f/g
          r=sqrt(s**2+1.0d0)
          e(i+1)=g*r
          c=1.0d0/r  
          s=s*c
        endif
        g=d(i+1)-p
        r=(d(i)-g)*s+2.0d0*c*b
        p=s*r
        d(i+1)=g+p
        g=c*r-b
        do k=1,n
          f=z(k,i+1)
          z(k,i+1)=s*z(k,i)+c*f
          z(k,i)=c*z(k,i)-s*f
        enddo
      enddo
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.0d0
      goto 1
    endif
  enddo
endif
return
end
