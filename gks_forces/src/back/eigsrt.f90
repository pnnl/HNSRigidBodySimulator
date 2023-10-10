subroutine eigsrt(d,v,n)
real*8 d(n)
real*8 v(n,n)
real*8 p
integer           :: n,i,j,k
do i=1,n-1
  k=i
  p=d(i)
  do j=i+1,n
    if(d(j).ge.p)then
      k=j
      p=d(j)
    endif
  enddo
  if(k.ne.i)then
    d(k)=d(i)
    d(i)=p
    do j=1,n
      p=v(j,i)
      v(j,i)=v(j,k)
      v(j,k)=p
    enddo
  endif
enddo
return
end
