subroutine diag(n,A,U,lam)
integer n
real*8 A(n,n)
real*8 U(n,n)
real*8 lam(n)
real*8, allocatable :: w(:)
allocate(w(n))
U(:,:)=A(:,:)
call tred2(U,n,lam,w)
call tqli(lam,w,n,U)
call eigsrt(lam,U,n)
deallocate(w)
return
end
