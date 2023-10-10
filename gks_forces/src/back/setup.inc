integer   DEG,NPC
parameter (DEG=4,NPC=6)
integer   NP
parameter (NP=10)
integer   NcP
parameter (NcP=1000)
integer   NmP
parameter (NmP=1000)
integer   NlP
parameter (NlP=1000)
integer   N,NAT
integer   Nc(NP),Nm(NP),Nl(NP)
real*8    rcB(3,NcP,NP),    ZQ(NcP,NP)
real*8    rlB(3,NlP,NP),    sig(NlP,NP),      eps(NlP,NP)
real*8    rmB(3,NmP,NP),    m(NmP,NP)
real*8    mTOT(NP),           In(3,NP)
real*8    del,                tolz
real*8    alpha(NcP,NP)
real*8    MU(3,NcP,NP)
real*8    pi
real*8    etaCM,etaB(3)
character*2 ch(NmP,NP)
common/setup1/N,NAT
common/setup2/Nc,Nm,Nl
common/setup3/rcB,ZQ,rlB,sig,eps,rmB,m,mTOT,In,alpha,MU
common/setup4/del,tolz,pi
common/setup5/ch
common/setup6/etaCM,etaB
