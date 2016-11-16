program main

use, intrinsic :: iso_fortran_env, only: iostat_end

use datam
implicit none

integer nh

!Decode ole file, extract and open Excel stream


character*4 :: ext, tout
character*50 :: filee, outfile, infile
character*50 :: rep

real*8 :: a,b,c, bet,bt, sx,sy, sxy,sx2,at,eta, L10
real*8 :: BET0, ratio, ST,ST2,ST3,SS, SS2

real*8, dimension(50) ::  trank,x,y
character*1, dimension(50) ::  tstatus
character*3 d


integer i,j,n, prank, invrank, jk

integer, dimension (50) :: AJrank

real*8, dimension(2) :: fvec, xh, xopt
real*8 :: vt(2,1), vp(1,2)
external :: testf
integer :: info
real*8 :: tol, F, sigb, sige, R, R2, P, Q
real*8 :: S1, S2, WETA, WBETA,F11, F12, F21, F22, T, V
real*8, dimension(2,2) :: matif, fisher
real*8 :: dteta, dtbeta
real*8 :: v10i(1,2), v10
real*8 :: maxb, minb, maxe, mine, max10, min10, KK


integer io_status

!DEC$ IF DEFINED(_WIN32)
rep= 'C:\Calcul\Etudes\Open\weibull\'
!DEC$ ELSEIF DEFINED(__linux)
rep = '/home/bruno/Etude/weibull/'
!DEC$ ELSE
... Oops
!DEC$ ENDIF


ext = '.dat'
tout = '.txt'


!filee = 'DPE5035' ;
filee = 'Book1';



infile =  TRIM(rep)//TRIM(filee)//ext
outfile = filee//tout ;

write(*,*) 'infile = ', infile
open(5, FILE=infile);

write(*,*) infile
write(*,*) outfile

! total number of samples
n = 0;

!Read Data
do i = 1,50
read(5,*, iostat=io_status) a,b,c,d
if(io_status == iostat_end) exit
tTime(i) =a
tnumber(i)=b
tstatus(i)=d
n = n + 1
enddo

! close input file
close(5)

!open the  ouput stream
open(6,FILE=outfile)

! find all case where status is Failed
nk = 0
do i = 1,n
if (tstatus(i)=='D') then
 nk = nk+1;
 k(nk) = dble(i)
 endif
enddo

! find all case where status is Suspended
ns = 0
do i = 1,n
if (tstatus(i)=='S') then
 ns = ns+1;
 s(ns)= dble(i)
 endif
enddo

!rank correction if one S is founded

prank = 0 ;
jk = 0
if ( ns > 0) then
do i = 1,n
    if (tstatus(i)== 'D')  then
        jk = jk + 1 ;
        invRank = n-i+1;
        AJrank(jk) = ((invRank*prank)+ n+1.)/ (invRank+1.)
        k(jk) = AJrank(jk)
    endif
enddo

end if

continue

if (nk>1) then

! Rank regression by Bernard's approximation
trank = 0.0D+00
Do i = 1,nk
   trank(i)= (k(i)-0.3)/(real(n)+0.4) ;
   !write(*,*) i, trank(i)
enddo

!regression on Y
do i = 1,nk
j = k(i)
x(i)=ALOG(real(tTime(j)));
y(i) = ALOG(ALOG( 1./(1.-real(trank(i)) ) ));
enddo

!tab are set to 0
sx  = 0.0
sy  = 0.0
sxy = 0.0
sx2 = 0.0


! calculation of summations
do i=1,nk
sx = sx + x(i);
sy = sy + y(i);
sxy = sxy+ x(i)*y(i);
sx2 = sx2 + x(i)*x(i);
enddo



! our slope result
bt = (sxy -sx*sy/nk)/(sx2-sx*sx/nk)
at = sy/nk-bt*sx/nk

!our distribution results
eta = exp(-at/bt)
bet = bt

C = (0.1054)**(1.0/bet)
L10 = C*eta

write(*,100) "RRY", "beta=",bet," eta=",eta, "L10=", L10
write(6,100) "RRY", "beta=",bet," eta=",eta, "L10=", L10

100 format(A10, A10,F7.3, A10, F7.3, A10, F7.2)

! regression on X
do i = 1,nk
y(i)=log(real(tTime(k(i))))
x(i) = log(log( 1/(1-real(trank(i)) ) ));
enddo

! tab are set to 0
sx=0.0 ; sy = 0.0 ; sxy=0.0 ; sx2=0.0 ;


! calculation of summations
do i=1,nk
sx = sx + x(i);
sy = sy + y(i);
sxy = sxy+ x(i)*y(i);
sx2 = sx2 + x(i)*x(i);
enddo

!our slope result
bt = (sxy -sx*sy/nk)/(sx2-sx*sx/nk);
at = sy/nk-bt*sx/nk;



!our distribution results
bet = 1/bt ;
eta = exp(at);

C = (0.1054)**(1.0/bet);
L10 = C*eta

write(*,100) "RRX", "beta=",bet," eta=",eta, "L10=", L10
write(6,100) "RRX", "beta=",bet," eta=",eta, "L10=", L10

endif



!---------------------------------------------------------Approximate MLE-------------------------------------

bet0 = 3.;
ratio = 20;


do while (ratio > 1)
st = 0.0 ; st2 = 0. ; st3=0.; ss= 0. ;  ss2 = 0. ;

do i = 1,nk
j =k(i);
st = st +  tTime(j)**bet * log(tTime(j));
st2 = st2 + tTime(j)**bet;
st3 = st3 + log(tTime(j));

enddo



do i = 1,ns
j =s(i);
ss = ss + (n-nk)*tTime(j)**bet * log(real(tTime(j)))
ss2 = st2 +(n-nk)* tTime(j)**bet
enddo

bet =  st/st2 -1.0/nk*st3;
bet = 1.0/bet;

ratio = abs((bet -bet0)/bet0)*100;
bet0 = bet;
enddo

eta = (1.0/nk)*st2;

if (ns > 0) then
    eta = eta + 1.0/ns*ss2;
endif
eta = eta**(1.0/bet);

C = (0.1054)**(1.0/bet);
L10 = C*eta


!write(*,*) "MLE ", "beta=",bet," eta=",eta, "L10=", L10
!write(6,*) "MLE ", "beta=",bet," eta=",eta, "L10=", L10

! --------------------------------------Approximate mle is used to initialize the simplex solver

nh = 2
xh(1) = bet
xh(2) = eta
xopt = 0.0

call MINIMIZE_2D(xh,testf,xopt,F,Info)

C = (0.1054)**(1.0/xopt(1));
L10 = C*xopt(2)

write(*,100) " MLE", "beta=",xopt(1)," eta=",xopt(2), "L10=", L10
write(6,100) " MLE", "beta=",xopt(1)," eta=",xopt(2), "L10=", L10


Wbeta = xopt(1)
Weta = xopt(2)


! ------------------------------------Fisher Matrix ----------------------------------------------


!  ------------------------------------22 term
S1= 0. ; S2 = 0.
if (nk>0) then
do i =1,nk
    S1 = S1 + (Wbeta+1)*(tTime(k(i))/Weta)**Wbeta ;
enddo
end if

if (ns>0) then
do i =1,ns
    S2 = S2 + (Wbeta+1)*(tTime(s(i))/Weta)**Wbeta ;
enddo
end if
f22 = -nk*Wbeta/Weta**2 + Wbeta/Weta**2 * S1 + Wbeta/Weta**2 * S2 ;


!---------------------------------- 21 & 12 terms
S1= 0. ; S2 = 0.
if (nk>0) then
do i =1,nk
    T = tTime(k(i));
    V = 1.0 + Wbeta*log(T/Weta) ;
    S1 = S1 +(T/Weta)**Wbeta*V;
enddo
end if

if (ns>0) then
do i =1,ns
    T = tTime(s(i));
    V = 1.0 + Wbeta*log(T/Weta) ;
    S2 = S2 + (T/Weta)**Wbeta * V;
enddo
end if
f21 = DBLE(nk)/Weta - S1/Weta - S2 /Weta;
f12 = f21;

!------------------------------------ 11 term

S1= 0. ; S2 = 0.
if (nk>0) then
do i =1,nk
    T= tTime(k(i));
    S1 = S1 + (T/Weta)**Wbeta*log(T/Weta)*log(T/Weta);
enddo
end if

if (ns>1) then

do i =1,ns
    T= tTime(s(i));
    S2 = S2 + (T/Weta)**Wbeta*log(T/Weta)*log(T/Weta);
enddo
end if
f11 = DBLE(nk)/Wbeta**2 + S1 + S2;

matif(1,1) = F11
matif(1,2) = F12
matif(2,1) = F21
matif(2,2) = F22

call invmat(matif, fisher, 2)
!
write(*,*)
write(*,*) 'Matrice de Fisher :'
write(*,*)
write(*,'(2F8.3)') fisher(1,1),fisher(1,2)
write(*,'(2F8.3)') fisher(1,2),fisher(2,2)
!
write(6,*)
write(6,*) 'Matrice de Fisher :'
write(6,*)
write(6,'(2F8.3)') fisher(1,1),fisher(1,2)
write(6,'(2F8.3)') fisher(1,2),fisher(2,2)


!----------------------loi normale pour beta  et eta

sigb = dsqrt(fisher(1,1));
sige = dsqrt(fisher(2,2));

! Niveau de confiance 90% deux cotés
R = 0.9/2.;
R2 = 0.9;


P = R ; Q = 1-P ;

dtEta = (-log(R2))**(1.0/Wbeta) ;
dtBeta = -Weta/Wbeta**2*(-log(R2))**(1.0/Wbeta) * log(-log(R2)) ;

vt(1,1)  = dtBeta
vt(2,1) = dtEta

vp(1,1) = dtBeta
vp(1,2) =  dtEta

v10i = matmul(vp, fisher)
v10 = v10i(1,1)*vt(1,1) + v10i(1,2)*vt(2,1)


kk= 1.6448536


maxb = Wbeta * dexp(kk*sigb/Wbeta) ;
minb = Wbeta / dexp(kk*sigb/Wbeta) ;

maxe = Weta * dexp(kk*sige/Weta) ;
mine = Weta / dexp(kk*sige/Weta) ;

max10 = L10 * dexp(kk*dsqrt(v10)/L10) ;
min10 = L10 / dexp(kk*dsqrt(v10)/L10) ;



write(*,'(/,A20, F5.2)') 'Min et Max  de Beta et Eta confiance à :', 2*P
write(*,'(2F7.2)') minb,maxb
write(*,'(2F7.2)') mine,maxe
!
write(*,'(/,A20, F5.2)') 'Min et Max  de L10 à :', 2*P
write(*,'(2F7.2)') min10,max10
!
write(6,'(/,A20, F5.2)') 'Min et Max  de Beta et Eta confiance à :', 2*P
write(6,'(2F7.2)') minb,maxb
write(6,'(2F7.2)') mine,maxe
!
write(6,'(/,A20, F5.2)') 'Min et Max  de L10 à :', 2*P
write(6,'(2F7.2)') min10,max10





 end program




subroutine testf(x,f)

    use datam

       real*8, dimension(2) :: x
       real*8  :: f


       real*8 :: a,b, sumi, sumj, a0, a1,b1
       integer i,j

        a = x(1);
        b = x(2);

        Sumi = 0.0;
        Sumj = 0.0;
        do i=1,nk
            j= k(i);
            a0 = tNumber(j);
            a1=  (tTime(j)/b)**(a-1);
            b1 = exp(-(tTime(j)/b)**a);
            Sumi = Sumi + a0 * log(a/b*a1* b1);
        end do

        do i=1,ns
            j = s(i);
            Sumj = Sumj  + TNumber(j)*(tTime(j)/b)**a;
        end do

        f =-(Sumi-Sumj)


   end subroutine testf


