module kernel
integer, constant :: imin,imax,jmin,jmax,kmin,kmax
real, constant :: RE,GA,dx,dy,dz,dt
contains

attributes(global) subroutine BOUNDARY(U1,U2,V1,V2,W1,W2)
implicit none

real,dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),intent(inout) ::&
	U1,U2,V1,V2,W1,W2
integer :: i,j,k
i=blockDim%x*(blockIdx%x-1)+threadIdx%x
j=blockDim%y*(blockIdx%y-1)+threadIdx%y
k=blockDim%z*(blockIdx%z-1)+threadIdx%z

!x=imax No-Slip
IF(i.eq.imin) then
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN 
U1(imax+1,j,k)=-U1(imax,j,k)
U2(imax+1,j,k)=-U2(imax,j,k)
V1(imax+1,j,k)=1.0
V2(imax+1,j,k)=1.0
W1(imax+1,j,k)=-W1(imax,j,k)
W2(imax+1,j,k)=-W2(imax,j,k)

!x=imin No-Splip
U1(imin-1,j,k)=0.0
U2(imin-1,j,k)=0.0
V1(imin-1,j,k)=-V1(imin,j,k)
V2(imin-1,j,k)=-V2(imin,j,k)
W1(imin-1,j,k)=0.0
W2(imin-1,j,k)=0.0
END IF
END IF
END IF

!y=jmin No-Sip
IF(j.eq.jmin) then
IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((k.GE.kmax).AND.(k.LE.kmax)) THEN
U1(i,jmin-1,k)=-U1(i,jmin,k)
U2(i,jmin-1,k)=-U2(i,jmin,k)
V1(i,jmin-1,k)=0.0
V2(i,jmin-1,k)=0.0
W1(i,jmin-1,k)=0.0
W2(i,jmin-1,k)=0.0

!y=jmax No-Slip
U1(i,jmax+1,k)=0.0
U2(i,jmax+1,k)=0.0
V1(i,jmax+1,k)=-V1(i,jmax,k)
V2(i,jmax+1,k)=-V2(i,jmax,k)
W1(i,jmax+1,k)=-W1(i,jmax,k)
W2(i,jmax+1,k)=-W2(i,jmax,k)

END IF
END IF
END IF

IF(k.eq.kmin) then
IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
!z=kmin
U1(i,j,kmin-1)=-U1(i,j,kmin)
U2(i,j,kmin-1)=-U2(i,j,kmin)
V1(i,j,kmin-1)=-V1(i,j,kmin)
V2(i,j,kmin-1)=-V2(i,j,kmin)
W1(i,j,kmin-1)=0.0
W2(i,j,kmin-1)=0.0

!z=kmax
U1(i,j,kmax+1)=0.0
U2(i,j,kmax+1)=0.0
V1(i,j,kmax+1)=0.0
V2(i,j,kmax+1)=0.0
W1(i,j,kmax+1)=-W1(i,j,kmax)
W2(i,j,kmax+1)=-W2(i,j,kmax)
END IF
END IF
END IF

return
end subroutine BOUNDARY

attributes(global) subroutine CALCF(F,U1,V1,W1)
implicit none

real,dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),intent(inout) :: U1,V1,W1,F
real :: d2udx2,d2udy2,d2udz2,du2dx,duvdy,duwdz,FLX

integer :: i,j,k

i=blockDim%x*(blockIdx%x-1)+threadIdx%x
j=blockDim%y*(blockIdx%y-1)+threadIdx%y
k=blockDim%z*(blockIdx%z-1)+threadIdx%z

IF((i.GE.imin).AND.(i.LE.(imax-1))) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN

d2udx2=(u1(i+1,j,k) - 2.0*u1(i,j,k) + u1(i-1,j,k))/(dx*dx)
d2udy2=(u1(i,j+1,k) - 2.0*u1(i,j,k) + u1(i,j-1,k))/(dy*dy)
d2udz2=(u1(i,j,k+1) - 2.0*u1(i,j,k) + u1(i,j,k-1))/(dz*dz)

du2dx=(0.25/dx)*( (u1(i,j,k)+u1(i+1,j,k))**2 &
		- (u1(i-1,j,k)+u1(i,j,k))**2) &
		+(0.25*Ga/dx)* &
		( ABS(u1(i,j,k)+u1(i+1,j,k))*(u1(i,j,k)-u1(i+1,j,k))&
		 - ABS(u1(i-1,j,k)+u1(i,j,k))*(u1(i-1,j,k)-u1(i,j,k)) )

duvdy=(0.25/dy)*( (v1(i,j,k)+v1(i+1,j,k))*(u1(i,j,k)+u1(i,j+1,k)) &
		- (v1(i,j-1,k)+v1(i+1,j-1,k))*(u1(i,j-1,k)+u1(i,j,k)) ) &
		+(0.25*Ga/dy)* &
		( ABS(v1(i,j,k)+v1(i+1,j,k))*(u1(i,j,k)-u1(i,j+1,k)) &
		- ABS(v1(i,j-1,k)+v1(i+1,j-1,k))*(u1(i,j-1,k)-u1(i,j,k)) )

duwdz=(0.25/dz)*( (w1(i,j,k)+w1(i+1,j,k))*(u1(i,j,k)+u1(i,j,k+1)) &
		- (w1(i,j,k-1)+w1(i+1,j,k-1))*(u1(i,j,k-1)+u1(i,j,k)) ) &
		+(0.25*Ga/dz)* &
		( ABS(w1(i,j,k)+w1(i+1,j,k))*(u1(i,j,k)-u1(i,j,k+1)) &
		- ABS(w1(i,j,k-1)+w1(i+1,j,k-1))*(u1(i,j,k-1)-u1(i,j,k)) )

FLX=0.0

F(i,j,k)=u1(i,j,k)+dt*( (1./RE)*(d2udx2 + d2udy2 + d2udz2)&
		-du2dx -duvdy -duwdz + FLX )


END IF
END IF
END IF
return
end subroutine


attributes(global) subroutine CALCG(G,U1,V1,W1)
implicit none

real,dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),intent(inout) :: U1,V1,W1,G
real :: d2vdx2,d2vdy2,d2vdz2,dvudx,dv2dy,dvwdz,FLY

integer :: i,j,k

i=blockDim%x*(blockIdx%x-1)+threadIdx%x
j=blockDim%y*(blockIdx%y-1)+threadIdx%y
k=blockDim%z*(blockIdx%z-1)+threadIdx%z

IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.(jmax-1))) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN

d2vdx2=(v1(i+1,j,k) - 2.0*v1(i,j,k) + v1(i-1,j,k))/(dx*dx)
d2vdy2=(v1(i,j+1,k) - 2.0*v1(i,j,k) + v1(i,j-1,k))/(dy*dy)
d2vdz2=(v1(i,j,k+1) - 2.0*v1(i,j,k) + v1(i,j,k-1))/(dz*dz)

dv2dy=(0.25/dy)*( (v1(i,j,k)+v1(i,j+1,k))**2 &
		- (v1(i,j-1,k)+v1(i,j,k))**2) &
		+(0.25*Ga/dy)* &
		( ABS(v1(i,j,k)+v1(i,j+1,k))*(v1(i,j,k)-v1(i,j+1,k))&
		 - ABS(v1(i,j-1,k)+v1(i,j,k))*(v1(i,j-1,k)-v1(i,j,k)) )

dvudx=(0.25/dx)*( (u1(i,j,k)+u1(i,j+1,k))*(v1(i,j,k)+v1(i+1,j,k)) &
		- (u1(i-1,j,k)+u1(i-1,j+1,k))*(v1(i-1,j,k)+v1(i,j,k)) ) &
		+(0.25*Ga/dx)* &
		( ABS(u1(i,j,k)+u1(i,j+1,k))*(v1(i,j,k)-v1(i+1,j,k)) &
		- ABS(u1(i-1,j,k)+u1(i-1,j+1,k))*(v1(i-1,j,k)-v1(i,j,k)) )

dvwdz=(0.25/dz)*( (w1(i,j,k)+w1(i,j+1,k))*(v1(i,j,k)+v1(i,j,k+1)) &
		- (w1(i,j,k-1)+w1(i,j+1,k-1))*(v1(i,j,k-1)+v1(i,j,k)) ) &
		+(0.25*Ga/dz)* &
		( ABS(w1(i,j,k)+w1(i,j+1,k))*(v1(i,j,k)-v1(i,j,k+1)) &
		- ABS(w1(i,j,k-1)+w1(i,j+1,k-1))*(v1(i,j,k-1)-v1(i,j,k)) )

FLY=0.0

G(i,j,k)=v1(i,j,k)+dt*( (1./RE)*(d2vdx2 + d2vdy2 + d2vdz2)&
		-dv2dy -dvudx -dvwdz + FLY )
END IF
END IF
END IF
return
end subroutine

attributes(global) subroutine CALCL(L,U1,V1,W1)
implicit none

real,dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),intent(inout) :: U1,V1,W1,L
real :: d2wdx2,d2wdy2,d2wdz2,dwudx,dwvdy,dw2dz,FLZ

integer :: i,j,k

i=blockDim%x*(blockIdx%x-1)+threadIdx%x
j=blockDim%y*(blockIdx%y-1)+threadIdx%y
k=blockDim%z*(blockIdx%z-1)+threadIdx%z

IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.(kmax-1))) THEN

d2wdx2=(w1(i+1,j,k) - 2.0*w1(i,j,k) + w1(i-1,j,k))/(dx*dx)
d2wdy2=(w1(i,j+1,k) - 2.0*w1(i,j,k) + w1(i,j-1,k))/(dy*dy)
d2wdz2=(w1(i,j,k+1) - 2.0*w1(i,j,k) + w1(i,j,k-1))/(dz*dz)

dw2dz=(0.25/dz)*( (w1(i,j,k)+w1(i,j,k+1))**2 &
		- (w1(i,j,k-1)+w1(i,j,k))**2) &
		+(0.25*Ga/dz)* &
		( ABS(w1(i,j,k)+w1(i,j,k+1))*(w1(i,j,k)-w1(i,j,k+1))&
		 - ABS(w1(i,j,k-1)+w1(i,j,k))*(w1(i,j,k-1)-w1(i,j,k)) )

dwudx=(0.25/dx)*( (u1(i,j,k)+u1(i,j,k+1))*(w1(i,j,k)+w1(i+1,j,k)) &
		- (u1(i-1,j,k)+u1(i-1,j,k+1))*(w1(i-1,j,k)+w1(i,j,k)) ) &
		+(0.25*Ga/dx)* &
		( ABS(u1(i,j,k)+u1(i,j,k+1))*(w1(i,j,k)-w1(i+1,j,k)) &
		- ABS(u1(i-1,j,k)+u1(i-1,j,k+1))*(w1(i-1,j,k)-w1(i,j,k)) )

dwvdy=(0.25/dy)*( (v1(i,j,k)+v1(i,j,k+1))*(w1(i,j,k)+w1(i,j+1,k)) &
		- (v1(i,j-1,k)+v1(i,j-1,k+1))*(w1(i,j-1,k)+w1(i,j,k)) ) &
		+(0.25*Ga/dy)* &
		( ABS(v1(i,j,k)+v1(i,j,k+1))*(w1(i,j,k)-w1(i,j+1,k)) &
		- ABS(v1(i,j-1,k)+v1(i,j-1,k+1))*(w1(i,j-1,k)-w1(i,j,k)) )

FLZ=0.0

L(i,j,k)=w1(i,j,k)+dt*( (1./RE)*(d2wdx2 + d2wdy2 + d2wdz2)&
		-dw2dz -dwudx -dwvdy + FLZ )
END IF
END IF
END IF

return
end subroutine 

attributes(global) subroutine FGLBOUND(F,G,L,U2,V2,W2)
implicit none

real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: F,G,L,U2,V2,W2
integer :: i,j,k
i=blockDim%x*(blockIdx%x-1)+threadIdx%x
j=blockDim%y*(blockIdx%y-1)+threadIdx%y
k=blockDim%z*(blockIdx%z-1)+threadIdx%z

IF(i.eq.imin) then
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN 
		F(imin-1,J,k)=U2(imin-1,J,k)
		F(imax,J,k)=U2(imax,J,k)
END IF
END IF
END IF

IF(j.eq.jmin) then
IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((k.GE.kmax).AND.(k.LE.kmax)) THEN
		G(I,jmin-1,k)=V2(I,jmin-1,k)
		G(I,jmax,k)=V2(I,jmax,k)
END IF
END IF
END IF

IF(k.eq.kmin) then
IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
		L(I,J,kmin-1)=W2(I,J,kmin-1)
		L(I,J,kmax)=W2(I,J,kmax)
END IF
END IF
END IF

return
end subroutine

attributes(global) subroutine PRESSURE(p1,p2,F,G,L)
implicit none

real,dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),intent(inout) :: P1,P2,F,G,L
integer :: i,j,k,pt,pn
real :: pdt,HT,d2pdx2,d2pdy2,d2pdz2,dFdx,dGdy,dLdz

i=blockDim%x*(blockIdx%x-1)+threadIdx%x
j=blockDim%y*(blockIdx%y-1)+threadIdx%y
k=blockDim%z*(blockIdx%z-1)+threadIdx%z

pdt=0.000001
HT=1./dt

IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN

d2pdx2=(p2(i+1,j,k)-2.*p2(i,j,k)+p2(i-1,j,k))/(dx*dx)
d2pdy2=(p2(i,j+1,k)-2.*p2(i,j,k)+p2(i,j-1,k))/(dy*dy)
d2pdz2=(p2(i,j,k+1)-2.*p2(i,j,k)+p2(i,j,k-1))/(dz*dz)

dFdx=(F(i,j,k)-F(i-1,j,k))/dx
dGdy=(G(i,j,k)-G(i,j-1,k))/dy
dLdz=(L(i,j,k)-L(i,j,k-1))/dz

p2(i,j,k)=p1(i,j,k)+pdt*( d2pdx2+d2pdy2+d2pdz2-HT*(dFdx+dGdy+dLdz ) )
end if
end if
end if

IF(i.eq.imin) then
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN 
p2(imin-1,j,k)=p2(imin,j,k)
p2(imax+1,j,k)=p2(imax,j,k)
END IF
END IF
END IF

IF(j.eq.jmin) then
IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN
p2(i,jmin-1,k)=p2(i,jmin,k)
p2(i,jmax+1,k)=p2(i,jmax,k)
END IF
END IF
END IF

IF(k.eq.kmin) then
IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
p2(i,j,kmin-1)=p2(i,j,kmin)
p2(i,j,kmax+1)=p2(i,j,kmax)
END IF
END IF
END IF

return
end subroutine PRESSURE

attributes(global) subroutine U2V2W2(P2,F,G,L,U2,V2,W2)
implicit none

real,dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),intent(inout) ::&
	P2,F,G,L,U2,V2,W2
integer :: i,j,k

i=blockDim%x*(blockIdx%x-1)+threadIdx%x
j=blockDim%y*(blockIdx%y-1)+threadIdx%y
k=blockDim%z*(blockIdx%z-1)+threadIdx%z


IF((i.GE.imin).AND.(i.LE.(imax-1))) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN
	U2(i,j,k)=F(i,j,k)-(dt/dx)*(p2(i+1,j,k)-p2(i,j,k))
END IF
END IF
END IF

IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.(jmax-1))) THEN
IF((k.GE.kmin).AND.(k.LE.kmax)) THEN
	V2(i,j,k)=G(i,j,k)-(dt/dy)*(p2(i,j+1,k)-p2(i,j,k))
END IF
END IF
END IF

IF((i.GE.imin).AND.(i.LE.imax)) THEN
IF((j.GE.jmin).AND.(j.LE.jmax)) THEN
IF((k.GE.kmin).AND.(k.LE.(kmax-1))) THEN
	W2(i,j,k)=L(i,j,k)-(dt/dz)*(p2(i,j,k+1)-p2(i,j,k))
END IF
END IF
END IF
return
end subroutine

attributes(host) subroutine TIME_REMAINING(t2,t1,Tmax,time_r,t)
implicit none

real, value :: t2,t1
integer, value :: Tmax,t
character, intent(inout) :: time_r*12
character*7 :: ST
real :: time

time=(float(Tmax)/float(t)-1.0)*(t2-t1)

if(time.le.60.0) then
write(ST,'(I4)') int(time)
time_r=ST//' seg'
else if(((60.0).lt.time).and.(time.le.(60.0*60.0))) then
write(ST,'(F7.2)') time/60.0
time_r=ST//' min'
else if(((60.0*60.0).lt.time).and.(time.le.(60.0*60.0*24.0))) then
write(ST,'(F7.2)') time/3600.0
time_r=ST//' hrs'
else
write(ST,'(F7.2)') time/(3600.0*24.0)
time_r=ST//' days'
end if

return
end subroutine

end module kernel
!%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%
program CUDANS3D
use cudafor
use kernel
implicit none

integer, parameter :: nx=32 , ny=32 , nz=32
real, parameter :: LX=1.0 , LY=1.0 , LZ=1.0

integer :: i,j,k,t,Tmax,imin_h,imax_h,kmin_h,kmax_h,jmin_h,jmax_h,Nplot,np,NMaxP,cont
real :: dx_h,dy_h,dz_h,Re_h,dt_h,t1,t2
real, dimension(:), allocatable :: x,y,z
real, device, dimension(:), allocatable :: x_d,y_d,z_d

real, dimension(:,:,:), allocatable :: U1,U2,V1,V2,W1,W2,P1,P2,F,G,L
real, device, dimension(:,:,:), allocatable :: U1_d,U2_d,V1_d,V2_d,W1_d,W2_d,&
						P1_d,P2_d,F_d,G_d,L_d
character*12 :: time_r
 
integer(kind=cuda_stream_kind) :: stream1,stream2,stream3
type(dim3) :: grid,tBlock
integer :: istat

tBlock=dim3(16,16,4)
grid=dim3((nx+tBlock%x-1)/tBlock%x,(ny+tBlock%y-1)/tBlock%y,(nz+tBlock%z-1)/tBlock%z)

imin=1 ; jmin=1 ; kmin=1
imax=nx ; jmax=ny ; kmax=nz

imin_h=imin ; jmin_h=jmin ; kmin_h=kmin
imax_h=imax ; jmax_h=jmax ; kmax_h=kmax

dx=LX/nx ; dy=LY/ny ; dz=LZ/nz
dx_h=dx ; dy_h=dy ; dz_h=dz

RE=1.0 ; GA=0.7
RE_h=RE

NMaxP=100
Nplot=500
dt=0.5*(RE_h/2.0)/( (1.0/dx_h)**2 + (1.0/dy_h)**2 + (1.0/dz_h)**2) ; Tmax=nint(0.5/dt)
dt_h=dt
print*,'■  dt=', dt_h
print*,'■  Tmax=', Tmax

ALLOCATE(x(imin-1:imax+1),y(jmin-1:jmax+1),z(kmin-1:kmax+1))
ALLOCATE(u1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),u2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(v1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),v2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(w1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),w2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(p1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),p2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(F(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),G(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(L(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))

ALLOCATE(x_d(imin-1:imax+1),y_d(jmin-1:jmax+1),z_d(kmin-1:kmax+1))
ALLOCATE(u1_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),u2_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(v1_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),v2_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(w1_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),w2_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(p1_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),p2_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(F_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),G_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(L_d(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))

x(imin_h-1)=0.
x(imin_h)=0.5*dx_h
DO i=imin_h+1,imax_h
x(i)=0.5*dx_h+(i-1)*dx_h
END DO
x(imax_h+1)=LX

y(jmin_h-1)=0.
y(jmin_h)=0.5*dy_h
DO j=jmin_h+1,jmax_h
y(j)=0.5*dy_h+(j-1)*dy_h
END DO
y(jmax_h+1)=Ly

z(kmin_h-1)=0.
z(kmin_h)=0.5*dz_h
DO k=kmin_h+1,kmax_h
z(k)=0.5*dz_h+(k-1)*dz_h
END DO
z(kmax_h+1)=Lz

U1=0.0 ; V1=0.0 ; W1=0.0 ; P1=0.0
U2=U1 ; V2=V1 ; W2=W1 ; P2=P1

u1_d=u1 ; v1_d=v1 ; w1_d=w1 ; p1_d=p1
u2_d=u2 ; v2_d=v2 ; w2_d=w2 ; p2_d=p2

print*, '   STEP        U            V            W            P      TIME REMAINING'
print*, ' ======== ============ ============ ============ ============ ============'
call CPU_time(t1)
cont=1
DO t=1, Tmax
call BOUNDARY<<<grid,tBlock>>>(U1_d,U2_d,V1_d,V2_d,W1_d,W2_d)

istat = cudaStreamCreate(stream1)
istat = cudaStreamCreate(stream2)
istat = cudaStreamCreate(stream3)
call CALCF<<<grid,tBlock,0,stream1>>>(F_d,U1_d,V1_d,W1_d)
call CALCG<<<grid,tBlock,0,stream2>>>(G_d,U1_d,V1_d,W1_d)
call CALCL<<<grid,tBlock,0,stream3>>>(L_d,U1_d,V1_d,W1_d)
istat = cudaStreamDestroy(stream1)
istat = cudaStreamDestroy(stream2)
istat = cudaStreamDestroy(stream3)

call FGLBOUND<<<grid,tBlock>>>(F_d,G_d,L_d,U2_d,V2_d,W2_d)
do np=1,NMaxP
	call PRESSURE<<<grid,tBlock>>>(p1_d,p2_d,F_d,G_d,L_d)
	P1_d=P2_d
enddo

call U2V2W2<<<grid,tBlock>>>(P2_d,F_d,G_d,L_d,U2_d,V2_d,W2_d)

IF(mod(t,Nplot).EQ.0) THEN

call CPU_time(t2)
call TIME_REMAINING(t2,t1,Tmax,time_r,t)

 2    format(1x,I8,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,A)
	w2(nx/2,ny/2,nz/2)=w2_d(nx/2,ny/2,nz/2)
	u2(nx/2,ny/2,nz/2)=u2_d(nx/2,ny/2,nz/2)
	v2(nx/2,ny/2,nz/2)=v2_d(nx/2,ny/2,nz/2)
	P2(nx/2,ny/2,nz/2)=P2_d(nx/2,ny/2,nz/2)
	print 2,t,U2(nx/2,ny/2,nz/2),V2(nx/2,ny/2,nz/2),W2(nx/2,ny/2,nz/2),P2(nx/2,ny/2,nz/2),time_r
END IF

u1_d=u2_d
v1_d=v2_d
w1_d=w2_d
ENDDO

p2=p2_d
u1=u1_d
v1=v1_d
w1=w1_d

open(10,file='data.dat')
write(10,*) ' TITLE="NS3D" '
write(10,*) ' VARIABLES="X","Y","Z","U","V","W","PRESS","MAGU","GRADP" '
write(10,*) ' ZONE T="1", I=', imax_h,' J=', jmax_h,' K=', kmax_h

do k=kmin_h,kmax_h
do j=jmin_h,jmax_h
do i=imin_h,imax_h
	write(10,25) x(i),y(j),z(k),u1(i,j,k),v1(i,j,k),w1(i,j,k),p2(i,j,k),&
		sqrt(u1(i,j,k)**2 + v1(i,j,k)**2 + w1(i,j,k)**2),&
		sqrt(((p2(i+1,j,k)-p2(i-1,j,k))/(2.0*dx_h) )**2 +&
		 ((p2(i,j+1,k)-p2(i,j-1,k))/(2.0*dy_h) )**2 +&
		 ((p2(i,j,k+1)-p2(i,j,k-1))/(2.0*dz_h) )**2)
enddo
enddo
enddo 
close(10)
25    format(1x,E12.5,2x,E12.5,2x,E12.5,2x,E12.5,2x,E12.5,2x,E12.5,2x,E12.5,2x,E12.5,2x,E12.5)
print*, ' =================== END PROGRAM ============================'
end program
