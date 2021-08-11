subroutine debug_print_commons()
  implicit none 
  integer, parameter :: nmax=100
  character*10 labels(nmax)
  character rnf(nmax)
  integer ivn(nmax)
  real*8 su, se, s1, s2
  integer ndim,npar,np,mp,np2,ndim2,npar2,nlm
  integer iflag,imcmc,iboot,iseed,ierr
  common /LABELS/ labels,rnf,ivn
  common /BOUNDS/ su(nmax),se(nmax),s1(nmax),s2(nmax)
  common /NPRODCTS/ ndim,npar,np,mp,np2,ndim2,npar2,nlm  
  common /IFLAGS/ iflag,imcmc,iboot,iseed,ierr  

  write (*,*) 'DEBUG: ndim=',ndim,'npar=',npar
  write (*,*) 'DEBUG: labels=',labels(1:ndim)
  write (*,*) 'DEBUG: rnf=',rnf(1:ndim)
  write (*,*) 'DEBUG: ivn=',ivn(1:ndim)
  write (*,*) 'DEBUG: su=',su(1:ndim)
  write (*,*) 'DEBUG: se=',se(1:ndim)
  write (*,*) 'DEBUG: s1=',s1(1:ndim)
  write (*,*) 'DEBUG: s2=',s2(1:ndim)
end subroutine debug_print_commons

subroutine set_variables (nvars_total, x_labels, x_ini, is_fixed, x_est, x_min, x_max, x0) bind(c,name="ftn_set_variables")
  implicit none 
  integer, intent(in) :: nvars_total
  real*8, intent(in) :: x_ini(nvars_total), x_est(nvars_total), x_min(nvars_total), x_max(nvars_total)
  character*1, intent(in) :: x_labels(10,nvars_total)
  integer, intent(in) :: is_fixed(nvars_total)
  integer, intent(out) :: x0(nvars_total)

  integer, parameter :: nmax=100
  character*1 labels(10,nmax) ! note: storage-conforming to char*10(nmax)
  character rnf(nmax)
  integer ivn(nmax)
  real*8 su, se, s1, s2
  integer ndim,npar,np,mp,np2,ndim2,npar2,nlm
  integer iflag,imcmc,iboot,iseed,ierr
  common /LABELS/ labels,rnf,ivn
  common /BOUNDS/ su(nmax),se(nmax),s1(nmax),s2(nmax)
  common /NPRODCTS/ ndim,npar,np,mp,np2,ndim2,npar2,nlm  
  common /IFLAGS/ iflag,imcmc,iboot,iseed,ierr  

  integer :: i,j

  labels(:,1:nvars_total) = x_labels(:,:)
  su(1:nvars_total) = x_ini(:)
  se(1:nvars_total) = x_est(:)
  s1(1:nvars_total) = x_min(:)
  s2(1:nvars_total) = x_max(:)

  where (is_fixed(:)/=0) rnf(1:nvars_total) = 'S'
  where (is_fixed(:)==0) rnf(1:nvars_total) = 'R'
  j = 1
  do i=1,nvars_total
     if (is_fixed(i)==0) ivn(j)=i
     j = j+1
  end do

  ndim = nvars_total
  npar = count(is_fixed(:)==0)

  np=npar
  mp=np+1
  np2=2*npar
  ndim2=ndim*ndim
  npar2=npar*npar
  nlm=(npar*npar+3*npar)/2

  call debug_print_commons()
  
  ! parameter transformation
  call ftrafo(ndim,x0,su)

  !xs0 = x0
  !nct=0    ! ... counts no. of FUNC calls

  ! first function call 
  iflag=1
  imcmc=0
  iboot=0
  !covm=0.
  
  !call func(su,ndim,f_start)
  !nct=nct+1

  !PRINT*,'first call with function value: ',F_START
  !IFLAG=2   
end subroutine set_variables

subroutine func(u,npar,f)
  implicit none
  integer, intent(in) :: npar
  double precision, intent(in) :: u(npar)
  double precision, intent(out) :: f
 
  integer :: iflag,imcmc,iboot,iseed,ierr
  common /IFLAGS/ iflag,imcmc,iboot,iseed,ierr

  interface
     subroutine cpp_callback(f, iflag, u, npar) bind(c,name="cpp_callback")
       implicit none
       real*8, intent(out) :: f
       integer, intent(in) :: iflag
       integer, intent(in) :: npar
       double precision, intent(in) :: u(npar)
     end subroutine cpp_callback
  end interface

  call cpp_callback(f, iflag, u, npar)
end subroutine func


subroutine gnuplot_ini()
  implicit none
  write (*,*) 'DEBUG: Would have called GNUPLOT_INI here'
end subroutine gnuplot_ini


subroutine gnuplot()
  implicit none
  write (*,*) 'DEBUG: Would have called GNUPLOT here'
end subroutine gnuplot
