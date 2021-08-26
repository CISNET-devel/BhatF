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

subroutine set_variables (nvars_total, x_labels, x_ini, is_fixed, x_est, x_min, x_max) bind(c,name="ftn_set_variables")
  implicit none 
  integer, intent(in) :: nvars_total
  real*8, intent(in) :: x_ini(nvars_total), x_est(nvars_total), x_min(nvars_total), x_max(nvars_total)
  character*1, intent(in) :: x_labels(10,nvars_total)
  integer, intent(in) :: is_fixed(nvars_total)

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

  iflag=1
  imcmc=0
  iboot=0

  call debug_print_commons()
  
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


!!! ***STOPPED HERE!*** Calls migrad with transformed variables
subroutine ftn_migrad(nvar, x_final, f_final)  bind(c,name="ftn_migrad")
  implicit none
  ! Arguments:
  integer, intent(inout) :: nvar ! in:how many variables allocated; out:how many variables valid (not more than allocated)
  double precision, intent(out) :: x_final(nvar), f_final
  ! Common blocks:
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
  ! Local vars:
  integer nct
  double precision :: x0(nmax)

  ! parameter transformation
  call ftrafo(ndim,x0,su)
  ! First function call (FIXME: is it really needed?)
  iflag = 1
  call func(su,ndim,f_final)

  nct = 1
  iflag = 2 
  call migrad(ndim, npar, nct, x0, f_final)

  call btrafo(ndim, x0, su)
  iflag = 3 ! FIXME: the original code is not consistent in this
  call func(su,ndim,f_final)

  if (nvar>ndim) nvar=ndim
  x_final(1:nvar) = su(1:nvar)

end subroutine ftn_migrad


! This subroutine is needed because the C function can only accept string arrays,
! but the Fortran code supplies fixed-length string scalars.
! The subroutine utilizes the permitted argument association between
! `character(len=n)` (a string variable) and `character*1(n)` (a string array)
! (see http://www.sternwarte.uni-erlangen.de/~ai05/vorlesungen/astrocomp/f77-standard.pdf,
!  paragraphs 15.9.3.1; 15.9.3.3; 17.1.2)
subroutine message_cpp_interface(message, msglen)
  implicit none
       integer, intent(in) :: msglen
       character(len=1), intent(in) :: message(msglen)
  interface
     subroutine message_callback(message, msglen) bind(c,name="message_callback")
       implicit none
       integer, intent(in) :: msglen
       character(len=1), intent(in) :: message(msglen)
     end subroutine message_callback
  end interface
  call message_callback(message, msglen)
end subroutine message_cpp_interface


subroutine logmessage(wrtbuffer, nlines)
  ! Log the lines from `wrtbuffer` until an empty line
  implicit none 
  character(len=*), intent(in) :: wrtbuffer(:)
  integer, intent(in), optional :: nlines
  integer i, nl, sl

  nl = size(wrtbuffer)
  if (present(nlines)) nl = min(nlines, nl)
  do i=1,nl
     sl = len_trim(wrtbuffer(i))
     if (sl==0) exit
     call message_cpp_interface(wrtbuffer(i), sl)
  end do
end subroutine logmessage


subroutine gnuplot_ini()
  implicit none
  write (*,*) 'DEBUG: Would have called GNUPLOT_INI here'
end subroutine gnuplot_ini


subroutine gnuplot()
  implicit none
  write (*,*) 'DEBUG: Would have called GNUPLOT here'
end subroutine gnuplot
