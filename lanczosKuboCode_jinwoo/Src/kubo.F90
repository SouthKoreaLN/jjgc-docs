module kubo

   use mio
   use atoms,               only : in1, in2, nAt, inode1, inode2

   implicit none
   
   PRIVATE

!#ifdef MPI
!#define USE_MPI 1
!#include <sprng_f.h>
!#endif /* MPI */

   public :: KuboInitWF
   public :: KuboInitWFPDOS
   public :: KuboInitWFLayerDOS
   public :: KuboInitWFLayerAndSpeciesDOS
   public :: KuboDOS
   public :: KuboTEvol

contains

subroutine KuboInitWF(Psi)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID
   !use random,               only : RandNum, rand_t, RandSeed

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   logical PDOS,setSeed
   integer PDOSAtomNumber, seedValue

#ifdef DEBUG
   call MIO_Debug('KuboInitWF',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   call MIO_InputParameter('setSeed',setSeed,.false.)
   call MIO_InputParameter('seedValue',seedValue,123456)
   if (setSeed) then
       seed = seedValue
       print*, "seedValue considered from Gendata"
   else
       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   end if
   call random_seed(PUT = seed)

   !!$OMP PARALLEL DO 
   !call RandSeed(rng,nThread)
   !do i=1,nAt
   do i=inode1, inode2
      call random_number(r)
      !Psi(i) = exp(twopi*cmplx_i*RandNum(rng))/sqrt(real(nAt))
      Psi(i) = exp(twopi*cmplx_i*r)/sqrt(real(nAt))
   end do
   !!$OMP END PARALLEL DO

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWF',1)
#endif /* DEBUG */

end subroutine KuboInitWF

subroutine KuboInitWFPDOS(Psi,PDOSAtomNumber)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID
   !use random,               only : RandNum, rand_t, RandSeed

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer PDOSAtomNumber

#ifdef DEBUG
   call MIO_Debug('KuboInitWFPDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(PUT = seed)

   !!$OMP PARALLEL DO 
   !call RandSeed(rng,nThread)
   !do i=1,nAt
   do i=inode1, inode2
      Psi(i) = cmplx_0
   end do
   Psi(PDOSAtomNumber) = cmplx_1
   !!$OMP END PARALLEL DO

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWFPDOS',1)
#endif /* DEBUG */

end subroutine KuboInitWFPDOS

subroutine KuboInitWFLayerDOS(Psi,layerNumber,numberOfLayers,numberOfAtomsInLayer)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID
   use atoms,                only : layerIndex
   !use random,               only : RandNum, rand_t, RandSeed

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer :: SetSeed
   integer layerNumber, numberOfLayers, numberOfAtomsInLayer

#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   !seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call MIO_InputParameter('SeedSet',SetSeed,1235) 
   print*, "setting the RP seed to", SetSeed
   seed = SetSeed
   call random_seed(PUT = seed)

   print*, "inode1, inode2, nAt", inode1, inode2, nAt
   !!$OMP PARALLEL DO 
   !call RandSeed(rng,nThread)
   do i=1,nAt
   !do i=inode1, inode2
      if (layerIndex(i).eq.layerNumber) then
          call random_number(r)
          Psi(i) = exp(twopi*cmplx_i*r)/sqrt(real(numberOfAtomsInLayer))
      else
          Psi(i) = cmplx_0
      end if
   end do
   !!$OMP END PARALLEL DO

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerDOS',1)
#endif /* DEBUG */

end subroutine KuboInitWFLayerDOS

subroutine KuboInitWFLayerAndSpeciesDOS(Psi,layerNumber,speciesNumber,numberOfLayers,numberOfAtomsInLayer)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID
   use atoms,                only : layerIndex, Species
   !use random,               only : RandNum, rand_t, RandSeed

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer :: SetSeed
   integer layerNumber, numberOfLayers, numberOfAtomsInLayer, speciesNumber

#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerAndSpeciesDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   !seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call MIO_InputParameter('SeedSet',SetSeed,1235) 
   print*, "setting the RP seed to", SetSeed
   seed = SetSeed
   call random_seed(PUT = seed)

   print*, "inode1, inode2, nAt", inode1, inode2, nAt
   !!$OMP PARALLEL DO 
   !call RandSeed(rng,nThread)
   do i=1,nAt
   !do i=inode1, inode2
      if (layerIndex(i).eq.layerNumber .and. Species(i) .eq. speciesNumber) then
          call random_number(r)
          Psi(i) = exp(twopi*cmplx_i*r)/sqrt(real(numberOfAtomsInLayer))
      else
          Psi(i) = cmplx_0
      end if
   end do
   !!$OMP END PARALLEL DO

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerAndSpeciesDOS',1)
#endif /* DEBUG */

end subroutine KuboInitWFLayerAndSpeciesDOS

subroutine KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList,nRecurs,eps,nEn,Emin,Emax)

   use name,                 only : prefix
   use kubosubs,             only : KuboRecursion, KuboFrac, KuboFermi
   !use neigh,                only : nList

   !real(dp), parameter :: epsF=0.01_dp, EnF=10.0_dp
   !integer, parameter :: EptsF=20000

   complex(dp), intent(out) :: Psi(inode1:),Psin(inode1:),Psinm1(inode1:)
   real(dp), intent(out) :: a(nRecurs), b(nRecurs)
   complex(dp), intent(out) :: H(inode1:)
   real(dp), intent(in) :: H0(inode1:)
   complex(dp), intent(in) :: hopp(:,inode1:)
   integer, intent(in) :: NList(1:,inode1:)
   integer, intent(in) :: nRecurs, nEn
   real(dp), intent(in) :: eps
   real(dp) :: Emin, Emax
   real(dp) :: ac, bc

   character(len=100) :: flnm, flnm2

   logical :: calculateInterval

#ifdef DEBUG
   call MIO_Debug('KuboDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call MIO_print('Calculating DOS by recursion','kubo')
   call KuboRecursion(Psi,Psin,Psinm1,nRecurs,a,b,H,H0,hopp,nList,.true.,trim(prefix))
   ! DOS for all energy range to obtain Fermi energy
   !call KuboFrac(a,b,nRecurs,1.0_dp,epsF,EptsF,-EnF,EnF,'DAT.FERMI')
   !call KuboFermi(EptsF)
   ! The DOS with more resoulution in the range specified
   flnm = trim(prefix)//'.DOS'
   flnm2 = trim(prefix)//'.AB'
   !call MIO_InputParameter('CalculateInterval',calculateInterval,.false.)
   !if (calculateInterval) then
   !   call intervalle(a, b, ac, bc, nRecurs)
   !   Emin = ac-2.0d0*bc
   !   Emax = ac+2.0d0*bc
   !end if 
   call KuboFrac(a,b,nRecurs,1.0_dp,eps,nEn,Emin,Emax,flnm,flnm2)
   call MIO_Print('DOS written','kubo')
   call MIO_Print('')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboDOS',1)
#endif /* DEBUG */


end subroutine KuboDOS

subroutine KuboTEvol(Psi,ZUPsi,Psin,Psinm1,XpnPsi,XpnPsim1,c,dT,nT,az,bz,H,H0,hopp,NList,NeighD, &
                    ac,bc,nEn,Emin,Emax,eps,nRecurs,nPol,nWr)

   use constants,            only : cmplx_0
   use name,                 only : prefix
   use kubosubs,             only : KuboRecursion, KuboFrac, KuboEvol

   complex(dp), intent(inout) :: Psi(inode1:), ZUPsi(inode1:),Psin(inode1:)
   complex(dp), intent(inout) :: Psinm1(inode1:),XpnPsi(inode1:),XpnPsim1(inode1:)
   complex(dp), intent(in) :: c(nPol)
   real(dp), intent(in) :: dT, ac, bc, eps, Emin, Emax
   complex(dp), intent(out) :: H(inode1:)
   real(dp), intent(in) :: H0(inode1:), NeighD(1:,1:,inode1:)
   complex(dp), intent(in) :: hopp(:,inode1:)
   integer, intent(in) :: nT, nRecurs, nPol, nWr, NList(1:,inode1:), nEn
   real(dp), intent(inout) :: az(nRecurs), bz(nRecurs)

   logical :: prnt
   integer :: iT, i
   character(len=80) :: str, str2
   complex(dp) :: cnum, s

#ifdef DEBUG
   call MIO_Debug('KuboTEvol',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   prnt = .false.
   call MIO_Print('')
   call MIO_Print('Time evolution of wave function','kubo')
   call MIO_Print('')
   do iT=1,nT
      if (mod(iT,nWr)==0) then
         prnt = .true.
         write(str,'(a,i5,a)') '   *** Time Step ', iT, ' ***'
         call MIO_Print(trim(str),'kubo')
      end if
      call KuboEvol(Psi,ZUPsi,Psin,Psinm1,XpnPsi,XpnPsim1,c,H0,hopp,NList,NeighD,ac,bc,nPol,prnt)
      if (prnt) then
         cnum = cmplx_0
         !$OMP PARALLEL DO REDUCTION(+:cnum)
         do i=1,nAt
            cnum = cnum + ZUPsi(i)*conjg(ZUPsi(i))
         enddo
         !$OMP END PARALLEL DO
#ifdef MPI
         call MPIAllRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
         cnum = s
#endif /* MPI */
         !$OMP PARALLEL DO
         do i=1,nAt
            ZUPsi(i) = ZUPsi(i)/sqrt(abs(cnum))
         end do
         !$OMP END PARALLEL DO
         call KuboRecursion(ZUPsi,Psin,Psinm1,nRecurs,az,bz,H,H0,hopp,NList,.false.)
         !$OMP PARALLEL DO
         do i=1,nAt
            ZUPsi(i) = ZUPsi(i)*sqrt(abs(cnum))
         end do
         !$OMP END PARALLEL DO
         write(str,'(i0.4)') iT
         str2 = str
         str = trim(prefix)//'.Z2_'//trim(str)//'.dat'
         str2 = trim(prefix)//'.AB_'//trim(str2)//'.dat'
         call KuboFrac(az,bz,nRecurs,abs(cnum),eps,nEn,Emin,Emax,str,str2)
         prnt = .false.
      end if
   end do

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboTEvol',1)
#endif /* DEBUG */

end subroutine KuboTEvol

!!---------------------------------------------------------------------- 
!      FUNCTION pythag (a, b) 
!      DOUBLE PRECISION a, b, pythag 
!      DOUBLE PRECISION absa, absb 
!      absa = dabs (a) 
!      absb = dabs (b) 
!      IF (absa.gt.absb) THEN 
!         pythag = absa * dsqrt (1. + (absb / absa) **2) 
!      ELSE 
!         IF (absb.eq.0.) THEN 
!            pythag = 0. 
!         ELSE 
!            pythag = absb * dsqrt (1. + (absa / absb) **2) 
!         ENDIF 
!      ENDIF 
!      RETURN 
!      END FUNCTION pythag                           
!!-----------------------------------------------------------------------
!
!      SUBROUTINE intervalle (a, b, ac, bc, NRECURS)
!      IMPLICIT NONE
!      integer :: NRECURS
!      real(dp) :: a (NRECURS), b (NRECURS), ac, bc
!      !INTEGER NRECURS 
!      real(dp) :: am (NRECURS), bm (NRECURS), z (NRECURS, NRECURS)
!      integer ::  i, j
!      real(dp) :: emin, emax
!!  Fait appel a tqli.f, DOnc aussi a pythag.f (Numerical Recipes)       
!
!      z(:,:) = 0.
!
!      DO i = 1, NRECURS
!        z (i, i) = 1.
!      ENDDO
!
!      DO i = 1, NRECURS - 1
!      am (i) = a (i)
!      bm (i + 1) = b (i)
!      ENDDO
!      am (NRECURS) = a (NRECURS)
!
!      CALL tqli (am, bm, NRECURS, NRECURS, z)
!
!      emin = am (1)
!      emax = am (1)
!      DO i = 2, NRECURS
!      IF (am (i) .le.emin) emin = am (i)
!      IF (am (i) .ge.emax) emax = am (i)
!      ENDDO
!
!!  Intervalle [ac-2bc,ac+2bc] contenant tout le spectre (avec 10% de    
!!  marge) :                                                             
!      ac = (emax + emin) / 2.
!      bc = 1.1 * (emax - emin) / 4.
!      !print*, "inside subroutine place 5"
!
!      RETURN
!      END SUBROUTINE intervalle
!
!      SUBROUTINE tqli (d, e, n, np, z) 
!      use pythag
!      INTEGER n, np 
!      DOUBLE PRECISION d (np), e (np), z (np, np) 
!!U    USES pythag                                                       
!      INTEGER i, iter, k, l, m 
!      DOUBLE PRECISION b, c, dd, f, g, p, r, s, pythag 
!      DO 11 i = 2, n 
!         e (i - 1) = e (i) 
!   11 END DO 
!      e (n) = 0. 
!      DO 15 l = 1, n 
!         iter = 0 
!    1    DO 12 m = l, n - 1 
!            dd = dabs (d (m) ) + dabs (d (m + 1) ) 
!            IF (dabs (e (m) ) + dd.eq.dd) goto 2 
!   12    END DO 
!         m = n 
!    2    IF (m.ne.l) THEN 
!!         if(iter.eq.30)pause 'too many iterations in tqli'             
!            iter = iter + 1 
!            g = (d (l + 1) - d (l) ) / (2. * e (l) ) 
!            r = pythag (g, 1.d0) 
!            g = d (m) - d (l) + e (l) / (g + sign (r, g) ) 
!            s = 1. 
!            c = 1. 
!            p = 0. 
!            DO 14 i = m - 1, l, - 1 
!               f = s * e (i) 
!               b = c * e (i) 
!               r = pythag (f, g) 
!               e (i + 1) = r 
!               IF (r.eq.0.) THEN 
!                  d (i + 1) = d (i + 1) - p 
!                  e (m) = 0. 
!                  GOTO 1 
!               ENDIF 
!               s = f / r 
!               c = g / r 
!               g = d (i + 1) - p 
!               r = (d (i) - g) * s + 2. * c * b 
!               p = s * r 
!               d (i + 1) = g + p 
!               g = c * r - b 
!!     Omit lines from here ...                                          
!               DO 13 k = 1, n 
!                  f = z (k, i + 1) 
!                  z (k, i + 1) = s * z (k, i) + c * f 
!                  z (k, i) = c * z (k, i) - s * f 
!   13          END DO 
!!     ... to here when finding only eigenvalues.                        
!   14       END DO 
!            d (l) = d (l) - p 
!            e (l) = g 
!            e (m) = 0. 
!            GOTO 1 
!         ENDIF 
!   15 END DO 
!      RETURN 
!      END SUBROUTINE tqli                           
                                                                        

end module kubo
