module diag

   use mio

   implicit none

   PRIVATE

   complex(dp), pointer :: ZWork(:)=>NULL()
   real(dp), pointer :: DWork(:)=>NULL()
   integer :: lwork
   real(dp), save :: Efermi=0.0_dp, Emin=-50.0_dp, Emax=20.0_dp

   public :: DiagInit
   public :: DiagDOS, Diag3DBands, DiagPDOS, DiagChern
   public :: DiagBands, DiagBandsG, DiagBandsAroundK, DiagBandsRashba
   public :: DiagSpectralFunction, DiagSpectralFunctionKGrid
   public :: DiagSpectralFunctionKGridInequivalent
   public :: DiagSpectralFunctionKGridInequivalent_v2
   public :: DiagSpectralFunctionKGridInequivalentEnergyCut
   public :: DiagSpectralFunctionKGridInequivalentEnergyCut_v2
   public :: DiagSpectralFunctionKGridInequivalentEnergyCutNickDale

contains

subroutine DiagInit(N)

   integer, intent(in) :: N

   complex(dp) :: A(1,1), OPT(1)
   real(dp) :: W(1), W2(1)
   integer :: INFO

#ifdef DEBUG
   call MIO_Debug('DiagInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call ZHEEV('N','L',N,A,N,W,OPT,-1,W2,INFO)
   if (INFO /= 0) call MIO_Kill('Error in workspace query for diagonalization','diag','DiagInit')
   lwork = int(OPT(1))
   !allocate(ZWork(lwork))
   !allocate(DWork(3*N-2))
   call MIO_Allocate(ZWork,lwork,'ZWork','diag')
   call MIO_Allocate(DWork,3*N-2,'DWork','diag')

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagInit',1)
#endif /* DEBUG */

end subroutine DiagInit

subroutine DiagDOS()

   use cell,                 only : rcell, ucell
   use atoms,                only : nEl, nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : twopi
   use math

   integer, parameter :: intorder=5

   integer :: nk(3), ptot, i1, i2, i3, ik, Epts, u, is, uu
   complex(dp), pointer :: H(:,:,:)=>NULL()
   real(dp), pointer :: Kgrid(:,:)=>NULL(), Eig(:,:)=>NULL(), DOS(:,:)=>NULL(), E(:)=>NULL()
   real(dp), pointer :: EStore(:,:)=>NULL()
   type(cl_file) :: file, file2
   character(len=100) :: flnm, flnm2
   real(dp) :: eps, E1, E2, s, sp, Ep

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   complex(dp) :: HLoc(nAt, nAt)


#ifdef DEBUG
   call MIO_Debug('DiagDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_Print('Calculating DOS by diagonalization','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   DOS = 0.0d0
   call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(EStore,[ptot,nAt],'DOS','diag')
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do
   call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
   call MIO_Allocate(Eig,[nAt,nspin],'Eig','diag')
   flnm = trim(prefix)//'.diag.DOS'
   call file%Open(name=flnm,serial=.true.)
   u = file%GetUnit()
   Emax = -huge(1.0_dp)
   Emin = huge(1.0_dp)
   do i2=1,Epts
        E(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
   end do
   do is=1,nspin
      !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, ik, HLoc), REDUCTION(+:DOS), REDUCTION(min:Emin), REDUCTION(max:Emax), &
      !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
      do ik=1,ptot
         KptsLoc = Kgrid(:,ik)
         ELoc = 0.0_dp
         HLoc = 0.0_dp
         !if (modulo(ik,int(ptot/10)).eq.0) print*, "progress is: ", ik/int(ptot/10)*10, "percent"
         !write(*,*) 'k-pt', ik, ptot
         !call DiagHam(nAt,nspin,is,H(:,:,is),Eig(:,is),Kgrid(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         call DiagHam(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         !Eig(:,is) = Eig(:,is)*g0
         ELoc = ELoc*g0
         !Emin = min(Emin,Eig(1,is))
         !Emax = max(Emax,Eig(nAt,is))
         !do i1=1,nAt
         !   do i2=1,Epts
         !      E(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
         !   end do
         !end do
         Emin = min(Emin,ELoc(1))
         Emax = max(Emax,ELoc(nAt))
         !do i1=1,nAt
         !   do i2=1,Epts
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-ELoc(i1))**2/(2.0_dp*eps**2))
         !   end do
         !end do
         do i1=1,nAt
            EStore(ik,i1) = ELoc(i1)
         end do
      end do
      !$OMP END PARALLEL DO
   end do
   do is=1,nspin
      do ik = 1,ptot
         do i1=1,nAt
            do i2=1,Epts
               DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
            end do
         end do
      end do
   end do
   if (nspin==1) then
      DOS = 2.0_dp*DOS/(eps*sqrt(twopi)*ptot)
   else
      DOS = DOS/(eps*sqrt(twopi)*ptot)
   end if
   do i1=1,Epts
      write(u,*) E(i1), (DOS(i1,is),is=1,nspin)
   end do
   call file%Close()
   call MIO_Deallocate(Eig,'Eig','diag')
   call MIO_Deallocate(H,'H','diag')
   call MIO_Print('Emin: '//trim(num2str(Emin,4))//', Emax: '//trim(num2str(Emax,4)),'diag')
   eps = (E2-E1)/(Epts-1)
   sp = 0.0_dp
   Ep = Emin
   do ik=1,Epts-2*intorder
      s = 0.0_dp
      do is=1,nspin
         s = s + TrapezoidalInt(DOS(:intorder*2+ik,is),intorder*2+ik,eps,intorder)
      end do
      if (s>=nEl) then
         Efermi = (E(intorder*2+ik)+Ep)/2.0_dp
         exit
      else
         if (s/=sp) then
            Ep = E(intorder*2+ik)
            sp = s
         end if
      end if
   end do
   call MIO_Print('Efermi: '//trim(num2str(Efermi,5)),'diag')
   call MIO_Print('')

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagDOS',1)
#endif /* DEBUG */
   !

end subroutine DiagDOS

subroutine DiagPDOS()

   use cell,                 only : rcell, ucell
   use atoms,                only : nEl, nAt, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : twopi
   use math

   integer, parameter :: intorder=5

   integer :: nk(3), ptot, i1, i2, i3, ik, Epts, u, is, uu, i1bis, ivec, ivec2, PDOSLayerIndex, numberOfLayers
   complex(dp), pointer :: H(:,:,:)=>NULL()
   real(dp), pointer :: Kgrid(:,:)=>NULL(), Eig(:,:)=>NULL()
   real(dp), pointer :: DOS(:,:)=>NULL(), E(:)=>NULL()
   real(dp), pointer :: DOS_thread(:,:,:)=>NULL()
   real(dp), pointer :: EStore(:,:)=>NULL()
   type(cl_file) :: file, file2
   character(len=100) :: flnm, flnm2
   character(len=2) :: fmt, x1
   real(dp) :: eps, E1, E2, s, sp, Ep

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   complex(dp) :: HLoc(nAt, nAt)
   real(dp) :: pipj, vectormultip
   complex(dp), pointer :: HStore(:,:,:)=>NULL()

   integer :: omp_get_thread_num, omp_get_max_threads, index_ii


#ifdef DEBUG
   call MIO_Debug('DiagPDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_Print('Calculating PDOS by diagonalization','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(EStore,[ptot,nAt],'DOS','diag')
   call MIO_Allocate(HStore,[ptot,nAt,nAt],'DOS','diag')
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do
   call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
   call MIO_Allocate(Eig,[nAt,nspin],'Eig','diag')
   Emax = -huge(1.0_dp)
   Emin = huge(1.0_dp)
   do i2=1,Epts
        E(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
   end do
   do is=1,nspin
      !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, ik, HLoc), REDUCTION(min:Emin), REDUCTION(max:Emax), &
      !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
      do ik=1,ptot
         KptsLoc = Kgrid(:,ik)
         ELoc = 0.0_dp
         HLoc = 0.0_dp
         !if (modulo(ik,int(ptot/10)).eq.0) print*, "progress is: ", ik/int(ptot/10)*10, "percent"
         !write(*,*) 'k-pt', ik, ptot
         !call DiagHam(nAt,nspin,is,H(:,:,is),Eig(:,is),Kgrid(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         call DiagHamPDOS(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         !Eig(:,is) = Eig(:,is)*g0
         ELoc = ELoc*g0
         !Emin = min(Emin,Eig(1,is))
         !Emax = max(Emax,Eig(nAt,is))
         !do i1=1,nAt
         !   do i2=1,Epts
         !      E(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
         !   end do
         !end do
         Emin = min(Emin,ELoc(1))
         Emax = max(Emax,ELoc(nAt))
         !do i1=1,nAt
         !   do i2=1,Epts
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-ELoc(i1))**2/(2.0_dp*eps**2))
         !   end do
         !end do
         do i1=1,nAt
            EStore(ik,i1) = ELoc(i1)
            do i1bis=1,nAt
              HStore(ik,i1,i1bis)=HLoc(i1,i1bis)
            end do
         end do
      end do
      !$OMP END PARALLEL DO
   end do
   call MIO_InputParameter('numberOfLayers',numberOfLayers,2)
   do PDOSLayerIndex=1,numberOfLayers
      !fmt = '(I2.2)' ! an integer of width 5 with zeros at the left
      !write (x1,fmt) PDOSLayerIndex
      flnm = trim(prefix)//'.diag.DOS.Layer'//trim(num2str(PDOSLayerIndex))
      call file%Open(name=flnm,serial=.true.)
      u = file%GetUnit()
      call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
      call MIO_Allocate(DOS_thread,[Epts,nspin,omp_get_max_threads()],'DOS','diag')
      DOS = 0.0_dp
      DOS_thread = 0.0_dp
      do is=1,nspin
         do ik = 1,ptot ! Loop over k
            do i1=1,nAt ! Loop over atoms (same as number of bands here)
               if (layerIndex(i1).eq.PDOSLayerIndex) then
                  do i2=1,Epts
                     DOS(i2,is) = 0.0_dp
                     do ivec=1,nAt ! add this for the vector multiplication projection operator
                         !!$OMP PARALLEL DO PRIVATE(vectormultip, pipj), REDUCTION(+:DOS_thread), &
                         !$OMP PARALLEL DO PRIVATE(vectormultip, pipj), &
                         !$OMP& SHARED(E, EStore, eps, DOS_thread)
                         do ivec2=1,nAt ! add this for the vector multiplication projection operator
                     ! outer product between vector multiplication of two eigenvectors
                     ! and the eigenenergy vector. The outer product is thus over the
                     ! eigenenergies
                           vectormultip=conjg(HStore(ik,ivec,i1))*HStore(ik,ivec2,i1) ! 
                           pipj=EStore(ik,i1)*vectormultip
                           DOS_thread(i2,is,omp_get_thread_num()+1) = DOS_thread(i2,is,omp_get_thread_num()+1) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2)) * pipj
                           !DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2)) * pipj
                         end do
                         !$OMP END PARALLEL DO
                     end do
                  end do
               end if
            end do
         end do
      end do
      ! TO DO, add OMP parallelisaton here if it works without
      DO index_ii = 1,omp_get_max_threads()
         do is=1,nspin
            do i2=1,Epts
               DOS(i2,is) = DOS(i2,is) + DOS_thread(i2,is,index_ii)
            enddo
         enddo
      ENDDO
      call MIO_Deallocate(DOS_thread,'DOS','diag')
      if (nspin==1) then
         DOS = 2.0_dp*DOS/(eps*sqrt(twopi)*ptot)
      else
         DOS = DOS/(eps*sqrt(twopi)*ptot)
      end if
      do i1=1,Epts
         write(u,*) E(i1), (DOS(i1,is),is=1,nspin)
      end do
      call file%Close()
      call MIO_Deallocate(Eig,'Eig','diag')
      call MIO_Deallocate(H,'H','diag')
      call MIO_Print('Emin: '//trim(num2str(Emin,4))//', Emax: '//trim(num2str(Emax,4)),'diag')
      eps = (E2-E1)/(Epts-1)
      sp = 0.0_dp
      Ep = Emin
      do ik=1,Epts-2*intorder
         s = 0.0_dp
         do is=1,nspin
            s = s + TrapezoidalInt(DOS(:intorder*2+ik,is),intorder*2+ik,eps,intorder)
         end do
         if (s>=nEl) then
            Efermi = (E(intorder*2+ik)+Ep)/2.0_dp
            exit
         else
            if (s/=sp) then
               Ep = E(intorder*2+ik)
               sp = s
            end if
         end if
      end do
      call MIO_Print('Efermi: '//trim(num2str(Efermi,5)),'diag')
      call MIO_Print('')
      call MIO_Deallocate(DOS,'DOS','diag')
   end do

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagPDOS',1)
#endif /* DEBUG */
   !

end subroutine DiagPDOS

subroutine Diag3DBands()

   use cell,                 only : rcell, ucell
   use atoms,                only : nEl, nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : twopi
   use math

   integer, parameter :: intorder=5

   integer :: nk(3), ptot, i1, i2, i3, ik, Epts, u, is, uu, i
   complex(dp), pointer :: H(:,:,:)=>NULL()
   real(dp), pointer :: Kgrid(:,:)=>NULL(), Eig(:,:)=>NULL(), DOS(:,:)=>NULL(), E(:)=>NULL()
   type(cl_file) :: file
   character(len=100) :: flnm, flnm2
   real(dp) :: eps, E1, E2, s, sp, Ep

#ifdef DEBUG
   call MIO_Debug('Diag3DBands',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_Print('Calculating 3D Bands by diagonalization','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do
   call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
   call MIO_Allocate(Eig,[nAt,nspin],'Eig','diag')
   !flnm = trim(prefix)//'.diag.DOS'
   !call file%Open(name=flnm,serial=.true.)
   !u = file%GetUnit()
   !Emax = -huge(1.0_dp)
   !Emin = huge(1.0_dp)


   flnm2 = trim(prefix)//'.3Dbands'
   uu=98
   open(uu,FILE=flnm2,STATUS='replace')
   !write(uu,'(f16.8)') Efermi
   !write(uu,'(2f16.8)') 0.0_dp, d
   !write(uu,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
   write(uu,'(3i8)') nAt, nspin, ptot
   do is=1,nspin
      do ik=1,ptot
         if (modulo(ik,int(ptot/10)).eq.0) call MIO_Print('progress is: '//trim(num2str((ik/ptot/10.0_dp*10.0_dp),4))//' percent','diag')
         !write(*,*) 'k-pt', ik, ptot
         call DiagHam(nAt,nspin,is,H(:,:,is),Eig(:,is),Kgrid(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         !hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0))
         !lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
         !v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         !d = d + sqrt(dot_product(v,v))
         write(uu,'(f12.6,10f14.6,/,(10x,10f14.6))') Kgrid(:,ik),(Eig(i,is)*g0, i=1,nAt)
         !write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)

         Eig(:,is) = Eig(:,is)*g0
         Emin = min(Emin,Eig(1,is))
         Emax = max(Emax,Eig(nAt,is))
         !do i1=1,nAt
         !   do i2=1,Epts
         !      E(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
         !   end do
         !end do
      end do
   end do
   !if (nspin==1) then
   !   DOS = 2.0_dp*DOS/(eps*sqrt(twopi)*ptot)
   !else
   !   DOS = DOS/(eps*sqrt(twopi)*ptot)
   !end if
   !do i1=1,Epts
   !   write(u,*) E(i1), (DOS(i1,is),is=1,nspin)
   !end do
   !call file%Close()
   call MIO_Deallocate(Eig,'Eig','diag')
   call MIO_Deallocate(H,'H','diag')
   !call MIO_Print('Emin: '//trim(num2str(Emin,4))//', Emax: '//trim(num2str(Emax,4)),'diag')
   !eps = (E2-E1)/(Epts-1)
   !sp = 0.0_dp
   !Ep = Emin
   !do ik=1,Epts-2*intorder
   !   s = 0.0_dp
   !   do is=1,nspin
   !      s = s + TrapezoidalInt(DOS(:intorder*2+ik,is),intorder*2+ik,eps,intorder)
   !   end do
   !   if (s>=nEl) then
   !      Efermi = (E(intorder*2+ik)+Ep)/2.0_dp
   !      exit
   !   else
   !      if (s/=sp) then
   !         Ep = E(intorder*2+ik)
   !         sp = s
   !      end if
   !   end if
   !end do
   !call MIO_Print('Efermi: '//trim(num2str(Efermi,5)),'diag')
   !call MIO_Print('')

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('Diag3DBands',1)
#endif /* DEBUG */
   !

end subroutine Diag3DBands

subroutine DiagBands()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, uu, uuu, uuuu, neig
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, useSameNumberOfPoints
   real(dp) :: theta, volume

   real(dp) :: gcell(3,3), grcell(3,3)

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   !complex(dp), allocatable :: HLoc(:,:)
   complex(dp), pointer :: HLoc(:,:)=>NULL()
   logical :: useDifferentLatticeVectors, keepWaveFunction, sparseDiagSolver, useTAPW

#ifdef DEBUG
   call MIO_Debug('DiagBands',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   call MIO_InputParameter('Bands.SparseNeig',neig,100)
   call MIO_InputParameter('Bands.UseDifferentLatticeVectors',useDifferentLatticeVectors,.false.)
   call MIO_InputParameter('Bands.UseSameNumberOfPoints',useSameNumberOfPoints,.false.)
   if (useDifferentLatticeVectors) then
       ucell(1,2) = 0.0_dp
       rcell(2,1) = 0.0_dp
   end if
   !print*, "ucell at beginning of diagBands"
   !print*, ucell(:,1)
   !print*, ucell(:,2)
   !print*, ucell(:,3)
   !print*, "rcell at beginning of diagBands"
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         !print*, "0: ", path(:,ip)
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         !print*, "1: ", path(:,ip)
         !if (MoireBS) then
         !   !print*, "theta=", theta
         !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
         !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
         !end if
         !print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else 
         !call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               if (useSameNumberOfPoints) then
                  nPts(ip) = nPts0
                  ptsTot = ptsTot + nPts0
               else
                  nPts(ip) = nint(real(d*nPts0)/real(d0))
                  ptsTot = ptsTot + nPts(ip)
               end if
               !ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      !print*, "ptsTot, nPath", ptsTot, nPath, nPts(1)
      Kpts(:,1) = path(:,1)
      !ip = 0
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
      !do i=2,nPath
         do j=1,nPts(i)
            ip = ip + 1
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      !print*, "initial Gamma: ", Kpts(:,1)
      !print*, "final Gamma: ", Kpts(:,ip), ip
      !call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt,nspin,ptsTot+1],'E','diag')
      flnm = trim(prefix)//'.ham'
      uu=101
      open(uu,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.prob'
      uuuu=103
      open(uuuu,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.eig'
      uuu=102
      open(uuu,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot+1
      call MIO_InputParameter('keepWaveFunction',keepWaveFunction,.false.)
      call MIO_InputParameter('sparseDiagSolver',sparseDiagSolver,.false.)
      call MIO_InputParameter('useTAPW',useTAPW,.false.)
      if (.not. sparseDiagSolver) then
         call MIO_Allocate(HLoc,[nAt,nAt],'HLoc','diag')
        !call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      end if
      !if (useTAPW) then
      !   allocate(eigvecs(M, M_band, nkx * nky))
      !end if
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      E = 0.0_dp
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !!$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, is, ip, HLoc), &
      !!$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kpts)
      !do ip=1,ptsTot + 1
      !   do is=1,nspin
      !      !if (modulo(ip,int(ptsTot/10)).eq.0) print*, "progress is: ", ip/int(ptsTot/10)*10, "percent"
      !      ELoc = 0.0_dp
      !      HLoc = 0.0_dp
      !      KptsLoc = Kpts(:,ip)
      !      if (keepWaveFunction) then
      !         call DiagHamWF(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
      !         do i=1,nAt
      !             !write(uu,*) (HLoc(j,i),j=1,nAt)
      !             if (ip.eq.1) then ! for now, let's only print the firt k-point (cfr ptsTot+1)
      !                 write(uu,'(20000(F10.5,1X))') (real(HLoc(j,i)),j=1,nAt)
      !                 write(uuuu,'(20000(F10.5,1X))') (real(conjg(HLoc(j,i))*HLoc(j,i)),j=1,nAt)
      !                 write(uuu,*) ELoc(i)
      !             end if
      !         end do
      !      else if (sparseDiagSolver) then
      !         !call DiagHamSparse(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
      !         call DiagHamSparse(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
      !         !call DiagHamSparse2(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)


      !      else
      !         call DiagHam(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
      !      end if
      !      E(:,is,ip) = ELoc
      !   end do
      !end do
      !!$OMP END PARALLEL DO
      if (ptsTot == 1) then
         ! Special case: only one k-point, do not use OpenMP here
         do ip = 1, ptsTot
            do is = 1, nspin
               ELoc = 0.0_dp
               if (.not. sparseDiagSolver) then
                  HLoc = 0.0_dp
               end if
               KptsLoc = Kpts(:, ip)
               if (keepWaveFunction) then
                  call DiagHamWF(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
                  do i = 1, nAt
                     if (ip == 1) then
                        write(uu, '(20000(F10.5,1X))') (real(HLoc(j, i)), j = 1, nAt)
                        write(uuuu, '(20000(F10.5,1X))') (real(conjg(HLoc(j, i)) * HLoc(j, i)), j = 1, nAt)
                        write(uuu, *) ELoc(i)
                     end if
                  end do
               else if (sparseDiagSolver) then
                  call DiagHamSparse(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else if (useTAPW) then
                  !call DiagH0TAPW(nAt, nspin, is, ELoc, eigvec, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
                  call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
                  !eigvecs(:, :, ik) = ZWorkLoc(:, 1:M_band)
               else
                  call DiagHam(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
               end if
               E(:, is, ip) = ELoc
            end do
         end do
         do ip=1,ptsTot
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
            d = d + sqrt(dot_product(v,v))
            write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
            !write(u,*) d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
            if (sum(nPts(:nPath))==ip-1) then
               nPath = nPath+1
               call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
                 '   '//trim(num2str(d,6)),'diag')
            end if
         end do
      else if (sparseDiagSolver) then
         ! Many k-points, but parallelisation happens inside DiagHamSparse so no OMP here (some obsolete code here, too lazy to
         ! remove the parts that will not be considered because we are only doing the sparseDiagSolver part here
         do ip = 1, ptsTot
            do is = 1, nspin
               ELoc = 0.0_dp
               if (.not. sparseDiagSolver) then
                  HLoc = 0.0_dp
               end if
               KptsLoc = Kpts(:, ip)
               if (keepWaveFunction) then
                  call DiagHamWF(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
                  do i = 1, nAt
                     if (ip == 1) then
                        write(uu, '(20000(F10.5,1X))') (real(HLoc(j, i)), j = 1, nAt)
                        write(uuuu, '(20000(F10.5,1X))') (real(conjg(HLoc(j, i)) * HLoc(j, i)), j = 1, nAt)
                        write(uuu, *) ELoc(i)
                     end if
                  end do
               else if (sparseDiagSolver) then
                  call DiagHamSparse(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else if (useTAPW) then
                  !call DiagH0TAPW(nAt, nspin, is, ELoc, eigvec, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
                  call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else
                  call DiagHam(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
               end if
               E(:, is, ip) = ELoc
            end do
         end do
         do ip=1,ptsTot
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
            d = d + sqrt(dot_product(v,v))
            write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
            !write(u,*) d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
            if (sum(nPts(:nPath))==ip-1) then
               nPath = nPath+1
               call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
                 '   '//trim(num2str(d,6)),'diag')
            end if
         end do
      else 
         ! General case: more than one k-point, use OpenMP
         !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, is, ip, HLoc), &
         !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kpts)
         do ip = 1, ptsTot + 1
            do is = 1, nspin
               ELoc = 0.0_dp
               if (.not. sparseDiagSolver) then
                  HLoc = 0.0_dp
               end if
               KptsLoc = Kpts(:, ip)
               if (keepWaveFunction) then
                  call DiagHamWF(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
                  do i = 1, nAt
                     if (ip == 1) then
                        write(uu, '(20000(F10.5,1X))') (real(HLoc(j, i)), j = 1, nAt)
                        write(uuuu, '(20000(F10.5,1X))') (real(conjg(HLoc(j, i)) * HLoc(j, i)), j = 1, nAt)
                        write(uuu, *) ELoc(i)
                     end if
                  end do
               else if (sparseDiagSolver) then
                  call DiagHamSparse(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else if (useTAPW) then
                  !call DiagH0TAPW(nAt, nspin, is, ELoc, eigvec, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
                  print*, "calling TAPW routine for", ip
                  call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else
                  call DiagHam(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
               end if
               E(:, is, ip) = ELoc
            end do
         end do
         !$OMP END PARALLEL DO
         do ip=1,ptsTot + 1
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
            d = d + sqrt(dot_product(v,v))
            write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
            !write(u,*) d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
            if (sum(nPts(:nPath))==ip-1) then
               nPath = nPath+1
               call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
                 '   '//trim(num2str(d,6)),'diag')
            end if
         end do
      end if
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      !call MIO_Deallocate(H,'H','diag')
      !call MIO_Deallocate(Kpts,'Htsp','diag')
      !call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBands',1)
#endif /* DEBUG */

end subroutine DiagBands

subroutine DiagChern()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, uu, uuu, uuuu
   real(dp), pointer :: path(:,:)=>NULL(), Kgrid(:,:)=>NULL(), Chern(:,:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   integer :: nk(3), i1, i2, i3, ik
   real(dp) :: d0, v(3), d, hv, lc, totalChern
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, useSameNumberOfPoints
   real(dp) :: theta, volume, dArea

   real(dp) :: gcell(3,3), grcell(3,3)

   real(dp) :: KptsLoc(3)
   real(dp) :: ChernLoc(nAt)
   complex(dp) :: HLoc(nAt, nAt)
   logical :: useDifferentLatticeVectors, keepWaveFunction

#ifdef DEBUG
   call MIO_Debug('DiagChern',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   call MIO_InputParameter('Bands.UseDifferentLatticeVectors',useDifferentLatticeVectors,.false.)
   call MIO_InputParameter('Bands.UseSameNumberOfPoints',useSameNumberOfPoints,.false.)
   if (useDifferentLatticeVectors) then
       ucell(1,2) = 0.0_dp
       rcell(2,1) = 0.0_dp
   end if
   !print*, "ucell at beginning of diagBands"
   !print*, ucell(:,1)
   !print*, ucell(:,2)
   !print*, ucell(:,3)
   !print*, "rcell at beginning of diagBands"
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)
   !if (MIO_InputFindBlock('Bands.Path',nPath)) then
      !call MIO_Print('Band calculation','diag')
      !call MIO_Allocate(path,[3,nPath],'path','diag')
      !call MIO_InputBlock('Bands.Path',path)
      !do ip=1,nPath
      !   !print*, "0: ", path(:,ip)
      !   path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
      !   !print*, "1: ", path(:,ip)
      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      !   !print*, "2: ", path(:,ip)
      !end do
      !if (nPath==1) then
      !   call MIO_Allocate(nPts,1,'nPts','diag')
      !   nPts(1) = 1
      !   ptsTot = 1
      !else 
      !   !call MIO_Allocate(nPts,nPath-1,'nPts','diag')
      !   call MIO_Allocate(nPts,nPath,'nPts','diag')
      !   nPts(1) = nPts0
      !   ptsTot = nPts0
      !   if (nPath > 2) then
      !      v = path(:,2) - path(:,1)
      !      d0 = sqrt(dot_product(v,v))
      !      do ip=2,nPath-1 ! coz 4 points, 3 segments
      !         v = path(:,ip+1) - path(:,ip)
      !         d = sqrt(dot_product(v,v))
      !         if (useSameNumberOfPoints) then
      !            nPts(ip) = nPts0
      !            ptsTot = ptsTot + nPts0
      !         else
      !            nPts(ip) = nint(real(d*nPts0)/real(d0))
      !            ptsTot = ptsTot + nPts(ip)
      !         end if
      !         !ptsTot = ptsTot + nPts(ip)
      !      end do
      !      nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
      !   end if
      !end if
      !!call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      !call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      !!print*, "ptsTot, nPath", ptsTot, nPath, nPts(1)
      !Kpts(:,1) = path(:,1)
      !!ip = 0
      !ip = 1
      !d = 0.0_dp
      !do i=1,nPath-1
      !!do i=2,nPath
      !   do j=1,nPts(i)
      !      ip = ip + 1
      !      !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
      !      Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
      !      v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
      !      d = d + sqrt(dot_product(v,v))
      !   end do
      !end do

      call MIO_Print('Calculating Chern number by diagonalization','diag')
      call MIO_InputParameter('KGrid',nk,[1,1,1])
      ptsTot = nk(1)*nk(2)*nk(3)
      call MIO_Allocate(Kgrid,[3,ptsTot],'Kgrid','diag')
      ik = 0
      do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
         ik = ik+1
         Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
           rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
      end do; end do; end do

      !print*, "initial Gamma: ", Kpts(:,1)
      !print*, "final Gamma: ", Kpts(:,ip), ip
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(Chern,[nAt,nspin,ptsTot+1],'Chern','diag')
      !flnm = trim(prefix)//'.ham'
      !uu=101
      !open(uu,FILE=flnm,STATUS='replace')
      !flnm = trim(prefix)//'.prob'
      !uuuu=103
      !open(uuuu,FILE=flnm,STATUS='replace')
      !flnm = trim(prefix)//'.eig'
      !uuu=102
      !open(uuu,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.chern'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot
      call MIO_InputParameter('keepWaveFunction',keepWaveFunction,.false.)
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      Chern = 0.0_dp
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !$OMP PARALLEL DO PRIVATE(ChernLoc, KptsLoc, is, ip, HLoc), &
      !$OMP& SHARED(Chern, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kgrid)
      do ip=1,ptsTot
         do is=1,nspin
            !if (modulo(ip,int(ptsTot/10)).eq.0) print*, "progress is: ", ip/int(ptsTot/10)*10, "percent"
            ChernLoc = 0.0_dp
            !HLoc = 0.0_dp
            KptsLoc = Kgrid(:,ip)
            !if (keepWaveFunction) then
            !   call DiagHamWF(nAt,nspin,is,HLoc,ChernLoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !   do i=1,nAt
            !       !write(uu,*) (HLoc(j,i),j=1,nAt)
            !       if (ip.eq.1) then ! for now, let's only print the firt k-point (cfr ptsTot+1)
            !           write(uu,'(20000(F10.5,1X))') (real(HLoc(j,i)),j=1,nAt)
            !           write(uuuu,'(20000(F10.5,1X))') (real(conjg(HLoc(j,i))*HLoc(j,i)),j=1,nAt)
            !           write(uuu,*) ELoc(i)
            !       end if
            !   end do
            !else
               call DiagHamChern(nAt,nspin,is,HLoc,ChernLoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !end if
            Chern(:,is,ip) = ChernLoc
         end do
      end do
      !$OMP END PARALLEL DO
      do ip=1,ptsTot
         !v = Kgrid(:,ip) - Kgrid(:,max(ip-1,1))
         !!!v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
         !d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') Kgrid(:,ip),((Chern(i,is,ip), i=1,nAt),is=1,nspin)
         !write(u,*) d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
         !if (sum(nPts(:nPath))==ip-1) then
         !   nPath = nPath+1
         !   call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
         !     '   '//trim(num2str(d,6)),'diag')
         !end if
      end do
      dArea = norm(CrossProd(rcell(:,1),rcell(:,2)))/(nk(1)*nk(2))
      do i=1,nAt
         totalChern = 0.0_dp
         do ip=1,ptsTot
            !dArea = abs(rcell(1,1)*rcell(2,2) - rcell(2,1)*rcell(1,2))/(nk(1)*nk(2))
            totalChern = totalChern+Chern(i,1,ip)*dArea
         end do
         print*, "Chern number for band ", i, " equals ", totalChern
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(Chern,'Chern','diag')
      call MIO_Deallocate(H,'H','diag')
      !call MIO_Deallocate(Kpts,'Htsp','diag')
      !call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   !end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagChern',1)
#endif /* DEBUG */

end subroutine DiagChern

subroutine DiagBandsRashba()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS
   real(dp) :: theta, volume

   real(dp) :: gcell(3,3), grcell(3,3)

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt*2)
   complex(dp) :: HLoc(nAt*2, nAt*2)
   logical :: useDifferentLatticeVectors

#ifdef DEBUG
   call MIO_Debug('DiagBandsRashba',0)
#endif /* DEBUG */
#ifdef TIMER
   !call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   call MIO_InputParameter('Bands.UseDifferentLatticeVectors',useDifferentLatticeVectors,.false.)
   if (useDifferentLatticeVectors) then
       ucell(1,2) = 0.0_dp
       rcell(2,1) = 0.0_dp
   end if
   !print*, "ucell at beginning of diagBands"
   !print*, ucell(:,1)
   !print*, ucell(:,2)
   !print*, ucell(:,3)
   !print*, "rcell at beginning of diagBands"
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         !print*, "0: ", path(:,ip)
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         !print*, "1: ", path(:,ip)
         !if (MoireBS) then
         !   !print*, "theta=", theta
         !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
         !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
         !end if
         !print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else 
         !call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      !print*, "ptsTot, nPath", ptsTot, nPath, nPts(1)
      Kpts(:,1) = path(:,1)
      !ip = 0
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
      !do i=2,nPath
         do j=1,nPts(i)
            ip = ip + 1
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      !print*, "initial Gamma: ", Kpts(:,1)
      !print*, "final Gamma: ", Kpts(:,ip), ip
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt*2,nspin,ptsTot+1],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt*2, nspin, ptsTot+1
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      E = 0.0_dp
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, is, ip, HLoc), &
      !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kpts)
      do ip=1,ptsTot + 1
         do is=1,nspin
            !if (modulo(ip,int(ptsTot/10)).eq.0) print*, "progress is: ", ip/int(ptsTot/10)*10, "percent"
            ELoc = 0.0_dp
            HLoc = 0.0_dp
            KptsLoc = Kpts(:,ip)
            call DiagHamRashba(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            E(:,is,ip) = ELoc
         end do
      end do
      !$OMP END PARALLEL DO
      do ip=1,ptsTot + 1
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         !v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt*2),is=1,nspin)
         if (sum(nPts(:nPath))==ip-1) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      !call MIO_Deallocate(Kpts,'Htsp','diag')
      !call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   !call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBandsRashba',1)
#endif /* DEBUG */

end subroutine DiagBandsRashba

subroutine DiagBandsAroundK()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS
   real(dp) :: theta, volume, vn(3), grapheneK(3)

   real(dp) :: gcell(3,3), grcell(3,3)

#ifdef DEBUG
   call MIO_Debug('DiagBandsAroundK',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !print*, [aG,0.0_dp,0.0_dp]
   !print*, (/aG,0.0_dp,0.0_dp/)
   !gcell(:,1) = [aG,0.0_dp,0.0_dp]
   !gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   !gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   gcell(:,1) = (/aG,0.0_dp,0.0_dp/)
   gcell(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp/)
   gcell(:,3) = (/0.0_dp,0.0_dp,40.0_dp/)

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)

   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   !print*, "rcell"
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)
   !print*, "grcell"
   !print*, grcell(:,1)
   !print*, grcell(:,2)
   !print*, grcell(:,3)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         !print*, "0: ", path(:,ip)
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         grapheneK = 2.0/3.0*grcell(:,1) + 1.0/3.0*grcell(:,2) + 0.0*grcell(:,3) 
         path(:,ip) = path(:,ip) + grapheneK
         !print*, "1: ", path(:,ip)
         !if (MoireBS) then
         !   !print*, "theta=", theta
         !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
         !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
         !end if
         !print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else 
         !call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      !print*, "ptsTot, nPath", ptsTot, nPath, nPts(1)
      Kpts(:,1) = path(:,1)
      !ip = 0
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
      !do i=2,nPath
         do j=1,nPts(i)
            ip = ip + 1
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      !print*, "initial Gamma: ", Kpts(:,1)
      !print*, "final Gamma: ", Kpts(:,ip), ip
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot+1
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      do ip=1,ptsTot + 1
         do is=1,nspin
            call DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !do i=1,nAt
            !!   if (E(i,is)<=Efermi) hv = max(hv,E(i,is))
            !!   if (E(i,is)>Efermi) lc = min(lc,E(i,is))
            !!end do
            hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0))
            lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
         end do
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         !v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         if (sum(nPts(:nPath))==ip-1) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      !call MIO_Deallocate(Kpts,'Htsp','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBandsAroundK',1)
#endif /* DEBUG */

end subroutine DiagBandsAroundK

subroutine DiagBandsG()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS
   real(dp) :: theta, volume

   real(dp) :: gcell(3,3), grcell(3,3), vn(3)

#ifdef DEBUG
   call MIO_Debug('DiagBandsG',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)

   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   !print*, "graphene cell: "
   !print*, grcell(:,1)
   !print*, grcell(:,2)
   !print*, grcell(:,3)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         !print*, "0: ", path(:,ip)
         !path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         path(:,ip) = path(1,ip)*grcell(:,1) + path(2,ip)*grcell(:,2) + path(3,ip)*grcell(:,3)
         !print*, "1: ", path(:,ip)
         !if (MoireBS) then
         !   !print*, "theta=", theta
         !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
         !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
         !end if
         !print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else 
         !call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      !print*, "ptsTot, nPath", ptsTot, nPath, nPts(1)
      Kpts(:,1) = path(:,1)
      !ip = 0
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
      !do i=2,nPath
         do j=1,nPts(i)
            ip = ip + 1
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      !print*, "initial Gamma: ", Kpts(:,1)
      !print*, "final Gamma: ", Kpts(:,ip), ip
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot+1
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      do ip=1,ptsTot + 1
         do is=1,nspin
            call DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !do i=1,nAt
            !!   if (E(i,is)<=Efermi) hv = max(hv,E(i,is))
            !!   if (E(i,is)>Efermi) lc = min(lc,E(i,is))
            !!end do
            hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0))
            lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
         end do
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         !v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         if (sum(nPts(:nPath))==ip-1) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      !call MIO_Deallocate(Kpts,'Htsp','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBandsG',1)
#endif /* DEBUG */

end subroutine DiagBandsG

subroutine DiagSpectralFunction()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   !use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   !real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()

   real(dp) :: GVec(3) !, G1(3), G2(3)

   real(dp) :: eps, factor
   integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   !real(dp) :: area, volume, grcell(3,3), aG
   !real(dp) :: gcell(3,3), vn(3)

   !integer :: cellSize

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunction',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)

   !call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   !gcell(:,1) = [aG,0.0_dp,0.0_dp]
   !gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   !gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   !vn = CrossProd(gcell(:,1),gcell(:,2))
   !volume = dot_product(gcell(:,3),vn)
   !area = norm(vn)
   !grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   !grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   !grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   !print*, grcell(:,1)
   !print*, grcell(:,2)
   !print*, grcell(:,3)

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      !call MIO_Allocate(pathG,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      !call MIO_InputBlock('Spectral.Path',pathG)
      do ip=1,nPath
         !print*, "0: ", path(:,ip)
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         !pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
         !print*, "1: ", path(:,ip)
         !if (MoireBS) then
         !   !print*, "theta=", theta
         !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
         !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
         !end if
         !print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      Kpts(:,1) = path(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      KptsG(:,1) = Kpts(:,1) + GVec 
      !print*, "yup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
     
      !print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      !G1 = matmul(rcell,[1,1,0])
      !G2 = matmul(rcell,[1,1,0])
      !call MIO_InputParameter('CellSize', cellSize, 1)
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            !KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            !call random_number(randu)
            !randi = FLOOR(22*randu) 
            !call random_number(randu)
            !randj = FLOOR(22*randu) 
            !GVec = matmul(rcell,[randi,randj,0]) ! we only want to translate them by one reciprocal lattice vector
            KptsG(:,ip) = Kpts(:,ip) + GVec ! Kpts is in SC, Kpts is for Graphene (PC)
            !print*, "yup"
            !print*, Kpts(:,ip), KptsG(:,ip), GVec
            !print*, KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !print*, mod(KptsG(2,ip),G2)
            !Kpts(1,ip) = KptsG(1,ip) - FLOOR(KptsG(1,ip)/G1) * G1
            !Kpts(2,ip) = KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !Kpts(3,ip) = KptsG(3,ip) 
            !Kpts(1,ip) = mod(KptsG(1,ip),G1(1))
            !Kpts(2,ip) = mod(KptsG(2,ip),G2(2))
            !Kpts(3,ip) = KptsG(3,ip) 
            !Kpts(1,ip) = KptsG(1,ip)/cellSize
            !Kpts(2,ip) = KptsG(2,ip)/cellSize
            !Kpts(3,ip) = KptsG(3,ip) 
            !v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !d = d + sqrt(dot_product(v,v))
         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1],[ptsTot,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'Ake','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      !print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !print*, gaussian
      Pkc = 0.0_dp ! spectral weight PkscI(k)
      !$OMP PARALLEL DO PRIVATE (is, ik, nPath, iee, ie)
      do ik=1,ptsTot ! k loop
         !do ip=1,ptsTot ! kc loop
            !call DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !call DiagSpectralWeight(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ip),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !call DiagSpectralWeightNishi(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ik),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            call DiagSpectralWeightWeiKu(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ik),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !do i=1,nAt
            !!   if (E(i,is)<=Efermi) hv = max(hv,E(i,is))
            !!   if (E(i,is)>Efermi) lc = min(lc,E(i,is))
            !!end do
            !hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0)) ! mask restrict search for E smaller than Efermi
            !lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
            !v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
            !d = d + sqrt(dot_product(v,v))
            !write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                   is = 1
                   !if (useGaussianBroadening) then
                   !   !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                   !   Ake(ik,iee) = Ake(ik,iee) + exp(-(E(ie,is) - Energy(iee))**2.0/(2.0_dp*eps**2))
                   !else
                   if(abs(E(ie,is) - Energy(iee)).lt.(0.005/g0)) then
                        Ake(ik,iee) = Ake(ik,iee) + abs(Pkc(ik,ie))**2
                   end if
                   !end if
                end do
            end do
            !do ie=1,nAt   ! epsilonIksc
            !   is = 1
            !   !if(abs(E(ie,is) - Energy(iee)).lt.(0.1/g0)) then
            !   iee = floor(E(ie,is)/(E2-E1)*Epts)+Epts/2
            !   if(iee.lt.1 .or. iee.gt.Epts) then
            !       cycle
            !   else
            !       Ake(ik,iee) = Ake(ik,iee) + Pkc(ik,ie)
            !   end if
            !end do
            !do ie=1,nAt   ! epsilonIksc
            !end do
            !print*, SIZE(real(Ake(ik,:))), SIZE(gaussian)
            AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
            !AkeGaussian(ik,:) = convolve(gaussian(:),real(Ake(ik,:)),Epts)
            !do i1=1,nAt
            !   do i2=1,Epts
            !      Energy(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
            !      !DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
            !      Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)*exp(-(Energy(i2)-E(i1,is))**2/(2.0_dp*eps**2)) 
            !      !Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)
            !   end do
            !end do
         !end do ! this loop has to be finished first as we are still adding Pkcs to Ake
         !do ip=1,ptsTot ! kc loop
         !   v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         !   d = d + sqrt(dot_product(v,v))
         !   do iee=1,Epts  ! epsilon
         !       ! write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         !       ! float of lenght 12 with 6 after the comma
         !       ! repeat float of lengt 14 with 6 after the commq 10 times
         !       ! go to next line
         !       ! 10 empty spaces
         !       ! repeat float of lengt 14 with 6 after the commq 10 times, as many times as needed because of ()
         !       ! Ill make it simpler, but maybe bigger file
         !       write(u,'(f12.6,f12.6,f12.6)') d, Energy(iee), REAL(Ake(ik,iee))
         !   end do
         !end do
      end do
      !$OMP END PARALLEL DO
      do ik=1,ptsTot ! k loop
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake(ik,iee))
            end do
         end if
      end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunction',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunction 

subroutine DiagSpectralFunctionKGrid()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()

   real(dp) :: GVec(3) , G1(3), G2(3)

   real(dp) :: eps, factor
   integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)

   integer :: cellSize

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGrid',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)

   call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   !print*, grcell(:,1)
   !print*, grcell(:,2)
   !print*, grcell(:,3)

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_Allocate(pathG,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      call MIO_InputBlock('Spectral.Path',pathG)
      do ip=1,nPath
         !print*, "0: ", path(:,ip)
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
         !print*, "1: ", path(:,ip)
         !if (MoireBS) then
         !   !print*, "theta=", theta
         !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
         !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
         !end if
         !print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      !Kpts(:,1) = path(:,1)
      KptsG(:,1) = path(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      !KptsG(:,1) = Kpts(:,1) + GVec 
      !Kpts(:,1) = KptsG(:,1) - GVec 
      call MIO_InputParameter('CellSize', cellSize, 1)
      !print*, "cell size =", cellSize
      Kpts(1,1) = KptsG(1,1)/cellSize
      Kpts(2,1) = KptsG(2,1)/cellSize
      Kpts(3,1) = KptsG(3,1) 
      !print*, "yup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
     
      !print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G1 = matmul(rcell,[1,1,0])
      G2 = matmul(rcell,[1,1,0])
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            !call random_number(randu)
            !randi = FLOOR(22*randu) 
            !call random_number(randu)
            !randj = FLOOR(22*randu) 
            !GVec = matmul(rcell,[randi,randj,0]) ! we only want to translate them by one reciprocal lattice vector
            !KptsG(:,ip) = Kpts(:,ip) + GVec ! Kpts is in SC, Kpts is for Graphene (PC)
            !print*, "yup"
            !print*, Kpts(:,ip), KptsG(:,ip), GVec
            !print*, KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !print*, mod(KptsG(2,ip),G2)

            !Kpts(:,ip) = KptsG(:,ip) - GVec ! this one gives the same as when starting from Kpts
            !Kpts(1,ip) = KptsG(1,ip) - FLOOR(KptsG(1,ip)/G1) * G1
            !Kpts(2,ip) = KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !Kpts(3,ip) = KptsG(3,ip) 
            !Kpts(1,ip) = mod(KptsG(1,ip),G1(1))
            !Kpts(2,ip) = mod(KptsG(2,ip),G2(2))
            !Kpts(3,ip) = KptsG(3,ip) 
            Kpts(1,ip) = KptsG(1,ip)/cellSize
            Kpts(2,ip) = KptsG(2,ip)/cellSize
            Kpts(3,ip) = KptsG(3,ip) 
            !v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !d = d + sqrt(dot_product(v,v))
         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptsTot,nAt,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(Ake1,[ptsTot,Epts],'Ake1','diag')
      call MIO_Allocate(AkeGaussian1,[ptsTot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(Ake2,[ptsTot,Epts],'Ake2','diag')
      call MIO_Allocate(AkeGaussian2,[ptsTot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      Ake1 = 0.0_dp
      AkeGaussian1 = 0.0_dp
      Ake2 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      !print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !print*, gaussian
      Pkc = 0.0_dp ! spectral weight PkscI(k)
      do ik=1,ptsTot ! k loop
         !do ip=1,ptsTot ! kc loop
            !call DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !call DiagSpectralWeight(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ip),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !call DiagSpectralWeightNishi(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ik),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            call DiagSpectralWeightWeiKu(nAt,nspin,is,Pkc(ik,:,:),E(:,is),Kpts(:,ik),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !print*, "1 :", ik, Pkc(ik,nAt,1)
            !print*, "2 :", ik, Pkc(ik,nAt,2)
            !do i=1,nAt
            !!   if (E(i,is)<=Efermi) hv = max(hv,E(i,is))
            !!   if (E(i,is)>Efermi) lc = min(lc,E(i,is))
            !!end do
            !hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0)) ! mask restrict search for E smaller than Efermi
            !lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
            !v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
            !d = d + sqrt(dot_product(v,v))
            !write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
            !do iee=1,Epts  ! epsilon
            !    do ie=1,nAt   ! epsilonIksc
            !       is = 1
            !       !if (ie.eq.1) then
            !           if(abs(E(ie,is) - Energy(iee)).lt.(0.005/g0)) then
            !                Ake(ik,iee) = Ake(ik,iee) + abs(Pkc(ik,ie))**2
            !           end if
            !       !end if
            !    end do
            !end do
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                   is = 1
                   !if (ie.eq.1) then
                       if(abs(E(ie,is) - Energy(iee)).lt.(0.005/g0)) then
                            Ake1(ik,iee) = Ake1(ik,iee) + abs(Pkc(ik,ie,1))**2
                            Ake2(ik,iee) = Ake2(ik,iee) + abs(Pkc(ik,ie,2))**2
                       end if
                   !end if
                end do
            end do
            Ake = Ake1 + Ake2
            !do ie=1,nAt   ! epsilonIksc
            !   is = 1
            !   !if(abs(E(ie,is) - Energy(iee)).lt.(0.1/g0)) then
            !   iee = floor(E(ie,is)/(E2-E1)*Epts)+Epts/2
            !   if(iee.lt.1 .or. iee.gt.Epts) then
            !       cycle
            !   else
            !       Ake(ik,iee) = Ake(ik,iee) + Pkc(ik,ie)
            !   end if
            !end do
            !do ie=1,nAt   ! epsilonIksc
            !end do
            !print*, SIZE(real(Ake(ik,:))), SIZE(gaussian)
            AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
            !AkeGaussian(ik,:) = convolve(gaussian(:),real(Ake(ik,:)),Epts)
            !do i1=1,nAt
            !   do i2=1,Epts
            !      Energy(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
            !      !DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
            !      Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)*exp(-(Energy(i2)-E(i1,is))**2/(2.0_dp*eps**2)) 
            !      !Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)
            !   end do
            !end do
         !end do ! this loop has to be finished first as we are still adding Pkcs to Ake
         !do ip=1,ptsTot ! kc loop
         !   v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         !   d = d + sqrt(dot_product(v,v))
         !   do iee=1,Epts  ! epsilon
         !       ! write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         !       ! float of lenght 12 with 6 after the comma
         !       ! repeat float of lengt 14 with 6 after the commq 10 times
         !       ! go to next line
         !       ! 10 empty spaces
         !       ! repeat float of lengt 14 with 6 after the commq 10 times, as many times as needed because of ()
         !       ! Ill make it simpler, but maybe bigger file
         !       write(u,'(f12.6,f12.6,f12.6)') d, Energy(iee), REAL(Ake(ik,iee))
         !   end do
         !end do
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake(ik,iee))
            end do
         end if
      end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGrid',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGrid

subroutine DiagSpectralFunctionKGridInequivalent()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1,uu2,uu3,uu4,uuu2,uuuu2,uuuuu2,uuuuuu2,uuuuuuu2,uuuuuuuu2
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   
   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   !complex(dp), pointer :: PkcLoc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake2C(:,:)=>NULL()
   complex(dp), pointer :: Ake2D(:,:)=>NULL()
   complex(dp), pointer :: Ake2E(:,:)=>NULL()
   complex(dp), pointer :: Ake2F(:,:)=>NULL()
   complex(dp), pointer :: Ake2G(:,:)=>NULL()
   complex(dp), pointer :: Ake2H(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()
 
   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)
   complex(dp) :: PkcLocC(nAt,4)
   complex(dp) :: PkcLocD(nAt,4)
   complex(dp) :: PkcLocE(nAt,4)
   complex(dp) :: PkcLocF(nAt,4)
   complex(dp) :: PkcLocG(nAt,4)
   complex(dp) :: PkcLocH(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id
   real(dp) :: ELoc1(nAt),ELoc2(nAt)

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, foldByOne, useCoordinates
   logical :: foldByZero
   
   real(dp) :: topBottomRatio
   logical :: WeiKu, Nishi, Lee, WeiKuOld, useGaussianBroadening
 

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)
   call MIO_InputParameter('Spectral.Lee',Lee,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       !gcell(1,:) = ucell(1,:)/mmm(1)/2.0
       !gcell(2,:) = ucell(2,:)/mmm(2)/2.0
       !gcell(3,:) = ucell(3,:)
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       gcell1 = gcell
       if (rotateRefSystem .eq. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           phi = -phi*180.0_dp/pi
           aa = phi*pi/180.0_dp
           !print*, "Angle= ", aa, phi
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eq. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       !print*, "Angle= ", aa, phi
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if


   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume
   gcell2 = gcell

   print*, grcell(:,1)
   print*, grcell(:,2)
   print*, grcell(:,3)

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_Allocate(pathG,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      call MIO_InputBlock('Spectral.Path',pathG)
      do ip=1,nPath
         !print*, "0: ", path(:,ip)
         if (useCoordinates) then
            path(:,ip) = path(:,ip)
            pathG(:,ip) = pathG(:,ip)
         else
            path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
            pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
         end if
         !print*, "1: ", path(:,ip)
         !print*, "2: ", path(:,ip)
      end do
      print*, "The k-points for the reference system equal in cartesian coordinates:", pathG
      print*, "The k-points for the moire system equal in cartesian coordinates:", path
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      !Kpts(:,1) = path(:,1)
      KptsG(:,1) = pathG(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      !print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0]) ! Lattice vectors of moire cell
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      !print*, ll, kk
      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1) 
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2) 
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1)
         kpts(2,1) = KptsG(2,1)
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      end if
 
      !print*, "yup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
      !rcellInv = matinv3(rcell)
     
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            !print*, ll, kk
            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1) 
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2) 
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip) 
               kpts(2,ip) = KptsG(2,ip) 
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            end if

         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectralA'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')
      
      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2C'
      uuu2=302
      open(uuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2D'
      uuuu2=402
      open(uuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2E'
      uuuuu2=502
      open(uuuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2F'
      uuuuuu2=602
      open(uuuuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2G'
      uuuuuuu2=702
      open(uuuuuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2H'
      uuuuuuuu2=802
      open(uuuuuuuu2,FILE=flnm,STATUS='replace')
      
      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')
      
      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')
      
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptsTot,nAt,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptsTot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptsTot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      !print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !call MIO_Allocate(Ake1Loc,Epts,'Ake1Loc','diag')
      !call MIO_Allocate(Ake2Loc,Epts,'Ake2Loc','diag')
      !Ake1Loc = 0.0_dp
      !Ake2Loc = 0.0_dp
      call MIO_Allocate(Ake1,[ptsTot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake1B,[ptsTot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2,[ptsTot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake2B,[ptsTot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake2C,[ptsTot,Epts],'Ake2C','diag')
      call MIO_Allocate(Ake2D,[ptsTot,Epts],'Ake2D','diag')
      call MIO_Allocate(Ake2E,[ptsTot,Epts],'Ake2E','diag')
      call MIO_Allocate(Ake2F,[ptsTot,Epts],'Ake2F','diag')
      call MIO_Allocate(Ake2G,[ptsTot,Epts],'Ake2G','diag')
      call MIO_Allocate(Ake2H,[ptsTot,Epts],'Ake2H','diag')
      call MIO_Allocate(Ake3,[ptsTot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake3B,[ptsTot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4,[ptsTot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake4B,[ptsTot,Epts],'Ake4B','diag')
      Ake1  = 0.0_dp
      Ake1B = 0.0_dp
      Ake2  = 0.0_dp
      Ake2B = 0.0_dp
      Ake2C = 0.0_dp
      Ake2D = 0.0_dp
      Ake2E = 0.0_dp
      Ake2F = 0.0_dp
      Ake2G = 0.0_dp
      Ake2H = 0.0_dp
      Ake3  = 0.0_dp
      Ake3B = 0.0_dp
      Ake4  = 0.0_dp
      Ake4B = 0.0_dp
      !print*, gaussian
      !Pkc = 0.0_dp ! spectral weight PkscI(k)
      !print*, "0",  Pkc
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, PkcLocC,PkcLocD,PkcLocE,PkcLocF,PkcLocG,PkcLocH, KptsLoc), REDUCTION(+:Ake1Loc, Ake2Loc), &  
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake) 
      do ik=1,ptsTot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         PkcLocB = 0.0_dp
         PkcLocC = 0.0_dp
         PkcLocD = 0.0_dp
         PkcLocE = 0.0_dp
         PkcLocF = 0.0_dp
         PkcLocG = 0.0_dp
         PkcLocH = 0.0_dp
         KptsLoc = Kpts(:,ik)
         !KptsLoc = KptsG(:,ik)-G11
            if (WeiKu) then
                if (WeiKuOld) then
                    call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                    call DiagSpectralWeightWeiKuInequivalentMoreOrbitals(nAt,nspin,is,PkcLocA,PkcLocB,PkcLocC,PkcLocD,PkcLocE,PkcLocF,PkcLocG,PkcLocH,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Lee) then
                call DiagSpectralWeightWeiKuInequivalentLee(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       !print*,"species and layer", Species(ie), layerIndex(ie)
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                             Ake4(ik,iee) = Ake4(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,4))**2
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2
                                 
                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2
                                 ! C sublattice
                                 Ake2C(ik,iee) = Ake2C(ik,iee) + abs(PkcLocC(ie,2))**2
                                 ! D sublattice
                                 Ake2D(ik,iee) = Ake2D(ik,iee) + abs(PkcLocD(ie,2))**2
                                 ! E sublattice
                                 Ake2E(ik,iee) = Ake2E(ik,iee) + abs(PkcLocE(ie,2))**2
                                 ! F sublattice
                                 Ake2F(ik,iee) = Ake2F(ik,iee) + abs(PkcLocF(ie,2))**2
                                 ! G sublattice
                                 Ake2G(ik,iee) = Ake2G(ik,iee) + abs(PkcLocG(ie,2))**2
                                 ! H sublattice
                                 Ake2H(ik,iee) = Ake2H(ik,iee) + abs(PkcLocH(ie,2))**2
                             end if
                          end if
                       else if (Lee) then
                          !if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                          if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + real(PkcLocA(ie,1))
                             Ake2(ik,iee) = Ake2(ik,iee) + real(PkcLocA(ie,2))
                             Ake3(ik,iee) = Ake3(ik,iee) + real(PkcLocA(ie,3))
                             Ake4(ik,iee) = Ake4(ik,iee) + real(PkcLocA(ie,4))
                             !Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2
                             !Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2
                             !Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                             Ake4(ik,iee) = Ake4(ik,iee) + PkcLocA(ie,4)
                          end if
                       end if
                end do
            end do
         !end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:) + Ake4(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      do ik=1,ptsTot ! K loop
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                !if (Lee) then
                !   write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake(ik,iee))**2
                !   write(u1,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake1(ik,iee))**2
                !   write(u2,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake2(ik,iee))**2
                !   write(u3,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake3(ik,iee))**2
                !else
                   !write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake(ik,iee))
                   ! A sublattice
                   write(u1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1(ik,iee))
                   write(u2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2(ik,iee))
                   write(u3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3(ik,iee))
                   write(u4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4(ik,iee))
                   ! B sublattice
                   write(uu1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1B(ik,iee))
                   write(uu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2B(ik,iee))
                   write(uu3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3B(ik,iee))
                   write(uu4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4B(ik,iee))
                   ! C sublattice
                   write(uuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2C(ik,iee))
                   ! D sublattice
                   write(uuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2D(ik,iee))
                   ! E sublattice
                   write(uuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2E(ik,iee))
                   ! F sublattice
                   write(uuuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2F(ik,iee))
                   ! G sublattice
                   write(uuuuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2G(ik,iee))
                   ! H sublattice
                   write(uuuuuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2H(ik,iee))
                !endif 
            end do
         end if
      end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalent

subroutine DiagSpectralFunctionKGridInequivalent_v2()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1,uu2,uu3,uu4
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   
   logical :: MoireBS, GaussConv, DirectBand
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   !complex(dp), pointer :: PkcLoc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: EAke(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()
 
   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id
   real(dp) :: ELoc1(nAt),ELoc2(nAt)

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, foldByOne, useCoordinates
   logical :: foldByZero
   
   real(dp) :: topBottomRatio
   logical :: WeiKu, Nishi, Lee, WeiKuOld, useGaussianBroadening
 

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent_v2',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   print*, rcell(:,1)
   print*, rcell(:,2)
   print*, rcell(:,3)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)
   call MIO_InputParameter('Spectral.Lee',Lee,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       !gcell(1,:) = ucell(1,:)/mmm(1)/2.0
       !gcell(2,:) = ucell(2,:)/mmm(2)/2.0
       !gcell(3,:) = ucell(3,:)
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       gcell1 = gcell
       if (rotateRefSystem .eq. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           phi = -phi*180.0_dp/pi
           aa = phi*pi/180.0_dp
           !print*, "Angle= ", aa, phi
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eq. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       !print*, "Angle= ", aa, phi
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if


   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume
   gcell2 = gcell

   print*, grcell(:,1)
   print*, grcell(:,2)
   print*, grcell(:,3)

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_Allocate(pathG,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      call MIO_InputBlock('Spectral.Path',pathG)
      do ip=1,nPath
         print*, "0: ", path(:,ip)
         if (useCoordinates) then
            path(:,ip) = path(:,ip)
            pathG(:,ip) = pathG(:,ip)
         else
            path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
            pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
         end if
         print*, "1: ", path(:,ip)
         print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      !Kpts(:,1) = path(:,1)
      KptsG(:,1) = pathG(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0]) ! Lattice vectors of moire cell
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      !print*, ll, kk
      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1) 
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2) 
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1)
         kpts(2,1) = KptsG(2,1)
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      end if
 
      !print*, "yup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
      !rcellInv = matinv3(rcell)
     
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            !print*, ll, kk
            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1) 
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2) 
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip) 
               kpts(2,ip) = KptsG(2,ip) 
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            end if

         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral_v2'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')
      
      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')
      
      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')
      
      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')
      
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptsTot,nAt,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptsTot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptsTot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.DirectBandVal',DirectBand,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !call MIO_Allocate(Ake1Loc,Epts,'Ake1Loc','diag')
      !call MIO_Allocate(Ake2Loc,Epts,'Ake2Loc','diag')
      !Ake1Loc = 0.0_dp
      !Ake2Loc = 0.0_dp
      call MIO_Allocate(Ake1,[ptsTot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake1B,[ptsTot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2,[ptsTot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake2B,[ptsTot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake3,[ptsTot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake3B,[ptsTot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4,[ptsTot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake4B,[ptsTot,Epts],'Ake4B','diag')
      Ake1  = 0.0_dp
      Ake1B = 0.0_dp
      Ake2  = 0.0_dp
      Ake2B = 0.0_dp
      Ake3  = 0.0_dp
      Ake3B = 0.0_dp
      Ake4  = 0.0_dp
      Ake4B = 0.0_dp
      !print*, gaussian
      !Pkc = 0.0_dp ! spectral weight PkscI(k)
      !print*, "0",  Pkc
      if (DirectBand) then
        call MIO_Allocate(EAke,[ptsTot,Epts],'EAke','diag')
        EAke = 0.0_dp
      end if
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), REDUCTION(+:Ake1Loc, Ake2Loc), &  
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake) 
      do ik=1,ptsTot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
         !KptsLoc = KptsG(:,ik)-G11
            if (WeiKu) then
                if (WeiKuOld) then
                    call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                    call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Lee) then
                call DiagSpectralWeightWeiKuInequivalentLee(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                          else if (DirectBand) then
                             if (iee == int(ie-(nAt/2-Epts/2))) then
                                 EAke(ik,iee) = ELoc(ie)
                                 
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio 
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio
                                 
                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio 
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                             end if
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2
                                 
                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2
                             end if
                          end if
                       else if (Lee) then
                          !if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                          if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + real(PkcLocA(ie,1))
                             Ake2(ik,iee) = Ake2(ik,iee) + real(PkcLocA(ie,2))
                             Ake3(ik,iee) = Ake3(ik,iee) + real(PkcLocA(ie,3))
                             !Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2
                             !Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2
                             !Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                          end if
                       end if
                end do
            end do
         !end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      do ik=1,ptsTot ! K loop
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else if (DirectBand) then
            do iee=1,Epts
                !write(u,'(3f12.6,f12.6,10f12.6)') KptsG(:,ik), REAL(EAke(ik,iee)), REAL(Ake1(ik,iee)), REAL(Ake1B(ik,iee)),&
                !                                                                 & REAL(Ake2(ik,iee)), REAL(Ake2B(ik,iee)),&
                !                                                                 & REAL(Ake3(ik,iee)), REAL(Ake3B(ik,iee)),&
                !                                                                 & REAL(Ake4(ik,iee)), REAL(Ake4B(ik,iee))
                write(u,'(f12.6,f12.6,10f12.6)') d,  REAL(EAke(ik,iee)),           REAL(Ake1(ik,iee)), REAL(Ake1B(ik,iee)),&
                                                                                 & REAL(Ake2(ik,iee)), REAL(Ake2B(ik,iee)),&
                                                                                 & REAL(Ake3(ik,iee)), REAL(Ake3B(ik,iee)),&
                                                                                 & REAL(Ake4(ik,iee)), REAL(Ake4B(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                !if (Lee) then
                !   write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake(ik,iee))**2
                !   write(u1,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake1(ik,iee))**2
                !   write(u2,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake2(ik,iee))**2
                !   write(u3,'(f12.6,f12.6,f22.6)') d, Energy(iee), abs(Ake3(ik,iee))**2
                !else
                   !write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake(ik,iee))
                   ! A sublattice
                   write(u1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1(ik,iee))
                   write(u2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2(ik,iee))
                   write(u3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3(ik,iee))
                   write(u4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4(ik,iee))
                   ! B sublattice
                   write(uu1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1B(ik,iee))
                   write(uu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2B(ik,iee))
                   write(uu3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3B(ik,iee))
                   write(uu4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4B(ik,iee))
                !endif 
            end do
         end if
      end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent_v2',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalent_v2

subroutine DiagSpectralFunctionKGridInequivalentEnergyCut()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1, uu2, uu3, uu4
   integer :: ik, iee, ie
   integer :: nk(3), ptot,  i1, i2, i3 !, ik, Epts, u, is, uu
   real(dp), pointer :: Kgrid(:,:)=>NULL()
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   real(dp) :: ELoc1(nAt),ELoc2(nAt)
   
   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   !complex(dp), pointer :: PkcLoc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()
 
   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   !integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, rotateOpposite, foldByOne
   logical :: foldByZero

   real(dp) :: K1(3)

   real(dp) :: gridCut
   real(dp) :: K1x1
   real(dp) :: K1x2
   real(dp) :: deltaKx
   real(dp) :: K1y1
   real(dp) :: K1y2
   real(dp) :: deltaKy
   real(dp) :: K1z1
   real(dp) :: K1z2
   real(dp) :: deltaKz
   real(dp) :: topBottomRatio

   logical :: WeiKu, Nishi, useCoordinates, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)

   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       !gcell(1,:) = ucell(1,:)/mmm(1)/2.0
       !gcell(2,:) = ucell(2,:)/mmm(2)/2.0
       !gcell(3,:) = ucell(3,:)
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       call MIO_InputParameter('Spectral.RotateOpposite',rotateOpposite,.false.)
       gcell1 = gcell
       !print*,"gcell before rotation ", gcell(:,1)
       !print*,"gcell before rotation ", gcell(:,2)
       !print*,"gcell before rotation ", gcell(:,3)
       if (rotateRefSystem .eq. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           if (rotateOpposite) then
               phi = phi*180.0_dp/pi
           else
               phi = -phi*180.0_dp/pi
           end if
           aa = phi*pi/180.0_dp
           !print*, "phi =", phi
           !print*, "Angle= ", aa, phi
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           !rot(:,1) = [cos(aa),sin(aa),0.0_dp]
           !rot(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eq. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       !print*, "Angle= ", aa, phi
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if


   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   gcell2 = gcell
   !print*,"gcell after rotation ", gcell(:,1)
   !print*,"gcell after rotation ", gcell(:,2)
   !print*,"gcell after rotation ", gcell(:,3)
   !print*,"grcell ", grcell(:,1)
   !print*,"grcell ", grcell(:,2)
   !print*,"grcell ", grcell(:,3)
   !print*,"ucell ", ucell(:,1)
   !print*,"ucell ", ucell(:,2)
   !print*,"ucell ", ucell(:,3)
   !print*,"rcell ", rcell(:,1)
   !print*,"rcell ", rcell(:,2)
   !print*,"rcell ", rcell(:,3)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)

   call MIO_Print('Calculating Spectral function around K1 (2/3,1/3)','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1]) 
   call MIO_InputParameter('KGridCut',gridCut,0.1_dp)
   !call MIO_InputParameter('Epsilon',eps,0.01_dp)
   !call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   !call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   !call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   !call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   !call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (useCoordinates) then
     if (MIO_InputFindBlock('Spectral.Path',nPath)) then
        call MIO_Allocate(path,[3,nPath],'path','diag')
        call MIO_InputBlock('Spectral.Path',path)
        K1 = path(:,1)
     end if
   else
      K1 = grcell(:,1)*2.0_dp/3.0_dp + grcell(:,2)*1.0_dp/3.0_dp + grcell(:,3)*0.0  ! Macro K valley of graphene
   end if
   print*, "K1: ", K1
   K1x1 = K1(1) - gridCut
   K1x2 = K1(1) + gridCut
   deltaKx = (K1x2 - K1x1)/nk(1)
   K1y1 = K1(2) - gridCut
   K1y2 = K1(2) + gridCut
   deltaKy = (K1y2 - K1y1)/nk(2)
   K1z1 = K1(3) - gridCut
   K1z2 = K1(3) + gridCut
   deltaKz = (K1z2 - K1z1)/nk(3)
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      !Kgrid(:,ik) = grcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
      !              grcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + &
      !              grcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
      Kgrid(1,ik) = (K1x1 + (i1 * deltaKx))
      Kgrid(2,ik) = (K1y1 + (i2 * deltaKy))
      Kgrid(3,ik) = (K1z1 + (i3 * deltaKz))
   end do; end do; end do

   !if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      !call MIO_Print('Spectral function calculation','diag')
      !call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      !call MIO_Allocate(path,[3,nPath],'path','diag')
      !call MIO_Allocate(pathG,[3,nPath],'path','diag')
      !call MIO_InputBlock('Spectral.Path',path)
      !call MIO_InputBlock('Spectral.Path',pathG)
      !do ip=1,nPath
      !   print*, "0: ", path(:,ip)
      !   path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
      !   pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
      !   print*, "1: ", path(:,ip)
      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      !   print*, "2: ", path(:,ip)
      !end do
      !if (nPath==1) then
      !   call MIO_Allocate(nPts,1,'nPts','diag')
      !   nPts(1) = 1
      !   ptsTot = 1
      !else
      !   call MIO_Allocate(nPts,nPath-1,'nPts','diag')
      !   nPts(1) = nPts0
      !   ptsTot = nPts0
      !   if (nPath > 2) then
      !      v = path(:,2) - path(:,1)
      !      d0 = sqrt(dot_product(v,v))
      !      do ip=2,nPath-1
      !         v = path(:,ip+1) - path(:,ip)
      !         d = sqrt(dot_product(v,v))
      !         nPts(ip) = nint(real(d*nPts0)/real(d0))
      !         ptsTot = ptsTot + nPts(ip)
      !      end do
      !   end if
      !end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      !call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      !call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      call MIO_Allocate(Kpts,[3,ptot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptot],'KptsGFrac','diag')
      !Kpts(:,1) = path(:,1)
      !KptsG(:,1) = pathG(:,1)
      KptsG(:,1) = Kgrid(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      !KptsG(:,1) = Kpts(:,1) + GVec 
      !Kpts(:,1) = KptsG(:,1) - GVec 
      !call MIO_InputParameter('CellSize', cellSize, 1)
      !print*, "cell size =", cellSize
      !print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0])
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      !Kpts(1,1) = KptsG(1,1)/cellSize
      !Kpts(2,1) = KptsG(2,1)/cellSize
      !Kpts(3,1) = KptsG(3,1) 
      !kpts(1,1) = mod(KptsG(1,1),abs(G10(1)))
      !kpts(2,1) = mod(KptsG(2,1),abs(G10(2)))
      !kpts(3,1) = KptsG(3,1) 
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      !print*, ll, kk
      !kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
      !kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1) 
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2) 
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1) 
         kpts(2,1) = KptsG(2,1) 
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      end if
 
      !print*, "tup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
      !rcellInv = matinv3(rcell)
      !do i=1,nPath-1
         !do j=1,nPts(i)
         !do j=1,ptot
         do ip=1,ptot
            !ip = ip + 1
            !KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            KptsG(:,ip) = Kgrid(:,ip)
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            !call random_number(randu)
            !randi = FLOOR(22*randu) 
            !call random_number(randu)
            !randj = FLOOR(22*randu) 
            !GVec = matmul(rcell,[randi,randj,0]) ! we only want to translate them by one reciprocal lattice vector
            !KptsG(:,ip) = Kpts(:,ip) + GVec ! Kpts is in SC, Kpts is for Graphene (PC)
            !print*, "yup"
            !print*, Kpts(:,ip), KptsG(:,ip), GVec
            !print*, KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !print*, mod(KptsG(2,ip),G2)
            ! 2 equations, 2 unknowns. Bring point back to SC reciprocal cell
            ! using G10 and G01. 
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            !print*, ll, kk
            !kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
            !kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1) 
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2) 
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip) 
               kpts(2,ip) = KptsG(2,ip) 
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            end if

            !KptsGFrac(2,ip) = grcell(1,1)*(KptsG(2,ip)-KptsG(1,ip)*grcell(2,1)/grcell(1,1))/(grcell(2,2)*grcell(1,1)-grcell(1,2)*grcell(2,1))
            !KptsGFrac(1,ip) = (KptsG(1,ip)-KptsGFrac(2,ip)*grcell(1,2))/grcell(1,1)
            !KptsGFrac(3,ip) = KptsG(3,ip)/grcell(3,3)
            !!!Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)
            !Kpts(:,ip) = KptsGFrac(1,ip)*rcell(:,1) + KptsGFrac(2,ip)*rcell(:,2) + KptsGFrac(3,ip)*rcell(:,3)

            !Kpts(:,ip) = KptsG(:,ip) - GVec ! this one gives the same as when starting from Kpts
            !Kpts(1,ip) = KptsG(1,ip) - FLOOR(KptsG(1,ip)/G1) * G1
            !Kpts(2,ip) = KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !Kpts(3,ip) = KptsG(3,ip) 
            !kpts(1,ip) = mod(kptsG(1,ip),abs(G10(1)))
            !kpts(2,ip) = mod(kptsG(2,ip),abs(G10(2)))
            !kpts(3,ip) = kptsG(3,ip) 
            !Kpts(1,ip) = matmul(rcellInv,KptsG(1,ip))
            !Kpts(2,ip) = matmul(rcellInv,KptsG(2,ip)) 
            !Kpts(3,ip) = matmul(rcellInv,KptsG(3,ip)) 
            !Kpts(1,ip) = KptsG(1,ip)/cellSize
            !Kpts(2,ip) = KptsG(2,ip)/cellSize
            !Kpts(3,ip) = KptsG(3,ip) 
            !v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !d = d + sqrt(dot_product(v,v))
         end do
      !end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectralA'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')

      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptot,nAt,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      !print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !call MIO_Allocate(Ake1Loc,Epts,'Ake1Loc','diag')
      !call MIO_Allocate(Ake2Loc,Epts,'Ake2Loc','diag')
      !Ake1Loc = 0.0_dp
      !Ake2Loc = 0.0_dp
      call MIO_Allocate(Ake1,[ptot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake2,[ptot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake3,[ptot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake4,[ptot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake1B,[ptot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2B,[ptot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake3B,[ptot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4B,[ptot,Epts],'Ake4B','diag')
      Ake1 = 0.0_dp
      Ake2 = 0.0_dp
      Ake3 = 0.0_dp
      Ake4 = 0.0_dp
      Ake1B = 0.0_dp
      Ake2B = 0.0_dp
      Ake3B = 0.0_dp
      Ake4B = 0.0_dp
      !print*, gaussian
      !Pkc = 0.0_dp ! spectral weight PkscI(k)
      !print*, "0",  Pkc
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), REDUCTION(+:Ake1Loc, Ake2Loc), &  
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, gcell1, gcell2, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake) 
      !do ik=1,ptsTot ! K loop
      do ik=1,ptot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                   call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                   call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                             Ake4(ik,iee) = Ake4(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,4))**2
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio 
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio
                                 
                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio 
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                                 !end if
                             end if
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                             Ake4(ik,iee) = Ake4(ik,iee) + PkcLocA(ie,4)
                          end if
                       end if
                end do
            end do
         !end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:) + Ake4(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      !print*, "lets write it all out"
      do ik=1,ptot ! K loop
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                !write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake(ik,iee))
                write(u1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1(ik,iee))
                write(u2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2(ik,iee))
                write(u3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3(ik,iee))
                write(u4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4(ik,iee))
                write(uu1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1B(ik,iee))
                write(uu2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2B(ik,iee))
                write(uu3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3B(ik,iee))
                write(uu4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4B(ik,iee))
            end do
         end if
      end do
      !end do; end do; end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Kgrid,'Kgrid','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   !end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalentEnergyCut

subroutine DiagSpectralFunctionKGridInequivalentEnergyCut_v2()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1, uu2, uu3, uu4
   integer :: ik, iee, ie
   integer :: nk(3), ptot,  i1, i2, i3 !, ik, Epts, u, is, uu
   real(dp), pointer :: Kgrid(:,:)=>NULL()
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   real(dp) :: ELoc1(nAt),ELoc2(nAt)
   
   logical :: MoireBS, GaussConv, DirectBand
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   !complex(dp), pointer :: PkcLoc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: EAke(:,:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()
 
   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   !integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, rotateOpposite, foldByOne
   logical :: foldByZero

   real(dp) :: K1(3)

   real(dp) :: gridCut
   real(dp) :: K1x1
   real(dp) :: K1x2
   real(dp) :: deltaKx
   real(dp) :: K1y1
   real(dp) :: K1y2
   real(dp) :: deltaKy
   real(dp) :: K1z1
   real(dp) :: K1z2
   real(dp) :: deltaKz
   real(dp) :: topBottomRatio

   logical :: WeiKu, Nishi, useCoordinates, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut_v2',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   print*, rcell(:,1)
   print*, rcell(:,2)
   print*, rcell(:,3)

   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       !gcell(1,:) = ucell(1,:)/mmm(1)/2.0
       !gcell(2,:) = ucell(2,:)/mmm(2)/2.0
       !gcell(3,:) = ucell(3,:)
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       call MIO_InputParameter('Spectral.RotateOpposite',rotateOpposite,.false.)
       gcell1 = gcell
       print*,"gcell before rotation ", gcell(:,1)
       print*,"gcell before rotation ", gcell(:,2)
       print*,"gcell before rotation ", gcell(:,3)
       if (rotateRefSystem .eq. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           if (rotateOpposite) then
               phi = phi*180.0_dp/pi
           else
               phi = -phi*180.0_dp/pi
           end if
           aa = phi*pi/180.0_dp
           print*, "phi =", phi
           !print*, "Angle= ", aa, phi
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           !rot(:,1) = [cos(aa),sin(aa),0.0_dp]
           !rot(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eq. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       !print*, "Angle= ", aa, phi
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if


   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   gcell2 = gcell
   print*,"gcell after rotation ", gcell(:,1)
   print*,"gcell after rotation ", gcell(:,2)
   print*,"gcell after rotation ", gcell(:,3)
   print*,"grcell ", grcell(:,1)
   print*,"grcell ", grcell(:,2)
   print*,"grcell ", grcell(:,3)
   print*,"ucell ", ucell(:,1)
   print*,"ucell ", ucell(:,2)
   print*,"ucell ", ucell(:,3)
   print*,"rcell ", rcell(:,1)
   print*,"rcell ", rcell(:,2)
   print*,"rcell ", rcell(:,3)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)

   call MIO_Print('Calculating Spectral function around K1 (2/3,1/3)','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1]) 
   call MIO_InputParameter('KGridCut',gridCut,0.1_dp)
   !call MIO_InputParameter('Epsilon',eps,0.01_dp)
   !call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   !call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   !call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   !call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   !call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (useCoordinates) then
     if (MIO_InputFindBlock('Spectral.Path',nPath)) then
        call MIO_Allocate(path,[3,nPath],'path','diag')
        call MIO_InputBlock('Spectral.Path',path)
        K1 = path(:,1)
     end if
   else
      K1 = grcell(:,1)*2.0_dp/3.0_dp + grcell(:,2)*1.0_dp/3.0_dp + grcell(:,3)*0.0  ! Macro K valley of graphene
   end if
   print*, "K1: ", K1
   K1x1 = K1(1) - gridCut
   K1x2 = K1(1) + gridCut
   deltaKx = (K1x2 - K1x1)/nk(1)
   K1y1 = K1(2) - gridCut
   K1y2 = K1(2) + gridCut
   deltaKy = (K1y2 - K1y1)/nk(2)
   K1z1 = K1(3) - gridCut
   K1z2 = K1(3) + gridCut
   deltaKz = (K1z2 - K1z1)/nk(3)
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      !Kgrid(:,ik) = grcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
      !              grcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + &
      !              grcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
      Kgrid(1,ik) = (K1x1 + (i1 * deltaKx))
      Kgrid(2,ik) = (K1y1 + (i2 * deltaKy))
      Kgrid(3,ik) = (K1z1 + (i3 * deltaKz))
   end do; end do; end do

   !if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      !call MIO_Print('Spectral function calculation','diag')
      !call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      !call MIO_Allocate(path,[3,nPath],'path','diag')
      !call MIO_Allocate(pathG,[3,nPath],'path','diag')
      !call MIO_InputBlock('Spectral.Path',path)
      !call MIO_InputBlock('Spectral.Path',pathG)
      !do ip=1,nPath
      !   print*, "0: ", path(:,ip)
      !   path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
      !   pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
      !   print*, "1: ", path(:,ip)
      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      !   print*, "2: ", path(:,ip)
      !end do
      !if (nPath==1) then
      !   call MIO_Allocate(nPts,1,'nPts','diag')
      !   nPts(1) = 1
      !   ptsTot = 1
      !else
      !   call MIO_Allocate(nPts,nPath-1,'nPts','diag')
      !   nPts(1) = nPts0
      !   ptsTot = nPts0
      !   if (nPath > 2) then
      !      v = path(:,2) - path(:,1)
      !      d0 = sqrt(dot_product(v,v))
      !      do ip=2,nPath-1
      !         v = path(:,ip+1) - path(:,ip)
      !         d = sqrt(dot_product(v,v))
      !         nPts(ip) = nint(real(d*nPts0)/real(d0))
      !         ptsTot = ptsTot + nPts(ip)
      !      end do
      !   end if
      !end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      !call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      !call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      call MIO_Allocate(Kpts,[3,ptot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptot],'KptsGFrac','diag')
      !Kpts(:,1) = path(:,1)
      !KptsG(:,1) = pathG(:,1)
      KptsG(:,1) = Kgrid(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      !KptsG(:,1) = Kpts(:,1) + GVec 
      !Kpts(:,1) = KptsG(:,1) - GVec 
      !call MIO_InputParameter('CellSize', cellSize, 1)
      !print*, "cell size =", cellSize
      print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0])
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      !Kpts(1,1) = KptsG(1,1)/cellSize
      !Kpts(2,1) = KptsG(2,1)/cellSize
      !Kpts(3,1) = KptsG(3,1) 
      !kpts(1,1) = mod(KptsG(1,1),abs(G10(1)))
      !kpts(2,1) = mod(KptsG(2,1),abs(G10(2)))
      !kpts(3,1) = KptsG(3,1) 
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      !print*, ll, kk
      !kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
      !kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1) 
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2) 
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1) 
         kpts(2,1) = KptsG(2,1) 
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      end if
 
      !print*, "tup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
      !rcellInv = matinv3(rcell)
      !do i=1,nPath-1
         !do j=1,nPts(i)
         !do j=1,ptot
         do ip=1,ptot
            !ip = ip + 1
            !KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            KptsG(:,ip) = Kgrid(:,ip)
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            !call random_number(randu)
            !randi = FLOOR(22*randu) 
            !call random_number(randu)
            !randj = FLOOR(22*randu) 
            !GVec = matmul(rcell,[randi,randj,0]) ! we only want to translate them by one reciprocal lattice vector
            !KptsG(:,ip) = Kpts(:,ip) + GVec ! Kpts is in SC, Kpts is for Graphene (PC)
            !print*, "yup"
            !print*, Kpts(:,ip), KptsG(:,ip), GVec
            !print*, KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !print*, mod(KptsG(2,ip),G2)
            ! 2 equations, 2 unknowns. Bring point back to SC reciprocal cell
            ! using G10 and G01. 
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            !print*, ll, kk
            !kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
            !kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1) 
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2) 
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip) 
               kpts(2,ip) = KptsG(2,ip) 
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            end if

            !KptsGFrac(2,ip) = grcell(1,1)*(KptsG(2,ip)-KptsG(1,ip)*grcell(2,1)/grcell(1,1))/(grcell(2,2)*grcell(1,1)-grcell(1,2)*grcell(2,1))
            !KptsGFrac(1,ip) = (KptsG(1,ip)-KptsGFrac(2,ip)*grcell(1,2))/grcell(1,1)
            !KptsGFrac(3,ip) = KptsG(3,ip)/grcell(3,3)
            !!!Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)
            !Kpts(:,ip) = KptsGFrac(1,ip)*rcell(:,1) + KptsGFrac(2,ip)*rcell(:,2) + KptsGFrac(3,ip)*rcell(:,3)

            !Kpts(:,ip) = KptsG(:,ip) - GVec ! this one gives the same as when starting from Kpts
            !Kpts(1,ip) = KptsG(1,ip) - FLOOR(KptsG(1,ip)/G1) * G1
            !Kpts(2,ip) = KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !Kpts(3,ip) = KptsG(3,ip) 
            !kpts(1,ip) = mod(kptsG(1,ip),abs(G10(1)))
            !kpts(2,ip) = mod(kptsG(2,ip),abs(G10(2)))
            !kpts(3,ip) = kptsG(3,ip) 
            !Kpts(1,ip) = matmul(rcellInv,KptsG(1,ip))
            !Kpts(2,ip) = matmul(rcellInv,KptsG(2,ip)) 
            !Kpts(3,ip) = matmul(rcellInv,KptsG(3,ip)) 
            !Kpts(1,ip) = KptsG(1,ip)/cellSize
            !Kpts(2,ip) = KptsG(2,ip)/cellSize
            !Kpts(3,ip) = KptsG(3,ip) 
            !v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !d = d + sqrt(dot_product(v,v))
         end do
      !end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral_v2'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')

      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptot,nAt,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      call MIO_InputParameter('Spectral.DirectBandVal',DirectBand,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !call MIO_Allocate(Ake1Loc,Epts,'Ake1Loc','diag')
      !call MIO_Allocate(Ake2Loc,Epts,'Ake2Loc','diag')
      !Ake1Loc = 0.0_dp
      !Ake2Loc = 0.0_dp
      call MIO_Allocate(Ake1,[ptot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake2,[ptot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake3,[ptot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake4,[ptot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake1B,[ptot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2B,[ptot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake3B,[ptot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4B,[ptot,Epts],'Ake4B','diag')
      Ake1 = 0.0_dp
      Ake2 = 0.0_dp
      Ake3 = 0.0_dp
      Ake4 = 0.0_dp
      Ake1B = 0.0_dp
      Ake2B = 0.0_dp
      Ake3B = 0.0_dp
      Ake4B = 0.0_dp
      
      if (DirectBand) then
        call MIO_Allocate(EAke,[ptot,Epts],'EAke','diag')
        EAke = 0.0_dp
      end if

      !print*, gaussian
      !Pkc = 0.0_dp ! spectral weight PkscI(k)
      !print*, "0",  Pkc
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), REDUCTION(+:Ake1Loc, Ake2Loc), &  
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, gcell1, gcell2, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake) 
      !do ik=1,ptsTot ! K loop
      do ik=1,ptot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                   call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                   call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                             Ake4(ik,iee) = Ake4(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,4))**2
                          else if (DirectBand) then
                             if (iee == int(ie-(nAt/2-Epts/2))) then
                                 EAke(ik,iee) = ELoc(ie)
                                 
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio 
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio
                                 
                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio 
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                             end if
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio 
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio
                                 
                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio 
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                                 !end if
                             end if
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                             Ake4(ik,iee) = Ake4(ik,iee) + PkcLocA(ie,4)
                          end if
                       end if
                end do
            end do
         !end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      print*, "lets write it all out"
      do ik=1,ptot ! K loop
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else if (DirectBand) then
            do iee=1,Epts
                write(u,'(3f12.6,f12.6,10f12.6)') KptsG(:,ik), REAL(EAke(ik,iee)), REAL(Ake1(ik,iee)), REAL(Ake1B(ik,iee)),&
                                                                                 & REAL(Ake2(ik,iee)), REAL(Ake2B(ik,iee)),&
                                                                                 & REAL(Ake3(ik,iee)), REAL(Ake3B(ik,iee)),&
                                                                                 & REAL(Ake4(ik,iee)), REAL(Ake4B(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                !write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake(ik,iee))
                write(u1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1(ik,iee))
                write(u2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2(ik,iee))
                write(u3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3(ik,iee))
                write(u4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4(ik,iee))
                write(uu1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1B(ik,iee))
                write(uu2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2B(ik,iee))
                write(uu3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3B(ik,iee))
                write(uu4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4B(ik,iee))
            end do
         end if
      end do
      !end do; end do; end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Kgrid,'Kgrid','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   !end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut_v2',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalentEnergyCut_v2

subroutine DiagSpectralFunctionKGridInequivalentEnergyCutNickDale()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3
   integer :: ik, iee, ie
   integer :: nk(3), ptot,  i1, i2, i3 !, ik, Epts, u, is, uu
   real(dp), pointer :: Kgrid(:,:)=>NULL()
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   real(dp) :: ELoc1(nAt),ELoc2(nAt)
   
   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   !complex(dp), pointer :: PkcLoc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()
 
   complex(dp) :: PkcLocA(nAt,3)
   complex(dp) :: PkcLocB(nAt,3)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   !integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, rotateOpposite, foldByOne
   logical :: foldByZero, lowerGridHalf

   real(dp) :: K1(3)

   real(dp) :: gridCut
   real(dp) :: gridCutX
   real(dp) :: gridCutY
   real(dp) :: K1x1
   real(dp) :: K1x2
   real(dp) :: deltaKx
   real(dp) :: K1y1
   real(dp) :: K1y2
   real(dp) :: deltaKy
   real(dp) :: K1z1
   real(dp) :: K1z2
   real(dp) :: deltaKz
   real(dp) :: topBottomRatio

   logical :: WeiKu, Nishi, useCoordinates, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   !print*, rcell(:,1)
   !print*, rcell(:,2)
   !print*, rcell(:,3)

   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       !gcell(1,:) = ucell(1,:)/mmm(1)/2.0
       !gcell(2,:) = ucell(2,:)/mmm(2)/2.0
       !gcell(3,:) = ucell(3,:)
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       call MIO_InputParameter('Spectral.RotateOpposite',rotateOpposite,.false.)
       gcell1 = gcell
       !print*,"gcell before rotation ", gcell(:,1)
       !print*,"gcell before rotation ", gcell(:,2)
       !print*,"gcell before rotation ", gcell(:,3)
       if (rotateRefSystem .eq. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           if (rotateOpposite) then
               phi = phi*180.0_dp/pi
           else
               phi = -phi*180.0_dp/pi
           end if
           aa = phi*pi/180.0_dp
           !print*, "phi =", phi
           !print*, "Angle= ", aa, phi
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           !rot(:,1) = [cos(aa),sin(aa),0.0_dp]
           !rot(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eq. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       !print*, "Angle= ", aa, phi
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if


   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   gcell2 = gcell
   !print*,"gcell after rotation ", gcell(:,1)
   !print*,"gcell after rotation ", gcell(:,2)
   !print*,"gcell after rotation ", gcell(:,3)
   !print*,"grcell ", grcell(:,1)
   !print*,"grcell ", grcell(:,2)
   !print*,"grcell ", grcell(:,3)
   !print*,"ucell ", ucell(:,1)
   !print*,"ucell ", ucell(:,2)
   !print*,"ucell ", ucell(:,3)
   !print*,"rcell ", rcell(:,1)
   !print*,"rcell ", rcell(:,2)
   !print*,"rcell ", rcell(:,3)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)

   call MIO_Print('Calculating Spectral function around K1 (2/3,1/3)','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1]) 
   call MIO_InputParameter('KGridCut',gridCut,0.1_dp)
   call MIO_InputParameter('KGridCutX',gridCutX,0.1_dp)
   call MIO_InputParameter('KGridCutY',gridCutY,0.1_dp)
   call MIO_InputParameter('KGridLowerGridHalf',lowerGridHalf,.false.)
   !call MIO_InputParameter('Epsilon',eps,0.01_dp)
   !call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   !call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   !call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   !call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   !call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (useCoordinates) then
     if (MIO_InputFindBlock('Spectral.Path',nPath)) then
        call MIO_Allocate(path,[3,nPath],'path','diag')
        call MIO_InputBlock('Spectral.Path',path)
        K1 = path(:,1)
     end if
   else
      K1 = grcell(:,1)*2.0_dp/3.0_dp + grcell(:,2)*1.0_dp/3.0_dp + grcell(:,3)*0.0
   end if
   print*, "K1: ", K1
   if (lowerGridHalf) then
      K1x1 = K1(1) - gridCutX
      K1x2 = K1(1) 
   else
      K1x1 = K1(1) 
      K1x2 = K1(1) + gridCutX
   end if
   deltaKx = (K1x2 - K1x1)/nk(1)
   K1y1 = K1(2) + gridCutY
   !K1y2 = K1(2) + gridCutX
   !deltaKy = (K1y2 - K1y1)/nk(2)
   K1z1 = K1(3) - gridCut
   K1z2 = K1(3) + gridCut
   deltaKz = (K1z2 - K1z1)/nk(3)
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      !Kgrid(:,ik) = grcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
      !              grcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + &
      !              grcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
      Kgrid(1,ik) = (K1x1 + (i1 * deltaKx))
      !Kgrid(2,ik) = (K1y1 + (i2 * deltaKy))
      Kgrid(2,ik) = K1y1 
      Kgrid(3,ik) = (K1z1 + (i3 * deltaKz))
   end do; end do; end do

   !if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      !call MIO_Print('Spectral function calculation','diag')
      !call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      !call MIO_Allocate(path,[3,nPath],'path','diag')
      !call MIO_Allocate(pathG,[3,nPath],'path','diag')
      !call MIO_InputBlock('Spectral.Path',path)
      !call MIO_InputBlock('Spectral.Path',pathG)
      !do ip=1,nPath
      !   print*, "0: ", path(:,ip)
      !   path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
      !   pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
      !   print*, "1: ", path(:,ip)
      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      !   print*, "2: ", path(:,ip)
      !end do
      !if (nPath==1) then
      !   call MIO_Allocate(nPts,1,'nPts','diag')
      !   nPts(1) = 1
      !   ptsTot = 1
      !else
      !   call MIO_Allocate(nPts,nPath-1,'nPts','diag')
      !   nPts(1) = nPts0
      !   ptsTot = nPts0
      !   if (nPath > 2) then
      !      v = path(:,2) - path(:,1)
      !      d0 = sqrt(dot_product(v,v))
      !      do ip=2,nPath-1
      !         v = path(:,ip+1) - path(:,ip)
      !         d = sqrt(dot_product(v,v))
      !         nPts(ip) = nint(real(d*nPts0)/real(d0))
      !         ptsTot = ptsTot + nPts(ip)
      !      end do
      !   end if
      !end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      !call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      !call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      call MIO_Allocate(Kpts,[3,ptot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptot],'KptsGFrac','diag')
      !Kpts(:,1) = path(:,1)
      !KptsG(:,1) = pathG(:,1)
      KptsG(:,1) = Kgrid(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      !KptsG(:,1) = Kpts(:,1) + GVec 
      !Kpts(:,1) = KptsG(:,1) - GVec 
      !call MIO_InputParameter('CellSize', cellSize, 1)
      !print*, "cell size =", cellSize
      !print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0])
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      !Kpts(1,1) = KptsG(1,1)/cellSize
      !Kpts(2,1) = KptsG(2,1)/cellSize
      !Kpts(3,1) = KptsG(3,1) 
      !kpts(1,1) = mod(KptsG(1,1),abs(G10(1)))
      !kpts(2,1) = mod(KptsG(2,1),abs(G10(2)))
      !kpts(3,1) = KptsG(3,1) 
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      !print*, ll, kk
      !kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
      !kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1) 
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2) 
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1) 
         kpts(2,1) = KptsG(2,1) 
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
      end if
 
      !print*, "tup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
      !rcellInv = matinv3(rcell)
      !do i=1,nPath-1
         !do j=1,nPts(i)
         !do j=1,ptot
         do ip=1,ptot
            !ip = ip + 1
            !KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            KptsG(:,ip) = Kgrid(:,ip)
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            !call random_number(randu)
            !randi = FLOOR(22*randu) 
            !call random_number(randu)
            !randj = FLOOR(22*randu) 
            !GVec = matmul(rcell,[randi,randj,0]) ! we only want to translate them by one reciprocal lattice vector
            !KptsG(:,ip) = Kpts(:,ip) + GVec ! Kpts is in SC, Kpts is for Graphene (PC)
            !print*, "yup"
            !print*, Kpts(:,ip), KptsG(:,ip), GVec
            !print*, KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !print*, mod(KptsG(2,ip),G2)
            ! 2 equations, 2 unknowns. Bring point back to SC reciprocal cell
            ! using G10 and G01. 
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            !print*, ll, kk
            !kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
            !kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1) 
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2) 
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip) 
               kpts(2,ip) = KptsG(2,ip) 
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 
            end if

            !KptsGFrac(2,ip) = grcell(1,1)*(KptsG(2,ip)-KptsG(1,ip)*grcell(2,1)/grcell(1,1))/(grcell(2,2)*grcell(1,1)-grcell(1,2)*grcell(2,1))
            !KptsGFrac(1,ip) = (KptsG(1,ip)-KptsGFrac(2,ip)*grcell(1,2))/grcell(1,1)
            !KptsGFrac(3,ip) = KptsG(3,ip)/grcell(3,3)
            !!!Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)
            !Kpts(:,ip) = KptsGFrac(1,ip)*rcell(:,1) + KptsGFrac(2,ip)*rcell(:,2) + KptsGFrac(3,ip)*rcell(:,3)

            !Kpts(:,ip) = KptsG(:,ip) - GVec ! this one gives the same as when starting from Kpts
            !Kpts(1,ip) = KptsG(1,ip) - FLOOR(KptsG(1,ip)/G1) * G1
            !Kpts(2,ip) = KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !Kpts(3,ip) = KptsG(3,ip) 
            !kpts(1,ip) = mod(kptsG(1,ip),abs(G10(1)))
            !kpts(2,ip) = mod(kptsG(2,ip),abs(G10(2)))
            !kpts(3,ip) = kptsG(3,ip) 
            !Kpts(1,ip) = matmul(rcellInv,KptsG(1,ip))
            !Kpts(2,ip) = matmul(rcellInv,KptsG(2,ip)) 
            !Kpts(3,ip) = matmul(rcellInv,KptsG(3,ip)) 
            !Kpts(1,ip) = KptsG(1,ip)/cellSize
            !Kpts(2,ip) = KptsG(2,ip)/cellSize
            !Kpts(3,ip) = KptsG(3,ip) 
            !v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !d = d + sqrt(dot_product(v,v))
         end do
      !end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptot
      flnm = trim(prefix)//'.spectral1'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptot,nAt,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      !print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !call MIO_Allocate(Ake1Loc,Epts,'Ake1Loc','diag')
      !call MIO_Allocate(Ake2Loc,Epts,'Ake2Loc','diag')
      !Ake1Loc = 0.0_dp
      !Ake2Loc = 0.0_dp
      call MIO_Allocate(Ake1,[ptot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake2,[ptot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake3,[ptot,Epts],'Ake3','diag')
      Ake1 = 0.0_dp
      Ake2 = 0.0_dp
      Ake3 = 0.0_dp
      !print*, gaussian
      !Pkc = 0.0_dp ! spectral weight PkscI(k)
      !print*, "0",  Pkc
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), REDUCTION(+:Ake1Loc, Ake2Loc), &  
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, gcell1, gcell2, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake) 
      !do ik=1,ptsTot ! K loop
      do ik=1,ptot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                   call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                   call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio 
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 !end if
                             end if
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                          end if
                       end if
                end do
            end do
         !end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      print*, "lets write it all out"
      do ik=1,ptot ! K loop
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake(ik,iee))
                write(u1,'(f12.6,f12.6,f22.6)') KptsG(:,ik), Energy(iee), REAL(Ake1(ik,iee))
                write(u2,'(f12.6,f12.6,f22.6)') KptsG(:,ik), Energy(iee), REAL(Ake2(ik,iee))
                write(u3,'(f12.6,f12.6,f22.6)') KptsG(:,ik), Energy(iee), REAL(Ake3(ik,iee))
            end do
         end if
      end do
      !end do; end do; end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Kgrid,'Kgrid','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   !end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalentEnergyCutNickDale

subroutine DiagHam(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   !complex(dp), intent(out) :: H(N,N)
   !real(dp), intent(out) :: E(N)
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   !print*, cell
   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   !call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      !print*, "info =", info
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if

end subroutine DiagHam

subroutine DiagHamSparse(N, ns, is, ELoc, KLoc, cell, H0, maxN, hopp, NList, Nneigh, neighCell,neig)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species
    use tbpar, only : U
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is, neig
    real(dp), intent(out) :: ELoc(N)
    real(dp), intent(in) :: KLoc(3), cell(3,3)
    real(dp), intent(inout) :: H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    real(dp) :: tol

    integer :: i, j, in, info
    real(dp) :: R(3), zz, resid_norm

    ! ARPACK parameters and variables
    integer :: nev, ncv, lworkl, ido, ierr
    character(1) :: bmat
    character(2) :: which
    complex(dp), allocatable :: resid(:), v(:,:), workd(:), d(:), workev(:)
    complex(dp), allocatable :: workl(:)
    complex(dp), allocatable :: z(:,:)
    double precision, allocatable :: rwork(:)
    integer, allocatable :: iparam(:), ipntr(:)
    logical, allocatable :: select(:)
    integer :: max_iter, nn, iter
    logical :: useShift
    double precision :: shift

    integer, allocatable :: row_ptr(:), col_ind(:)
    complex(dp), allocatable :: values(:)
    real(dp), allocatable :: rand_real(:), rand_imag(:)
    complex(dp) :: sigma

    ! PARDISO variables
    !integer, parameter :: mtype = 13  ! Complex unsymmetric matrix
    !integer :: pt(64), iparm(64), maxfct, mnum, phase, error, msglvl
    !integer :: nnn, nrhs
    !complex(dp), allocatable :: a(:)
    !integer, allocatable :: ia(:), ja(:)
    !complex(dp) :: ddum
    !integer :: idum

    INTEGER :: nnn, nnz
    !complex(dp) :: x(N) ! Pardiso needs X(N) even if sol is returned in b

    !complex(dp) :: b(N)
!C.. Internal solver memory pointer for 64-bit architectures
    INTEGER*8 pt(64)
!C.. Internal solver memory pointer for 32-bit architectures
!C.. INTEGER*4 pt(64)

!C.. This is OK in both cases
!    TYPE(MKL_PARDISO_HANDLE) pt(64)
!C.. All other variables
    INTEGER maxfct, mnum, mtype, phase, nrhs, msglvl
    INTEGER iparm(64)
    !allocate(row_ptr(N+1), col_ind(nnz_temp), values(nnz_temp))
    INTEGER, allocatable :: ia(:) ! row_ptr
    INTEGER, allocatable :: ja(:) ! col_ind
    INTEGER idum(1)
    complex(dp) ::  ddum(1)
    integer error
    complex(dp), allocatable :: a(:) ! values
    
    DATA nrhs /1/, maxfct /1/, mnum /1/

    complex(dp), allocatable :: EVectors(:,:)
    logical :: saveRitz, symmetric
    integer :: nconv
    complex(dp), allocatable :: ax(:)
    double precision, allocatable :: rd(:,:)
    double precision :: dznrm2, dlapy2


    
    ! Parameters
    !integer, parameter :: dp = c_double
    !integer, parameter :: int_kind = c_int
    
    ! Variables
    !type(C_PTR) :: mkl_handle
    !!integer(int_kind) :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !!integer(int_kind), allocatable :: row_ptr(:), col_ind(:)
    !integer :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !integer, allocatable :: row_ptr(:), col_ind(:)
    !real(dp), allocatable :: values(:), b(:), x(:), temp_value
    !type(SPARSE_MATRIX_DESCR) :: descr


    !integer :: nev, ncv, lworkl 
 
    nev=neig 
    ncv=nev*10
    lworkl= 3*NCV**2 + 5*NCV

    max_iter = 10000
    allocate(resid(N), v(N, ncv), workd(3*N), workl(lworkl), rwork(ncv), d(nev+1), iparam(11), ipntr(14), select(ncv), rand_real(N), rand_imag(N))
    allocate(z(N,nev))
    allocate(workev(2*ncv))
    allocate(ax(N))
    allocate(rd(ncv,3))

    ! Ensure correct size of nev and ncv
    if ( (nev < 1) .or. (nev >= ncv) .or. (ncv > N) ) then
        print *, 'Error: invalid parameters - nev=', nev, 'ncv=', ncv, 'N=', N
        stop
    end if

    ! Initialize arrays
    !v = 0.0d0
    !workd = 0.0d0
    !workl = 0.0d0
    !rwork = 0.0d0
    !d = 0.0d0

    bmat = 'I'

    iparam(1) = 1
    iparam(3) = max_iter
    call MIO_InputParameter('Diag.SparseSetShift',shift,0.0_dp)
    call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)
    call MIO_InputParameter('Diag.SparseSaveRitz',saveRitz,.false.)
    if (useShift) then
       which = 'LM'
       iparam(7) = 3
       sigma = cmplx(shift, 0.0_dp)
       print*, "We will find eigenvalues close to ", shift
    else
       which = 'SM'
       sigma = (0.0_dp, 0.0_dp)
       iparam(7) = 1
    end if
    ! Initialize ARPACK parameters

    call MIO_InputParameter('Diag.SparseSetTol',tol,0.1_dp)

    ido = 0
    info = 0

    !iparam(4) = 1
    !iparam(7) = 1

    nn = N

    ! Initialize the starting vector resid with random values
    !call random_number(rand_real)
    !call random_number(rand_imag)
    !resid = cmplx(rand_real, rand_imag)
    !resid = 0

    ! Check resid for initial state
    if (size(resid) /= N) then
        print *, 'Error: resid size mismatch: ', size(resid), ' expected: ', N
        stop
    end if

     ! Debug print for resid
    if (any(resid /= resid)) then
        print *, 'Error: resid contains NaN values initially.'
        stop
    end if

    resid_norm = sqrt(sum(abs(resid)**2))

    ! Debug prints
    !print *, 'Initial resid norm: ', resid_norm
    !print *, 'Initial iparam: ', iparam
    !print *, 'nn, nev, ncv: ', nn, nev, ncv

    ! Create CSR sparse matrix storage
    print *, "initialize the sparse matrix and put it in csr format"
    call initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values,sigma)
    print *, "done"
    !call test_sparse_matvec(nn, row_ptr, col_ind, values)
    symmetric = is_structurally_symmetric(values, row_ptr, col_ind, n)

    print*, "is it symmetric?", symmetric

    ! Debug prints for CSR matrix
    if (any(values /= values)) then
        print *, 'Error: values contains NaN values after initialization.'
        stop
    end if

    
    if (maxval(abs(values)) > 1e10) then
        print *, 'Warning: values contains extremely large values.'
    end if

    !if (any(col_ind < 1 .or. col_ind > N)) then
    !    print *, 'Error: col_ind contains invalid indices after initialization.'
    !    stop
    !end if



    !print *, 'iparam: ', iparam
    !print *, 'ipntr: ', ipntr
    !print *, 'ido: ', ido
    !print *, 'bmat: ', bmat
    !print *, 'nn: ', nn
    !print *, 'which: ', which
    !print *, 'nev: ', nev
    !print *, 'tol: ', tol
    call znaupd(ido, bmat, nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)

    resid_norm = sqrt(sum(abs(resid)**2))

    ! Debug prints
    !print *, 'Resid norm after znaupd: ', resid_norm

    ! Check for convergence and errors
    if (info == -5) then
        print *, 'Maximum number of iterations reached, INFO = ', info
        print *, 'Number of converged eigenvalues:', iparam(5)
        print *, 'iparam: ', iparam
        print *, 'ipntr: ', ipntr
        print *, 'ido: ', ido
        print *, 'bmat: ', bmat
        print *, 'nn: ', nn
        print *, 'which: ', which
        print *, 'nev: ', nev
        print *, 'tol: ', tol
    else if (info /= 0) then
        print *, 'Error with znaupd, INFO = ', info
        print *, 'iparam: ', iparam
        print *, 'ipntr: ', ipntr
        stop
    else
        print *, 'znaupd converged successfully'
    endif
 

    !if (useShift) then
    !    ! Sort each row
    !    do i = 1, N
    !        row_start = row_ptr(i)
    !        row_end = row_ptr(i + 1) - 1
    !        do j = row_start, row_end - 1
    !            do k = j + 1, row_end
    !                if (col_ind(j) > col_ind(k)) then
    !                    ! Swap indices
    !                    temp_index = col_ind(j)
    !                    col_ind(j) = col_ind(k)
    !                    col_ind(k) = temp_index
    !                    ! Swap values
    !                    temp_value = values(j)
    !                    values(j) = values(k)
    !                    values(k) = temp_value
    !                end if
    !            end do
    !        end do
    !    end do
    !    ! Set matrix descriptor
    !    descr.type = SPARSE_MATRIX_TYPE_GENERAL
    !    descr.mode = SPARSE_FILL_MODE_LOWER
    !    descr.diag = SPARSE_DIAG_NON_UNIT

    !    ! Create the matrix handle
    !    status = mkl_sparse_d_create_csr(mkl_handle, SPARSE_INDEX_BASE_ONE, N, N, row_ptr, row_ptr(2:N+1), col_ind, values)
    !    if (status /= SPARSE_STATUS_SUCCESS) then
    !        print *, "Error in matrix creation"
    !        stop
    !    end if

    !    ! Perform the analysis phase
    !    status = mkl_sparse_set_mv_hint(mkl_handle, SPARSE_OPERATION_NON_TRANSPOSE, descr, 1000)
    !    if (status /= SPARSE_STATUS_SUCCESS) then
    !        print *, "Error in setting hints"
    !        stop
    !    end if

    !    ! Optimize the matrix structure
    !    status = mkl_sparse_optimize(mkl_handle)
    !    if (status /= SPARSE_STATUS_SUCCESS) then
    !        print *, "Error in optimizing the matrix"
    !        stop
    !    end if
    !!   ! PARDISO initialization
    !!       nnn = nn
    !!       !nrhs = 1
    !!       a = values
    !!       ia = row_ptr
    !!       ja = col_ind
    !!       maxfct = 1
    !!       mnum = 1
    !!       msglvl = 0
    !!       error = 0
    !!   
    !!       !call pardisoinit(pt, mtype, iparm)
    !!       DO i = 1, 64
    !!          iparm(i) = 0
    !!          pt(i) = 0
    !!       END DO
    !!       iparm(1) = 1 ! no solver default
    !!       iparm(2) = 0 ! fill-in reordering from METIS
    !!       iparm(3) = 1 ! numbers of processors
    !!       iparm(4) = 0 ! no iterative-direct algorithm
    !!       iparm(5) = 0 ! no user fill-in reducing permutation
    !!       iparm(6) = 0 ! =0 solution on the first n components of x
    !!       iparm(7) = 0 ! not in use
    !!       iparm(8) = 9 ! numbers of iterative refinement steps
    !!       iparm(9) = 0 ! not in use
    !!       iparm(10) = 13 ! perturb the pivot elements with 1E-13
    !!       !iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS ! try changing this
    !!       iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS ! try changing this
    !!       iparm(12) = 0 ! not in use
    !!       !iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric) ! try changin this
    !!       iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric) ! try changin this
    !!       iparm(14) = 0 ! Output: number of perturbed pivots
    !!       iparm(15) = 0 ! not in use
    !!       iparm(16) = 0 ! not in use
    !!       iparm(17) = 0 ! not in use
    !!       iparm(18) = -1 ! Output: number of nonzeros in the factor LU
    !!       iparm(19) = -1 ! Output: Mflops for LU factorization
    !!       iparm(20) = 0 ! Output: Numbers of CG Iterations
    !!       error = 0 ! initialize error flag
    !!       msglvl = 0 ! print statistical information
    !!       !mtype = 11 ! real unsymmetric
    !!       !mtype = 13 ! complex and structurally nonsymmetric
    !!       mtype = 6 ! complex and structurally symmetric
 
    !!       !DO i = 1, 64
    !!       !  pt(i)%DUMMY = 0
    !!       !END DO
    !!   
    !!       phase = 11  ! Reordering and Symbolic Factorization
    !!       call pardiso(pt, maxfct, mnum, mtype, phase, nnn, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    !!       if (error /= 0) then
    !!           print *, 'PARDISO error during symbolic factorization:', error
    !!           stop
    !!       end if
    !!       print*, 'Reordering completed ... '
    !!   
    !!       phase = 22  ! Numerical factorization
    !!       call pardiso(pt, maxfct, mnum, mtype, phase, nnn, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    !!       if (error /= 0) then
    !!           print *, 'PARDISO error during numerical factorization:', error
    !!           stop
    !!       end if
    !!       print*, 'Factorization completed ... '

    !end if


    iter = 0

    if (useShift) then
       do while (ido /= 99)
           if (ido == -1 .or. ido == 1) then
               print*, "add the mkl solver here"
               !status = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0_dp, mkl_handle, descr, b, x)
               !if (status /= SPARSE_STATUS_SUCCESS) then
               !    print *, "Error in solving the system"
               !    stop
               !end if
               !! Solve the linear system (A - sigma*I) * y = x using PARDISO
               !!call zcopy ( n, workd(ipntr(1)),1, workd(ipntr(2)), 1)
               !!phase = 33  ! Back substitution and iterative refinement
               !!call pardiso(pt, maxfct, mnum, mtype, phase, nnn, a, ia, ja, idum, nrhs, iparm, msglvl, workd(ipntr(1)), workd(ipntr(2)), error)
               !!if (error /= 0) then
               !!    print *, 'PARDISO error during solve phase:', error
               !!    stop
               !!end if
           else
               print *, 'Error: ido has unexpected value ', ido
               stop
           end if
           !print*, 'Solve iteration completed ... '
       
           ! ARPACK iteration
           call znaupd(ido, bmat, nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)
           if (info /= 0) then
               print *, 'Error with znaupd during iteration, info = ', info
               stop
           end if
           !print*, 'znaupd iteration finished ...'
           !iter = iter + 1
!          ! resid_norm = sqrt(sum(abs(workd(ipntr(2):ipntr(2) + nn - 1))**2))
           !resid_norm = sqrt(sum(abs(resid)**2))
           !print *, 'Iteration:', iter, 'Residual norm:', resid_norm, 'IDO:', ido
       end do
       print*, 'Solve completed ... '
    else
       do while (ido /= 99)
           if (ido == -1 .or. ido == 1) then
               ! Print the input vector for debugging
               !print *, 'Input vector:', workd(ipntr(1):ipntr(1)+nn-1)
               
               ! Perform sparse matrix-vector multiplication
               call sparse_matvec(nn, row_ptr, col_ind, values, &
                                  workd(ipntr(1):ipntr(1)+nn-1), &
                                  workd(ipntr(2):ipntr(2)+nn-1))
       
               ! Print the output vector for debugging
               !print *, 'Output vector:', workd(ipntr(2):ipntr(2)+nn-1)
           else
               print *, 'Error: ido has unexpected value ', ido
               stop
           end if
       
           ! ARPACK iteration
           call znaupd(ido, bmat, nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)
           if (info /= 0) then
               print *, 'Error with znaupd during iteration, info = ', info
               stop
           end if
       end do
    end if

    !do while (ido /= 99)
    !   if (ido .eq. -1 .or. ido .eq. 1 ) then
    !       call ccopy( n, workd(ipntr(1)),1, workd(ipntr(2)), 1)
    !       !call cgttrs('N', n, 1, dl, dd, du, du2, ipiv, workd(ipntr(2)), n, ierr)
    !       call sparse_matvec(nn, row_ptr, col_ind, values, workd(ipntr(1):ipntr(1)+nn-1), workd(ipntr(2):ipntr(2)+nn-1))
    !       if ( ierr .ne. 0 ) then
    !           print*, ' '
    !           print*, ' ERROR with _gttrs in _NDRV2.'
    !           print*, ' '
    !       end if
    !  else
    !      print *, 'Error: ipntr(1) or ipntr(2) out of bounds'
    !      stop
    !  end if
    !  call znaupd(ido, bmat, nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)
    !  if (info /= 0) then
    !      print *, 'Error with znaupd during iteration, info = ', info
    !      print *, 'iparam: ', iparam
    !      print *, 'ipntr: ', ipntr
    !      stop
    !  end if
    !end do

    !print *, 'before zneupd', ido
    !d = 0.0d0

    ! Debug print for extreme values
    if (maxval(abs(resid)) > 1e10) then
        print *, 'Warning: resid contains extremely large values.'
    end if

    !call zneupd(.false., 'A', select, d, v, nn, sigma, resid, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, ierr)
    if (saveRitz) then
       allocate(EVectors(nn, nev))
       call zneupd(.true.,'A', select, d, z, nn, sigma, workev, bmat,nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr,workd,workl, lworkl, rwork, info )
       nconv = iparam(5)
       do j=1, nconv
          call av(n, v(1,j), ax)
          call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
          rd(j,1) = dble (d(j))
          rd(j,2) = dimag (d(j))
          rd(j,3) = dznrm2 (n, ax, 1)
          rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
       end do
       call dmout (6, nconv, 3, rd, ncv, -6, 'Ritz values (Real, Imag) and relative residuals')
       ELoc(:nev) = real(d(:nev))
       EVectors = reshape(z, (/nn, nev/))
       ! Open file for writing
       open(unit=10, file='eigenvectors.txt', status='replace')
       
       ! Write eigenvectors to file
       ! Write eigenvectors to file
       print*, "ucell:", cell
       do i = 1, nev
           do j = 1, nn
               write(10, '(ES20.12,ES20.12)', advance='no') real(EVectors(j, i)), aimag(EVectors(j, i))
           end do
           write(10, *)
       end do
       close(10)
    else
       call zneupd(.false.,'A', select, d, z, nn, sigma, workev, bmat,nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr,workd,workl, lworkl, rwork, info )

       ELoc(:nev) = real(d(:nev))
    end if
    ! Debug prints after zneupd
    if (any(d /= d)) then
        print *, 'Error: d contains NaN values after zneupd.'
    endif
    if (info /= 0) then
        print *, 'Error with zneupd, ierr = ', info
        stop
    end if
    print *, "eigenvalues: ", d


    !if (useShift) then
    !    ! Termination and release memory
    !    phase = -1  ! Release internal memory
    !    call pardiso(pt, maxfct, mnum, mtype, phase, nnn, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    !end if
      
     deallocate(resid, v, workd, workl, rwork, d, iparam, ipntr, select, row_ptr, col_ind, values, rand_real, rand_imag)
      
        !if (N /= 3) then
        !    call test_3x3_matrix()
        !end if


end subroutine DiagHamSparse

!subroutine DiagH0TAPW(N, ns, is, ELoc, eigvec, KLoc, cell_real, H0, maxN, hopp, NList, Nneigh, neighCell,neig)
subroutine DiagH0TAPW(N, ns, is, ELoc, KLoc, cell_real, H0, maxN, hopp, NList, Nneigh, neighCell,neig)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species, RAt
    use tbpar, only : U
    use cell,                 only : rcell, ucell
    use constants,            only : pi
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is, neig
    real(dp), intent(out) :: ELoc(N)
!    real(dp), intent(out) :: eigvec(N)
    real(dp), intent(in) :: KLoc(3), cell_real(3,3)
    real(dp), intent(inout) :: H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    real(dp) :: tol
    real(dp) :: aGtemp

    integer :: i, j, in, info
    real(dp) :: R(3), zz, resid_norm

    ! ARPACK parameters and variables
    integer :: nev, ncv, lworkl, ido, ierr
    character(1) :: bmat
    character(2) :: which
    complex(dp), allocatable :: resid(:), v(:,:), workd(:), d(:), workev(:)
    complex(dp), allocatable :: workl(:)
    complex(dp), allocatable :: z(:,:)
    double precision, allocatable :: rwork(:)
    integer, allocatable :: iparam(:), ipntr(:)
    logical, allocatable :: select(:)
    integer :: max_iter, nn, iter
    logical :: useShift
    double precision :: shift

    integer, allocatable :: row_ptr(:), col_ind(:)
    complex(dp), allocatable :: values(:)
    real(dp), allocatable :: rand_real(:), rand_imag(:)
    complex(dp) :: sigma

    ! PARDISO variables
    !integer, parameter :: mtype = 13  ! Complex unsymmetric matrix
    !integer :: pt(64), iparm(64), maxfct, mnum, phase, error, msglvl
    !integer :: nnn, nrhs
    !complex(dp), allocatable :: a(:)
    !integer, allocatable :: ia(:), ja(:)
    !complex(dp) :: ddum
    !integer :: idum

    INTEGER :: nnn, nnz
    !complex(dp) :: x(N) ! Pardiso needs X(N) even if sol is returned in b

    !complex(dp) :: b(N)
!C.. Internal solver memory pointer for 64-bit architectures
    INTEGER*8 pt(64)
!C.. Internal solver memory pointer for 32-bit architectures
!C.. INTEGER*4 pt(64)

!C.. This is OK in both cases
!    TYPE(MKL_PARDISO_HANDLE) pt(64)
!C.. All other variables
    INTEGER maxfct, mnum, mtype, phase, nrhs, msglvl
    INTEGER iparm(64)
    !allocate(row_ptr(N+1), col_ind(nnz_temp), values(nnz_temp))
    INTEGER, allocatable :: ia(:) ! row_ptr
    INTEGER, allocatable :: ja(:) ! col_ind
    INTEGER idum(1)
    complex(dp) ::  ddum(1)
    integer error
    complex(dp), allocatable :: a(:) ! values
    
    DATA nrhs /1/, maxfct /1/, mnum /1/

    complex(dp), allocatable :: EVectors(:,:)
    logical :: saveRitz, symmetric
    integer :: nconv
    complex(dp), allocatable :: ax(:)
    double precision, allocatable :: rd(:,:)
    double precision :: dznrm2, dlapy2

    ! === TAPW variables ===
    integer :: NG, Nlabel, M
    integer, allocatable :: label(:)
    real(dp), allocatable :: Gx(:), Gy(:)
    complex(dp), allocatable :: XArray(:,:), Hproj(:,:)
    real(dp), allocatable :: eigvals(:)
    complex(dp), allocatable :: ZWorkLoc(:)
    real(dp), allocatable :: DWorkLoc(:)
    integer :: lwork

    real(dp) :: k_ref(2)
    integer :: NGrange

    
    ! Parameters
    !integer, parameter :: dp = c_double
    !integer, parameter :: int_kind = c_int
    
    ! Variables
    !type(C_PTR) :: mkl_handle
    !!integer(int_kind) :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !!integer(int_kind), allocatable :: row_ptr(:), col_ind(:)
    !integer :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !integer, allocatable :: row_ptr(:), col_ind(:)
    !real(dp), allocatable :: values(:), b(:), x(:), temp_value
    !type(SPARSE_MATRIX_DESCR) :: descr


    !integer :: nev, ncv, lworkl 
 
    nev=neig 
    ncv=nev*10
    lworkl= 3*NCV**2 + 5*NCV

    max_iter = 10000
    allocate(resid(N), v(N, ncv), workd(3*N), workl(lworkl), rwork(ncv), d(nev+1), iparam(11), ipntr(14), select(ncv), rand_real(N), rand_imag(N))
    allocate(z(N,nev))
    allocate(workev(2*ncv))
    allocate(ax(N))
    allocate(rd(ncv,3))

    ! Ensure correct size of nev and ncv
    if ( (nev < 1) .or. (nev >= ncv) .or. (ncv > N) ) then
        print *, 'Error: invalid parameters - nev=', nev, 'ncv=', ncv, 'N=', N
        stop
    end if

    ! Initialize arrays
    !v = 0.0d0
    !workd = 0.0d0
    !workl = 0.0d0
    !rwork = 0.0d0
    !d = 0.0d0

    bmat = 'I'

    iparam(1) = 1
    iparam(3) = max_iter
    call MIO_InputParameter('Diag.SparseSetShift',shift,0.0_dp)
    call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)
    call MIO_InputParameter('Diag.SparseSaveRitz',saveRitz,.false.)
    if (useShift) then
       which = 'LM'
       iparam(7) = 3
       sigma = cmplx(shift, 0.0_dp)
       print*, "We will find eigenvalues close to ", shift
    else
       which = 'SM'
       sigma = (0.0_dp, 0.0_dp)
       iparam(7) = 1
    end if
    ! Initialize ARPACK parameters

    call MIO_InputParameter('Diag.SparseSetTol',tol,0.1_dp)

    ido = 0
    info = 0

    !iparam(4) = 1
    !iparam(7) = 1

    nn = N

    ! Initialize the starting vector resid with random values
    !call random_number(rand_real)
    !call random_number(rand_imag)
    !resid = cmplx(rand_real, rand_imag)
    !resid = 0

    ! Check resid for initial state
    if (size(resid) /= N) then
        print *, 'Error: resid size mismatch: ', size(resid), ' expected: ', N
        stop
    end if

     ! Debug print for resid
    if (any(resid /= resid)) then
        print *, 'Error: resid contains NaN values initially.'
        stop
    end if

    resid_norm = sqrt(sum(abs(resid)**2))

    ! Debug prints
    !print *, 'Initial resid norm: ', resid_norm
    !print *, 'Initial iparam: ', iparam
    !print *, 'nn, nev, ncv: ', nn, nev, ncv

    ! Create CSR sparse matrix storage
    print *, "initialize the sparse matrix and put it in csr format"
    call initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell_real, row_ptr, col_ind, values,sigma)
    print *, "done"
    !call test_sparse_matvec(nn, row_ptr, col_ind, values)
    symmetric = is_structurally_symmetric(values, row_ptr, col_ind, n)

    print*, "is it symmetric?", symmetric

    ! Debug prints for CSR matrix
    if (any(values /= values)) then
        print *, 'Error: values contains NaN values after initialization.'
        stop
    end if

    
    if (maxval(abs(values)) > 1e10) then
        print *, 'Warning: values contains extremely large values.'
    end if

    !if (any(col_ind < 1 .or. col_ind > N)) then
    !    print *, 'Error: col_ind contains invalid indices after initialization.'
    !    stop
    !end if

    ! === TAPW settings (set these appropriately) ===
    NG = 20         ! for example
    Nlabel = 4      ! typically 2 (e.g., sublattice A/B)
    

    aGtemp = sqrt(cell_real(1,1)**2 + cell_real(2,1)**2)
    k_ref = (/ 4.0_dp*pi/(3.0_dp*aGtemp), 0.0_dp /)   ! set appropriately

    NGrange = 4   ! e.g., use all G such that -4 ≤ n1,n2 ≤ +4
    
    !call generate_G_list_from_rcell(rcell, NGrange, Gx, Gy, NG)
    call generate_shifted_G_list(rcell, k_ref, NGrange, Gx, Gy, NG)

    M = NG * Nlabel

    
    allocate(label(N))
    do i = 1, N
       label(i) = Species(i)
    end do
    allocate(XArray(N, M))
    call build_X(XArray, Rat(1,:), Rat(2,:), label, Gx, Gy, N, NG, Nlabel)


    allocate(Hproj(M, M))
    call transform_sparse_hamiltonian(N, M, row_ptr, col_ind, values, XArray, Hproj)
    if (.not. allocated(values)) stop "values not allocated after transform"
    print *, "size(values) =", size(values)

    !call zheev('V', 'U', M, Hproj, M, eigvals, work, lwork, rwork, info)
    !call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
    allocate(eigvals(M))
    lwork = 2*M
    allocate(ZWorkLoc(2*M), DWorkLoc(3*M))
    
    call ZHEEV('V','U', M, Hproj, M, eigvals, ZWorkLoc, lwork, DWorkLoc, info)
    
    if (info /= 0) then
       print *, "Diagonalization failed: ZHEEV info =", info
       stop
    end if

    if (.not. allocated(eigvals)) stop "eigvals not allocated"
    !if (.not. allocated(ELoc)) then
    !    print *, "ELoc not allocated"
    !    stop
    !endif
    if (size(ELoc) < M) then
        print *, "ERROR: size(ELoc) = ", size(ELoc), " but M = ", M
        stop
    endif

    !ELoc = eigvals
    !ELoc(1:nev) = eigvals
    !eigvec = ZWorkLoc(:, 1:M_band)
    ELoc = 0.0_dp
    ELoc(1:M) = eigvals

! eigvals now contains eigenvalues of projected H


    print *, "Deallocating arrays..."
    
    if (allocated(eigvals)) then
        deallocate(eigvals)
        print *, "Deallocated eigvals"
    endif
    
    if (allocated(ZWorkLoc)) then
        deallocate(ZWorkLoc)
        print *, "Deallocated ZWorkLoc"
    endif
    
    if (allocated(DWorkLoc)) then
        deallocate(DWorkLoc)
        print *, "Deallocated DWorkLoc"
    endif
    
    if (allocated(XArray)) then
        deallocate(XArray)
        print *, "Deallocated XArray"
    endif
    
    if (allocated(Hproj)) then
        deallocate(Hproj)
        print *, "Deallocated Hproj"
    endif
    
    if (allocated(Gx)) then
        deallocate(Gx)
        print *, "Deallocated Gx"
    endif
    
    if (allocated(Gy)) then
        deallocate(Gy)
        print *, "Deallocated Gy"
    endif
    
    if (allocated(label)) then
        deallocate(label)
        print *, "Deallocated label"
    endif
    
    if (allocated(row_ptr)) then
        deallocate(row_ptr)
        print *, "Deallocated row_ptr"
    endif
    
    if (allocated(col_ind)) then
        deallocate(col_ind)
        print *, "Deallocated col_ind"
    endif
    
    if (allocated(values)) then
        deallocate(values)
        print *, "Deallocated values"
    endif
    
    if (allocated(rand_real)) then
        deallocate(rand_real)
        print *, "Deallocated rand_real"
    endif
    
    if (allocated(rand_imag)) then
        deallocate(rand_imag)
        print *, "Deallocated rand_imag"
    endif
    
    if (allocated(a)) then
        deallocate(a)
        print *, "Deallocated a"
    endif
    
    if (allocated(EVectors)) then
        deallocate(EVectors)
        print *, "Deallocated EVectors"
    endif
    
    if (allocated(ax)) then
        deallocate(ax)
        print *, "Deallocated ax"
    endif
    
    if (allocated(rd)) then
        deallocate(rd)
        print *, "Deallocated rd"
    endif
    
    if (allocated(XArray)) then
        deallocate(XArray)
        print *, "Deallocated XArray"
    endif
    
    if (allocated(Hproj)) then
        deallocate(Hproj)
        print *, "Deallocated Hproj"
    endif
    
    if (allocated(ZWorkLoc)) then
        deallocate(ZWorkLoc)
        print *, "Deallocated ZWorkLoc"
    endif
    
    if (allocated(DWorkLoc)) then
        deallocate(DWorkLoc)
        print *, "Deallocated DWorkLoc"
    endif
    
    print *, "Done deallocating."

end subroutine DiagH0TAPW

subroutine generate_G_list_from_rcell(rcell, NGrange, Gx, Gy, NG)
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
  integer, intent(in) :: NGrange
  real(wp), intent(in) :: rcell(3,3)
  real(wp), allocatable, intent(out) :: Gx(:), Gy(:)
  integer, intent(out) :: NG

  integer :: i, j, count, nmax
  real(wp) :: b1(2), b2(2)

  b1 = rcell(1:2,1)
  b2 = rcell(1:2,2)

  nmax = (2*NGrange + 1)**2
  allocate(Gx(nmax), Gy(nmax))
  count = 0

  do i = -NGrange, NGrange
     do j = -NGrange, NGrange
        count = count + 1
        Gx(count) = i*b1(1) + j*b2(1)
        Gy(count) = i*b1(2) + j*b2(2)
     end do
  end do

  NG = count

  ! Optionally trim arrays to actual size
  if (count < nmax) then
     call move_alloc(Gx, Gx, stat=i)
     call move_alloc(Gy, Gy, stat=i)
     allocate(Gx(count), Gy(count))
     Gx(:) = Gx(1:count)
     Gy(:) = Gy(1:count)
  end if
end subroutine generate_G_list_from_rcell

subroutine generate_shifted_G_list(rcell, k_ref, NGrange, Gx, Gy, NG)
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
  real(wp), intent(in) :: rcell(3,3)      ! Reciprocal lattice
  real(wp), intent(in) :: k_ref(2)        ! Reference k-point (kx, ky)
  integer, intent(in) :: NGrange          ! G-grid range (±NGrange)
  real(wp), allocatable, intent(out) :: Gx(:), Gy(:)
  integer, intent(out) :: NG              ! Total number of G-vectors

  real(wp) :: b1(2), b2(2)
  integer :: i, j, count, nmax

  b1 = rcell(1:2,1)
  b2 = rcell(1:2,2)

  nmax = (2*NGrange + 1)**2
  allocate(Gx(nmax), Gy(nmax))

  count = 0
  do i = -NGrange, NGrange
     do j = -NGrange, NGrange
        count = count + 1
        Gx(count) = k_ref(1) + i*b1(1) + j*b2(1)
        Gy(count) = k_ref(2) + i*b1(2) + j*b2(2)
     end do
  end do

  NG = count
end subroutine generate_shifted_G_list

subroutine transform_sparse_hamiltonian(N, M, row_ptr, col_ind, values, X, Hproj)
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none

  integer, intent(in) :: N, M
  integer, intent(in) :: row_ptr(N+1), col_ind(:)
  complex(wp), intent(in) :: values(:)
  complex(wp), intent(in) :: X(N,M)
  complex(wp), intent(out) :: Hproj(M,M)

  complex(wp), allocatable :: Y(:,:)
  integer :: i, k, j, nnz
  complex(wp) :: tmp

  ! Validate CSR size
  nnz = row_ptr(N+1) - 1
  if (size(values) < nnz) then
     print *, "ERROR: size(values) <", nnz
     stop "CSR format inconsistency: not enough values"
  endif
  if (maxval(col_ind(1:nnz)) > N) then
     print *, "ERROR: max(col_ind) >", N
     stop "col_ind contains out-of-bounds indices for X"
  endif

  ! Allocate and zero Y = H * X
  allocate(Y(N,M))
  Y = (0.0_wp, 0.0_wp)

  ! Perform sparse matrix multiplication Y = H * X
  do i = 1, N
     do k = row_ptr(i), row_ptr(i+1) - 1
        j = col_ind(k)
        if (j < 1 .or. j > N) then
           print *, "ERROR: j = col_ind(k) = ", j, " out of bounds at i=", i, " k=", k
           stop
        endif
        if (k < 1 .or. k > size(values)) then
           print *, "ERROR: k=", k, " out of bounds (values size=", size(values), ")"
           stop
        endif
        Y(i,:) = Y(i,:) + values(k) * X(j,:)
     end do
  end do

  ! Compute Hproj = X^† * Y
  call zgemm('C', 'N', M, M, N, cmplx(1.0_wp, 0.0_wp), X, N, Y, N, cmplx(0.0_wp, 0.0_wp), Hproj, M)

  deallocate(Y)
end subroutine transform_sparse_hamiltonian

subroutine build_X(X, xcoord, ycoord, label, Gx, Gy, N_orbit, N_G, N_label)
  implicit none
  integer, intent(in) :: N_orbit, N_G, N_label
  real*8, intent(in) :: xcoord(N_orbit), ycoord(N_orbit)          ! positions
  integer, intent(in) :: label(N_orbit)                 ! orbital label per atom (1..N_label)
  real*8, intent(in) :: Gx(N_G), Gy(N_G)                ! G vectors
  complex*16, intent(out) :: X(N_orbit, N_G * N_label)

  integer :: i_orb, i_G, i_label, col
  integer :: count(N_label)
  real*8 :: norm(N_label)
  complex*16 :: phase

  count = 0

  ! First count how many orbitals of each label
  do i_orb = 1, N_orbit
     count(label(i_orb)) = count(label(i_orb)) + 1
  end do

  ! Compute 1/sqrt(N_label) normalization
  do i_label = 1, N_label
     if (count(i_label) > 0) then
        norm(i_label) = 1.0d0 / sqrt(dble(count(i_label)))
     else
        norm(i_label) = 0.0d0
     end if
  end do

  ! Build X matrix
  do i_G = 1, N_G
     do i_label = 1, N_label
        col = (i_G - 1) * N_label + i_label
        do i_orb = 1, N_orbit
           if (label(i_orb) == i_label) then
              phase = cmplx(0.0d0, 1.0d0) * (Gx(i_G) * xcoord(i_orb) + Gy(i_G) * ycoord(i_orb))
              X(i_orb, col) = norm(i_label) * exp(phase)
           else
              X(i_orb, col) = (0.0d0, 0.0d0)
           end if
        end do
     end do
  end do
end subroutine build_X

logical function is_structurally_symmetric(a, ia, ja, n)
    implicit none
    integer, intent(in) :: ia(:), ja(:)
    complex(dp), intent(in) :: a(:)
    integer, intent(in) :: n
    integer :: i, j, k, m
    logical :: found

    ! Declare arrays to store transpose of the CSR matrix
    integer, allocatable :: trans_ia(:), trans_ja(:)
    integer :: nnz

    ! Initialize the return value
    is_structurally_symmetric = .true.

    ! Number of non-zeros in the matrix
    nnz = size(ja)

    ! Allocate transpose arrays
    allocate(trans_ia(n + 1))
    allocate(trans_ja(nnz))

    ! Initialize trans_ia array
    trans_ia = 0

    ! Count the number of entries in each column of the original matrix
    do i = 1, n
        do k = ia(i), ia(i+1) - 1
            j = ja(k)
            trans_ia(j+1) = trans_ia(j+1) + 1
        end do
    end do

    ! Convert counts to starting indices
    trans_ia(1) = 1
    do i = 2, n + 1
        trans_ia(i) = trans_ia(i) + trans_ia(i-1)
    end do

    ! Fill the trans_ja array
    do i = 1, n
        do k = ia(i), ia(i+1) - 1
            j = ja(k)
            m = trans_ia(j)
            trans_ja(m) = i
            trans_ia(j) = m + 1
        end do
    end do

    ! Restore the trans_ia array to correct starting indices
    do i = n, 1, -1
        trans_ia(i+1) = trans_ia(i)
    end do
    trans_ia(1) = 1

    ! Compare the structure of the original matrix with its transpose
    do i = 1, n
        do k = ia(i), ia(i+1) - 1
            j = ja(k)
            found = .false.
            do m = trans_ia(j), trans_ia(j+1) - 1
                if (trans_ja(m) == i) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                is_structurally_symmetric = .false.
                return
            end if
        end do
    end do

    ! Clean up
    deallocate(trans_ia)
    deallocate(trans_ja)
end function is_structurally_symmetric


real(dp) function norm2(x)
    implicit none
    complex(dp), intent(in) :: x(:)
    norm2 = sqrt(sum(abs(x)**2))
end function norm2

subroutine test_3x3_matrix()
    use constants, only : cmplx_i
    implicit none
    real(dp) :: H0_test(3, 3)
    real(dp) :: KLoc(3) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    real(dp) :: cell(3, 3) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    complex(dp) :: hopp(10, 3) = 0.0_dp
    integer :: NList(10, 3) = 0
    integer :: Nneigh(3) = 0
    integer :: neighCell(3, 10, 3) = 0
    real(dp) :: ELoc_test(3)
    real(dp) :: expected_eigenvalues(3)
    integer :: N_test = 3, maxN_test = 10
    integer :: i

    H0_test = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                       0.0_dp, 2.0_dp, 0.0_dp, &
                       0.0_dp, 0.0_dp, 3.0_dp], [3, 3])

    expected_eigenvalues = (/1.0_dp, 2.0_dp, 3.0_dp/)

    !print *, 'Testing 3x3 matrix eigenvalues:'
    ! Call the sparse matrix solver for the 3x3 test matrix
    !call DiagHamSparse(N_test, 1, 1, ELoc_test, KLoc, cell, H0_test, maxN_test, hopp, NList, Nneigh, neighCell)

    ! Print results
    !print *, 'Computed eigenvalues: ', ELoc_test
    !print *, 'Expected eigenvalues: ', expected_eigenvalues

    ! Check results
    do i = 1, 3
        if (abs(ELoc_test(i) - expected_eigenvalues(i)) > 1.0e-6) then
            print *, 'Error: Eigenvalue ', i, ' does not match. Computed: ', ELoc_test(i), ', Expected: ', expected_eigenvalues(i)
        else
            print *, 'Eigenvalue ', i, ' matches.'
        end if
    end do
end subroutine test_3x3_matrix

subroutine initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values,sigma)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species
    use tbpar, only : U
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
    real(dp), intent(in) :: KLoc(3), cell(3,3)
    real(dp), intent(inout) :: H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)
    complex(dp), allocatable, intent(out) :: values(:)
    integer :: i, j, k, l, nnz, in
    real(dp) :: R(3)
    integer :: nnz_temp
    double precision :: alpha
    logical :: useShift, found
    complex(dp) :: sigma

    call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)
    !call MIO_InputParameter('Diag.SparseAlpha',alpha,0.0_dp)
    !if (useShift) then
    !   do i = 1, N
    !       H0(i) = H0(i) + alpha
    !   end do
    !end if
  

    ! First pass to count non-zero elements
    nnz_temp = 0
    do i = 1, N
        if (useShift) then
           nnz_temp = nnz_temp + 1
        else
           if (H0(i) /= 0.0_dp) nnz_temp = nnz_temp + 1
        end if
        do j = 1, Nneigh(i)
            if (hopp(j,i) /= 0.0_dp) nnz_temp = nnz_temp + 1
        end do
    end do
    if (edgeHopp) then
        do i = 1, nQ
            do j = 1, nEdgeN(i)
                if (edgeH(j,i) /= 0.0_dp) nnz_temp = nnz_temp + 1
            end do
        end do
    end if

    ! Allocate space for CSR arrays
    !allocate(row_ptr(N+1), col_ind(nnz_temp), values(nnz_temp))
    allocate(row_ptr(N+1), col_ind(nnz_temp + 1000), values(nnz_temp + 1000))

    row_ptr(1) = 1
    l = 1  ! Use a different variable to increment the non-zero element index


    ! Populate CSR arrays directly
    do i = 1, N
        row_ptr(i) = l
        ! Diagonal element
        if (H0(i) /= 0.0_dp) then
            !if (l > nnz_temp) then
            !   print *, "CSR OVERFLOW: l =", l, " > nnz_temp =", nnz_temp
            !   stop "Sparse matrix allocation overflow in initialize_sparse_matrix"
            !endif 
            values(l) = H0(i) - sigma
            col_ind(l) = i
            l = l + 1
        end if

        ! Off-diagonal elements
        !if (useShift) then
        !   do j = 1, Nneigh(i)
        !      in = NList(j, i)
        !      if (in > i) then
        !          R = matmul(cell, neighCell(:, j, i))
        !          if (hopp(j, i) /= 0.0_dp) then
        !             !found = .false.
        !             !! Check if we already have an entry for this off-diagonal element
        !             !do k = row_ptr(i), l - 1
        !             !    if (col_ind(k) == in) then
        !             !        print*, "entering this loop means you are working on a small system, no?"
        !             !        values(k) = values(k) - hopp(j, i) * exp(cmplx_i * dot_product(KLoc, R))
        !             !        found = .true.
        !             !        exit
        !             !    end if
        !             !end do
        !             !! If no entry exists, add a new one
        !             !if (.not. found) then
        !                 values(l) = -hopp(j, i) * exp(cmplx_i * dot_product(KLoc, R))
        !                 col_ind(l) = in
        !                 l = l + 1
        !             !end if
        !          end if
        !      end if
        !   end do
        !else
           do j = 1, Nneigh(i)
               in = NList(j,i)
               R = matmul(cell,neighCell(:,j,i))
               if (hopp(j,i) /= 0.0_dp) then
                  found = .false.
                  ! Check if we already have an entry for this off-diagonal element
                  do k = row_ptr(i), l - 1
                      if (col_ind(k) == in) then
                          values(k) = values(k) - hopp(j, i) * exp(cmplx_i * dot_product(KLoc, R))
                          found = .true.
                          exit
                      end if
                  end do
                  ! If no entry exists, add a new one
                  if (.not. found) then
                      if (l > nnz_temp) then
                         print *, "CSR OVERFLOW: l =", l, " > nnz_temp =", nnz_temp
                         stop "Sparse matrix allocation overflow in initialize_sparse_matrix"
                      endif 
                      values(l) = -hopp(j, i) * exp(cmplx_i * dot_product(KLoc, R))
                      col_ind(l) = in
                      l = l + 1
                  end if
               end if
           end do
        !end if

        ! Edge hopping elements
        !if (edgeHopp) then
        !    do j = 1, nQ
        !        do k = 1, nEdgeN(j)
        !            in = NeI(k,j)
        !            R = matmul(cell, NedgeCell(:,k,j))
        !            if (edgeH(k,j) /= 0.0_dp) then
        !                values(l) = edgeH(k,j) * exp(cmplx_i * dot_product(KLoc, R))
        !                col_ind(l) = in
        !                l = l + 1
        !            end if
        !        end do
        !    end do
        !end if
    end do


    nnz = l - 1
    row_ptr(N+1) = nnz + 1

    ! Trim the allocated arrays to actual size
    !if (useShift) then
    !   values = values(1:nnz)
    !   col_ind = col_ind(1:nnz)
    !end if

    ! Debug print for sparse matrix
    print *, 'Sparse matrix row_ptr: ', row_ptr(1:min(N+1,10))
    print *, 'Sparse matrix col_ind: ', col_ind(1:min(nnz_temp,10))
    print *, 'Sparse matrix values: ', values(1:min(nnz,10))
    print *, "CSR max l =", l-1, "allocated nnz_temp =", nnz_temp
end subroutine initialize_sparse_matrix


!subroutine initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values)
!    use constants, only : cmplx_i
!    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
!    use scf, only : charge, Zch
!    use atoms, only : Species
!    use tbpar, only : U
!    implicit none
!    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
!    real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
!    complex(dp) :: HLoc(N,N)
!    complex(dp), intent(in) :: hopp(maxN,N)
!    integer, allocatable :: row_ptr(:), col_ind(:)
!    complex(dp), allocatable, intent(out) :: values(:)
!    integer :: i, j, k, nnz, in
!    real(dp) :: R(3)
!
!    HLoc = 0.0_dp
!    do i=1,N
!       HLoc(i,i) = H0(i)
!       if (ns==2) then
!          !zz = charge(1,i)*charge(2,i) ! Zch
!          if (is==1) then
!             HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
!          else
!             HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
!          end if
!       end if
!       do j=1,Nneigh(i)
!          in = NList(j,i)
!          R = matmul(cell,neighCell(:,j,i))
!          !print*, "nownow", KLoc
!          !print*, "niwniw", R
!          HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
!       end do
!    end do
!    if (edgeHopp) then
!       do i=1,nQ
!          do j=1,nEdgeN(i)
!             in = NeI(j,i)
!             R = matmul(cell,NedgeCell(:,j,i))
!             HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(cmplx_i*dot_product(KLoc,R))
!          end do
!       end do
!    end if
!    ! Count the number of non-zero elements
!    nnz = count(HLoc /= 0.0_dp)
!
!    ! Allocate space for CSR arrays
!    allocate(row_ptr(N+1), col_ind(nnz), values(nnz))
!
!    row_ptr(1) = 1
!    k = 1
!
!    ! Populate CSR arrays
!    do i = 1, N
!        do j = 1, N
!            if (HLoc(i, j) /= 0.0_dp) then
!                values(k) = HLoc(j, i)
!                col_ind(k) = j
!                k = k + 1
!            end if
!        end do
!        row_ptr(i+1) = k
!    end do
!    ! Parallelize the outer loop with OpenMP
!    !!$OMP PARALLEL DO PRIVATE(i, j, k_local) SHARED(row_ptr, col_ind, values, HLoc) REDUCTION(+:k)
!    !do i = 1, N
!    !    k_local = k  ! Local copy of k for each thread
!    !    do j = 1, N
!    !        if (HLoc(i, j) /= (0.0, 0.0)) then
!    !            values(k_local) = HLoc(j, i)
!    !            col_ind(k_local) = j
!    !            k_local = k_local + 1
!    !        end if
!    !    end do
!    !    row_ptr(i+1) = k_local
!    !    k = k_local  ! Update global k with the local k value
!    !end do
!    !!$OMP END PARALLEL DO
!
!    ! Debug print for sparse matrix
!    !print *, 'Sparse matrix row_ptr: ', row_ptr(1:min(N+1,10))
!    !print *, 'Sparse matrix col_ind: ', col_ind(1:min(nnz,10))
!    !print *, 'Sparse matrix values: ', values(1:min(nnz,10))
!end subroutine initialize_sparse_matrix

!subroutine sparse_matvec(nn, row_ptr, col_ind, values, x, y)
!    implicit none
!    integer, intent(in) :: nn, row_ptr(:), col_ind(:)
!    complex(dp), intent(in) :: values(:), x(:)
!    complex(dp), intent(out) :: y(nn)
!    integer :: i, j
!
!    y = 0.0_dp
!    do i = 1, nn
!        do j = row_ptr(i), row_ptr(i+1) - 1
!            y(i) = y(i) + values(j) * x(col_ind(j))
!        end do
!    end do
!end subroutine sparse_matvec

subroutine sparse_matvec(nn, row_ptr, col_ind, values, x, y)
    implicit none
    integer, intent(in) :: nn, row_ptr(:), col_ind(:)
    complex(dp), intent(in) :: values(:), x(:)
    complex(dp), intent(out) :: y(nn)
    integer :: i, j

    y = 0.0_dp

    ! Parallelize the outer loop with OpenMP
    !$OMP PARALLEL DO PRIVATE(i, j) SHARED(row_ptr, col_ind, values, x, y)
    do i = 1, nn
        do j = row_ptr(i), row_ptr(i+1) - 1
            y(i) = y(i) + values(j) * x(col_ind(j))
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine sparse_matvec

subroutine test_sparse_matvec(nn, row_ptr, col_ind, values)
    implicit none
    integer, intent(in) :: nn, row_ptr(:), col_ind(:)
    complex(dp), intent(in) :: values(:)
    complex(dp), allocatable :: x(:), y(:)
    real(dp), allocatable :: rand_real(:), rand_imag(:)

    ! Allocate and initialize test vectors
    allocate(x(nn), y(nn), rand_real(nn), rand_imag(nn))
    call random_number(rand_real)
    call random_number(rand_imag)
    x = cmplx(rand_real, rand_imag)

    ! Perform sparse matrix-vector multiplication
    call sparse_matvec(nn, row_ptr, col_ind, values, x, y)

    ! Print results for verification
    !print *, 'Test sparse matvec result: ', y(1:min(10, nn))
    !print *, 'Input vector x: ', x(1:min(10, nn))
    !print *, 'Row_ptr: ', row_ptr(1:min(size(row_ptr), 10))
    !print *, 'Col_ind: ', col_ind(1:min(size(col_ind), 10))
    !print *, 'Values: ', values(1:min(size(values), 10))
end subroutine test_sparse_matvec



subroutine DiagHamSparse2(N, ns, is, ELoc, KLoc, cell, H0, maxN, hopp, NList, Nneigh, neighCell)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species
    use tbpar, only : U
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
    real(dp), intent(out) :: ELoc(N)
    real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)

    integer :: i, j, k, iter, max_iter, nev, nnz
    real(dp) :: tol, alpha, beta
    complex(dp) :: mu
    real(dp), allocatable :: alpha_vec(:), beta_vec(:), diag(:), offdiag(:)
    complex(dp), allocatable :: q(:,:), v(:), w(:), r(:)
    integer, allocatable :: row_ptr(:), col_ind(:)
    complex(dp), allocatable :: values(:)
    real(dp), allocatable :: rand_real(:), rand_imag(:)

    max_iter = 1000
    tol = 1.0e-10
    nev = 100  ! Number of eigenvalues to compute

    ! Initialize sparse matrix storage
    call initialize_sparse_matrix2(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values)

    ! Allocate memory for Lanczos method
    allocate(q(N, nev+1), alpha_vec(nev), beta_vec(nev), v(N), w(N), r(N), rand_real(N), rand_imag(N))
    q = cmplx(0.0_dp, 0.0_dp)
    alpha_vec = 0.0_dp
    beta_vec = 0.0_dp

    ! Initialize the starting vector with random values
    call random_number(rand_real)
    call random_number(rand_imag)
    v = cmplx(rand_real, rand_imag)
    v = v / sqrt(sum(abs(v)**2))

    ! Shift-and-invert strategy parameters
    mu = cmplx(0.0_dp, 0.0_dp)  ! Fermi energy assumed to be zero

    ! Lanczos algorithm
    beta = 0.0_dp
    do iter = 1, nev
        if (iter > 1) then
            r = r - beta * q(:, iter-1)
        else
            r = v
        end if

        call shift_invert_sparse_matvec2(N, row_ptr, col_ind, values, r, w, mu)
        alpha = real(dot_product(r, w))
        r = w - alpha * r - beta * q(:, iter)
        beta = sqrt(sum(abs(r)**2))
        
        alpha_vec(iter) = alpha
        beta_vec(iter) = beta

        if (beta < tol) exit

        q(:, iter+1) = r / beta
    end do

    ! Tridiagonal matrix T construction
    diag = alpha_vec(1:iter)
    offdiag = beta_vec(2:iter)

    ! Compute eigenvalues of T using LAPACK routine
    call dstev('N', iter, diag, offdiag, diag, iter, diag, nnz)

    ! Copy the eigenvalues to output
    ELoc(1:nev) = diag(1:nev)

    ! Deallocate memory
    deallocate(q, alpha_vec, beta_vec, v, w, r, row_ptr, col_ind, values, rand_real, rand_imag)
end subroutine DiagHamSparse2

subroutine initialize_sparse_matrix2(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values)
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), ns, is
    real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    integer, intent(in) :: neighCell(3,maxN,N)
    integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)
    complex(dp), allocatable, intent(out) :: values(:)
    integer :: i, j, k, nnz

    nnz = 0
    allocate(row_ptr(N+1))
    row_ptr(1) = 1
    do i = 1, N
        do j = 1, Nneigh(i)
            nnz = nnz + 1
        end do
        row_ptr(i+1) = nnz + 1
    end do

    allocate(col_ind(nnz), values(nnz))

    nnz = 0
    do i = 1, N
        do j = 1, Nneigh(i)
            nnz = nnz + 1
            col_ind(nnz) = NList(j, i)
            values(nnz) = hopp(j, i)
        end do
    end do
end subroutine initialize_sparse_matrix2

subroutine shift_invert_sparse_matvec2(N, row_ptr, col_ind, values, x, y, mu)
    implicit none
    integer, intent(in) :: N, row_ptr(:), col_ind(:)
    complex(dp), intent(in) :: values(:), x(:), mu
    complex(dp), intent(out) :: y(N)
    complex(dp), allocatable :: temp(:)
    integer :: i, j

    allocate(temp(N))
    y = cmplx(0.0_dp, 0.0_dp)
    temp = cmplx(0.0_dp, 0.0_dp)

    ! Perform (A - mu*I) * x
    do i = 1, N
        temp(i) = -mu * x(i)
        do j = row_ptr(i), row_ptr(i+1)-1
            temp(i) = temp(i) + values(j) * x(col_ind(j))
        end do
    end do

    ! Solve (A - mu*I) * y = x
    ! Use a simple iterative solver like Conjugate Gradient or GMRES
    call iterative_solver2(N, row_ptr, col_ind, values, temp, y)

    deallocate(temp)
end subroutine shift_invert_sparse_matvec2

subroutine iterative_solver2(N, row_ptr, col_ind, values, b, x)
    implicit none
    integer, intent(in) :: N, row_ptr(:), col_ind(:)
    complex(dp), intent(in) :: values(:), b(:)
    complex(dp), intent(out) :: x(N)
    integer :: max_iter, iter
    real(dp) :: tol, alpha, beta, rsold, rsnew
    complex(dp), allocatable :: r(:), p(:), Ap(:)

    max_iter = 1000
    tol = 1.0e-10
    allocate(r(N), p(N), Ap(N))

    ! Initial guess x = 0
    x = cmplx(0.0_dp, 0.0_dp)

    ! r = b - A*x
    call sparse_matvec2(N, row_ptr, col_ind, values, x, Ap)
    r = b - Ap
    p = r
    rsold = real(dot_product(conjg(r), r))

    do iter = 1, max_iter
        call sparse_matvec2(N, row_ptr, col_ind, values, p, Ap)
        alpha = rsold / real(dot_product(conjg(p), Ap))
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = real(dot_product(conjg(r), r))
        if (sqrt(rsnew) < tol) exit
        beta = rsnew / rsold
        p = r + beta * p
        rsold = rsnew
    end do

    deallocate(r, p, Ap)
end subroutine iterative_solver2

subroutine sparse_matvec2(N, row_ptr, col_ind, values, x, y)
    implicit none
    integer, intent(in) :: N, row_ptr(:), col_ind(:)
    complex(dp), intent(in) :: values(:), x(:)
    complex(dp), intent(out) :: y(N)
    integer :: i, j

    y = cmplx(0.0_dp, 0.0_dp)
    do i = 1, N
        do j = row_ptr(i), row_ptr(i+1)-1
            y(i) = y(i) + values(j) * x(col_ind(j))
        end do
    end do
end subroutine sparse_matvec2






subroutine DiagHamWF(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   !complex(dp), intent(out) :: H(N,N)
   !real(dp), intent(out) :: E(N)
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   !print*, cell
   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   !call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   call ZHEEV('V','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      !print*, "info =", info
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if
   

end subroutine DiagHamWF

subroutine DiagHamChern(N,ns,is,HLoc,ChernLoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   !complex(dp), intent(out) :: H(N,N)
   !real(dp), intent(out) :: E(N)
   complex(dp), intent(out) :: HLoc(N,N)
   complex(dp) :: HLocM1dx(N,N), HLocM1dy(N,N), H_derivativeDX(N,N), H_derivativeDY(N,N)
   complex(dp) :: dHdx(N,N)
   complex(dp) :: dHdy(N,N)
   complex(dp) :: eigvec(N,N), Vxmn(N,N), Vymn(N,N), temp_matrix(N,N), ones_matrix(N,N), difference_matrix(N,N), eigval_matrix(N,N), eigval_repeated_matrix(N, N) 
   real(dp), intent(out) :: ChernLoc(N)
   complex(dp) :: ChernLocSum(N)
   real(dp) :: ELoc(N)
   real(dp) :: eigval(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   real(dp) :: KLocM1(3)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info, k
   real(dp) :: R(3), zz, eps

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   real(8) :: delta_kx, delta_ky

   delta_kx = 0.01d0
   delta_ky = 0.01d0


   call MIO_InputParameter('Epsilon',eps,0.001_dp)

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         KLocM1(1) = KLoc(1)-delta_kx
         KLocM1(2) = KLoc(2)
         KLocM1(3) = KLoc(3)
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLocM1,R))
      end do
   end do
   HLocM1dy = HLoc

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         KLocM1(1) = KLoc(1)
         KLocM1(2) = KLoc(2)-delta_ky
         KLocM1(3) = KLoc(3)
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLocM1,R))
      end do
   end do
   HLocM1dy = HLoc

   !print*, cell
   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
      end do
   end do

   ! Assuming HLoc and HLocM1 are matrices, you would first initialize a matrix for the derivative
   H_derivativeDX = 0.0_dp
   H_derivativeDY = 0.0_dp
   
   ! Calculation of Hamiltonians remains the same, so I won't repeat that here.
   
   ! Now, compute the finite difference derivative
   do i=1,N
      do j=1,N
         H_derivativeDX(i,j) = (HLoc(i,j) - HLocM1dx(i,j)) / (2.0_dp * delta_kx)
         dHdx(i,j) = cmplx(0.0, 1.0) * HLoc(i,j) * H_derivativeDX(i,j)
         H_derivativeDY(i,j) = (HLoc(i,j) - HLocM1dy(i,j)) / (2.0_dp * delta_ky)
         dHdy(i,j) = cmplx(0.0, 1.0) * HLoc(i,j) * H_derivativeDY(i,j)
      end do
   end do

   !call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   call ZHEEV('V','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      !print*, "info =", info
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if

   ! Calculating Vxmn and Vymn
   eigvec = HLoc
   eigval = ELoc
   print*, "eigval", eigval

   Vxmn = matmul(matmul(transpose(eigvec), dHdx), eigvec)
   Vymn = matmul(matmul(transpose(eigvec), dHdy), eigvec)
   
   ! Setting diagonals to zero
   do j = 1, N
       Vxmn(j,j) = 0.0_dp
       Vymn(j,j) = 0.0_dp
   end do
   print*, "Vxmn", Vxmn
   
   !forall(i=1:N, j=1:N) ones_matrix(i, j) = 1.0d0
   !! Calculating Chern
   !do j = 1, N
   !    temp_matrix = Vxmn * transpose(Vymn) / (matmul(ones_matrix, transpose(eigval)) - eigval)**2 + eps
   !    do k = 1, N
   !        Chern(i) = Chern(i) + 2.0 * aimag(temp_matrix(j, k))
   !    end do
   !end do
   forall(i=1:N, j=1:N) ones_matrix(i, j) = 1.0_dp

   ! Calculating Chern
   forall(i=1:N, j=1:N) eigval_matrix(i, j) = eigval(j)
   forall(i=1:N, j=1:N) eigval_repeated_matrix(i, j) = eigval(i)
   difference_matrix = eigval_matrix - eigval_repeated_matrix
   !print*, "diffma", difference_matrix
   temp_matrix = Vxmn * transpose(Vymn) / (difference_matrix**2 + eps)
   !print*, "temp_ma", temp_matrix
   
   !ChernLocSum = 0.0_dp
   do j = 1, N
       do k = 1, N
           !ChernLocSum = ChernLocSum + temp_matrix(j, k)
           ChernLoc = ChernLoc + 2.0 * aimag(temp_matrix(j, k))
       end do
   end do
   !print*, ChernLoc
   !ChernLoc = 2.0_dp*aimag(ChernLocSum)
   ! Chern[i] = 2*np.imag(np.sum(Vxmn*(Vymn.T/((eigval*np.ones((N_orbit,N_orbit))).T-eigval)**2+eps), axis=1))
   

end subroutine DiagHamChern

subroutine DiagHamPDOS(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   !complex(dp), intent(out) :: H(N,N)
   !real(dp), intent(out) :: E(N)
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   !print*, cell
   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('V','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   !call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      !print*, "info =", info
      call MIO_Kill('Error in diagonalization','diag','DiagHamPDOS')
   end if

end subroutine DiagHamPDOS

subroutine DiagHamRashba(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : neighD
   !use cell,                  only : aG

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   !complex(dp), intent(out) :: H(N,N)
   !real(dp), intent(out) :: E(N)
   complex(dp), intent(out) :: HLoc(N*2,N*2) ! double the output Hamiltonian size to account for spin
   real(dp), intent(out) :: ELoc(N*2) ! same for eigenergies
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N*2-2)

   real(dp) :: Rashbahopp
   real(dp) :: lambdaR, acc, dist

   call MIO_InputParameter('lambdaRasbha',lambdaR,0.0_dp)
   acc = 2.46_dp/sqrt(3.0_dp)
   !print*, cell
   HLoc = 0.0_dp
   !print*, "diaghere3", N
   do i=1,N
      HLoc(i,i) = H0(i)
      HLoc(i+N,i+N) = H0(i) ! Repeat initial H0 a second time
      !print*, "diaghere3a"
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      !print*, "diaghere3b", i
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         HLoc(in+N,i+N) = HLoc(in+N,i+N) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
         if (layerIndex(i).eq.layerIndex(in) .and. Species(i).ne.Species(in) .and. dist<acc*1.1_dp) then
             Rashbahopp = lambdaR*(neighD(1,j,i)+cmplx_i*neighD(2,j,i))
             HLoc(in+N,i) = HLoc(in+N,i) - Rashbahopp*exp(cmplx_i*dot_product(KLoc,R))
             HLoc(i,in+N) = HLoc(i,in+N) - Rashbahopp*exp(-cmplx_i*dot_product(KLoc,R))
         end if
      end do
      !print*, "diaghere3c"
   end do
   !print*, "diaghere4"
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   !print*, "diaghere1"
   call ZHEEV('N','L',N*2,HLoc,N*2,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   !call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      !print*, "info =", info
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if
   !print*, "diaghere2"

end subroutine DiagHamRashba

subroutine DiagHamArpack(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   !complex(dp), intent(out) :: H(N,N)
   !real(dp), intent(out) :: E(N)
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   !call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   !call matvecA('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      !print*, "info =", info
      call MIO_Kill('Error in diagonalization','diag','DiagHamArpack')
   end if

end subroutine DiagHamArpack

subroutine DiagSpectralWeightNishi(N,ns,is,Pkc,E,K,KG,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: Pkc(N)
   complex(dp) :: Pkcaux
   real(dp), intent(out) :: E(N)
   real(dp), intent(in) :: K(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   integer :: kk
   

   call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   !aBN = aG*1.018181818
   !bncell(:,1) = [aBN,0.0_dp,0.0_dp]
   !bncell(:,2) = [aBN/2.0_dp,sqrt(3.0_dp)*aBN/2.0_dp,0.0_dp]
   !bncell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   cellts = gcell
   celltsp = gcell

   !print*, "cellts", cellts
   !print*, "celltsp", celltsp

   Hts = 0.0_dp
   Htsp = 0.0_dp
   Pkc = 0.0_dp

   !Rts = matmul(cellts,neighCell(:,j,i)) ! Moire
   !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
   RtsVec(1,:) = [matmul(cellts,[1,0,0])] 
   RtsVec(2,:) = [matmul(cellts,[0,1,0])] 
   RtsVec(3,:) = [matmul(cellts,[1,1,0])] 
   RtsVec(4,:) = [matmul(cellts,[-1,0,0])]
   RtsVec(5,:) = [matmul(cellts,[0,-1,0])]
   RtsVec(6,:) = [matmul(cellts,[-1,-1,0])]

   RtspVec(1,:) = [matmul(cellts,[1,0,0])] 
   RtspVec(2,:) = [matmul(cellts,[0,1,0])] 
   RtspVec(3,:) = [matmul(cellts,[1,1,0])] 
   RtspVec(4,:) = [matmul(cellts,[-1,0,0])]
   RtspVec(5,:) = [matmul(cellts,[0,-1,0])]
   RtspVec(6,:) = [matmul(cellts,[-1,-1,0])] 

   !Rtsx = RtsVec(1,1)
   !Rtsy = RtsVec(2,2)
  
   call MIO_InputParameter('Spectral.numberOfTS',nTS,1)

   
   !if (frac) call AtomsSetCart()
   !do i=1,N 
   !    print*, Rat(1,i), Rat(2,i)
   !    print*, floor((Rat(1,i)-0.1)/Rtsx), floor(Rat(2,i)/Rtsy)
   !    RtsVec(i,:) = matmul(cellts,[floor((Rat(1,i)-0.1)/Rtsx),floor(Rat(2,i)/Rtsy),0]) ! Moire
   !    RtspVec(i,:) = matmul(celltsp,[floor((Rat(1,i)-0.1)/Rtsx),floor(Rat(2,i)/Rtsy),0]) ! Moire
   !end do

   !kk = 0
   !do i=0,4 
   !    do j=0,4 
   !        kk = kk+1
   !        RtsVec(kk,:) = matmul(cellts,[i,j,0]) ! Moire
   !        RtspVec(kk,:) = matmul(celltsp,[i,j,0]) ! Moire
   !    end do
   !end do

   nTS = 6
   
   !print*, "hohoho", RtsVec
   !print*, "hehehe", RtsVec(1,:)
   do i1=1,nTS!SIZE(RtsVec)
      do i2=1,nTS!SIZE(RtspVec)
         Rts = RtsVec(i1,:)
         Rtsp = RtspVec(i2,:)
         do i=1,N
            Hts(i,i) = H0(i)
            Htsp(i,i) = H0(i)
            !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
            do j=1,Nneigh(i)
               in = NList(j,i)
               !Rts = matmul(cellts,[floor(Rat(1,in)/RtsVec(1,1)),floor(Rat(2,in)/RtsVec(1,2)),0]) ! Moire
               !Rtsp = matmul(celltsp,[floor(Rat(1,i)/RtsVec(1,1)),floor(Rat(2,i)/RtsVec(1,2)),0]) ! Moire
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm-Rts)) 
               Htsp(in,i) = Htsp(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm-Rtsp))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,E,ZWork,lwork,DWork,info) 
         call ZHEEV('V','L',N,Htsp,N,E,ZWork,lwork,DWork,info) 
         !print*, "jk1"
         !print*, Htsp(:,1)
         !print*, "jk2"
         !print*, Hts(:,1)

         do i=1,N
            !Pkcaux = 0.0_dp
            !do i2=1,N
            !do j=1,Nneigh(i) ! maybe go back to full system
            do in=1,N
               !in = NList(j,i)
               !Rts = matmul(cellts,[1,0,0]) ! Moire
               !Rts = matmul(cellts,neighCell(:,j,i)) ! Moire
               !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
               !Rts = matmul(cellts,[floor(Rat(1,in)/RtsVec(1,1)),floor(Rat(2,in)/RtsVec(1,2)),0]) ! Moire
               !Rtsp = matmul(celltsp,[floor(Rat(1,i)/RtsVec(1,1)),floor(Rat(2,i)/RtsVec(1,2)),0]) ! Moire
               Pkcaux = CONJG(Hts(in,i)) * Htsp(in,i)
               !Pkcaux = CONJG(Hts(i,in)) * Htsp(i,in)
               Pkc(i) = Pkc(i) + exp(cmplx_i*dot_product(KG,Rts-Rtsp)) * Pkcaux
            end do
         end do
      end do
   end do
   Nc = 1.0_dp
   nc2 = 5.0_dp*5.0_dp
   Pkc = Nc/nc2 * Pkc
   !print*, "jk7"
   !print*, Pkc
   
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagSpectralWeight')
   end if

         !Pkc = Pkc + exp(cmplx_i*dot_product(KG,Rts-Rtsp)) * (dot_product(CONJG(Hts(:,N)), Htsp(:,N)))  
   ! steps
   !Uts = Hts
   !Utsp = Htsp 
   !print*, "jk1"
   !print*, Htsp(:,1)
   !print*, "jk2"
   !print*, dot_product(CONJG(Hts(:,1)), Htsp(:,1))
   !print*, "jk3"
   !print*, matmul(cellts,[1,1,1])
   !print*, "jk4"
   !print*, matmul(celltsp,[1,1,1])
   !print*, "jk5"
   !Rts = matmul(cellts,[1,1,1])
   !Rtsp = matmul(celltsp,[1,1,1])
   !print*, dot_product(KG,Rts-Rtsp)
   !print*, "jk6"
   !print*, KG

end subroutine DiagSpectralWeightNishi

subroutine DiagSpectralWeightWeiKu(N,ns,is,PkcLoc,E,K,KG,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   !complex(dp), intent(out) :: PkcLoc(N)
   complex(dp), intent(out) :: PkcLoc(N,2)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: E(N)
   real(dp), intent(in) :: K(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   integer :: kk
  
   integer :: cellSize

   integer :: at1, at2
   

   call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   !aBN = aG*1.018181818
   !bncell(:,1) = [aBN,0.0_dp,0.0_dp]
   !bncell(:,2) = [aBN/2.0_dp,sqrt(3.0_dp)*aBN/2.0_dp,0.0_dp]
   !bncell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   cellts = gcell
   celltsp = gcell

   !print*, "cellts", cellts
   !print*, "celltsp", celltsp

   Hts = 0.0_dp
   !Htsp = 0.0_dp
   PkcLoc = 0.0_dp

   !Rts = matmul(cellts,neighCell(:,j,i)) ! Moire
   !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
   !RtsVec(1,:) = [matmul(cellts,[1,0,0])] 
   !RtsVec(2,:) = [matmul(cellts,[0,1,0])] 
   !RtsVec(3,:) = [matmul(cellts,[1,1,0])] 
   !RtsVec(4,:) = [matmul(cellts,[-1,0,0])]
   !RtsVec(5,:) = [matmul(cellts,[0,-1,0])]
   !RtsVec(6,:) = [matmul(cellts,[-1,-1,0])]

   !RtspVec(1,:) = [matmul(cellts,[1,0,0])] 
   !RtspVec(2,:) = [matmul(cellts,[0,1,0])] 
   !RtspVec(3,:) = [matmul(cellts,[1,1,0])] 
   !RtspVec(4,:) = [matmul(cellts,[-1,0,0])]
   !RtspVec(5,:) = [matmul(cellts,[0,-1,0])]
   !RtspVec(6,:) = [matmul(cellts,[-1,-1,0])] 

   !Rtsx = RtsVec(1,1)
   !Rtsy = RtsVec(2,2)
  
   !call MIO_InputParameter('Spectral.numberOfTS',nTS,1)
   call MIO_InputParameter('CellSize', cellSize, 1)

   
   if (.not. frac) call AtomsSetFrac()
   do i=1,N 
       !print*, Rat(1,i), Rat(2,i), floor((Rat(1,i)-0.001)*cellSize), floor(Rat(2,i)*cellSize)
       RtsVec(i,:) = matmul(cellts,[floor((Rat(1,i)-0.001)*cellSize),floor(Rat(2,i)*cellSize),0]) ! Moire
   end do
   !print*, "next"

   !kk = 0
   !do i=0,4 
   !    do j=0,4 
   !        kk = kk+1
   !        RtsVec(kk,:) = matmul(cellts,[i,j,0]) ! Moire
   !        RtspVec(kk,:) = matmul(celltsp,[i,j,0]) ! Moire
   !    end do
   !end do

   !nTS = 6
   
   !print*, "hohoho", RtsVec
   !print*, "hehehe", RtsVec(1,:)
   !do i1=1,nTS!SIZE(RtsVec)
   !   do i2=1,nTS!SIZE(RtspVec)
   !      Rts = RtsVec(i1,:)
   !      Rtsp = RtspVec(i2,:)
         do i=1,N
            Hts(i,i) = H0(i)
            !Htsp(i,i) = H0(i)
            !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
            do j=1,Nneigh(i)
               in = NList(j,i)
               !Rts = matmul(cellts,[floor(Rat(1,in)/RtsVec(1,1)),floor(Rat(2,in)/RtsVec(1,2)),0]) ! Moire
               !Rtsp = matmul(celltsp,[floor(Rat(1,i)/RtsVec(1,1)),floor(Rat(2,i)/RtsVec(1,2)),0]) ! Moire
               Tm = matmul(cell,neighCell(:,j,i))
               !Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm-RtsVec(i,:))) 
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm))
               !Htsp(in,i) = Htsp(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,E,ZWork,lwork,DWork,info) 
         !call ZHEEV('V','L',N,Htsp,N,E,ZWork,lwork,DWork,info) 
         !print*, "jk1"
         !print*, Htsp(:,1)
         !print*, "jk2"
         !print*, Hts(:,1)
       
         !print*, "here", (Hts)

         !do i=1,N
         !   Pkcaux = 0.0_dp
         !   !do i2=1,N
         !   !do j=1,Nneigh(i) ! maybe go back to full system
         !   do in=1,N
         !      !in = i
         !      !in = NList(j,i)
         !      !Rts = matmul(cellts,[1,0,0]) ! Moire
         !      !Rts = matmul(cellts,neighCell(:,j,i)) ! Moire
         !      !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
         !      !Rts = matmul(cellts,[floor(Rat(1,in)/RtsVec(1,1)),floor(Rat(2,in)/RtsVec(1,2)),0]) ! Moire
         !      !Rtsp = matmul(celltsp,[floor(Rat(1,i)/RtsVec(1,1)),floor(Rat(2,i)/RtsVec(1,2)),0]) ! Moire
         !      !Pkcaux = CONJG(Hts(i,in)) * Htsp(i,in)
         !      !Pkcaux = Pkcaux + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(in,i) 
         !      Pkcaux = Pkcaux + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(1,in) 
         !      Pkcaux = Pkcaux + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(3,in) 
         !   !end do
         !   PkcLoc(i) = PkcLoc(i) + Pkcaux
         !end do
         !do i=1,N
         !   Pkcaux = 0.0_dp
         !   Pkcaux1 = 0.0_dp
         !   Pkcaux3 = 0.0_dp
         !   !do i2=1,N
         !   !do j=1,Nneigh(i) ! maybe go back to full system
         !   do in=1,N
         !      !in = i
         !      !in = NList(j,i)
         !      !Rts = matmul(cellts,[1,0,0]) ! Moire
         !      !Rts = matmul(cellts,neighCell(:,j,i)) ! Moire
         !      !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
         !      !Rts = matmul(cellts,[floor(Rat(1,in)/RtsVec(1,1)),floor(Rat(2,in)/RtsVec(1,2)),0]) ! Moire
         !      !Rtsp = matmul(celltsp,[floor(Rat(1,i)/RtsVec(1,1)),floor(Rat(2,i)/RtsVec(1,2)),0]) ! Moire
         !      !Pkcaux = CONJG(Hts(i,in)) * Htsp(i,in)
         !      !if (i.eq.in) then
         !           Pkcaux = Pkcaux + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(i,in) 
         !      !end if
         !      !Pkcaux1 = Pkcaux1 + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(in,1) 
         !      !Pkcaux3 = Pkcaux3 + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(in,3) 
         !   end do
         !   PkcLoc(i) = PkcLoc(i) + Pkcaux
         !   !PkcLoc(1) = PkcLoc(1) + Pkcaux1
         !   !PkcLoc(2) = PkcLoc(2) + Pkcaux3
         !end do
         !do i=1,2
         !print*, Species(in), Species(22), Species(21)
         !if (frac) call AtomsSetCart()
         !do i=1,N
         !   if (Rat(1,i).lt.1.5_dp) then
         !      at1 = i
         !   else if (Rat(2,i).lt.1.5_dp .and. Rat(1,i).lt.2.5) then
         !      at2 = i
         !   end if 
         !end do
         !print*, at1, at2
         do j=1,N ! These are the eigenvectors with band index J
             do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
                 if (Species(in).eq.1) then
                     PkcLoc(j,1) = PkcLoc(j,1) + exp(-cmplx_i*dot_product(KG, RtsVec(in,:))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
                     !print*, "hi1"
                 else if (Species(in).eq.2) then
                     PkcLoc(j,2) = PkcLoc(j,2) + exp(-cmplx_i*dot_product(KG, RtsVec(in,:))) * Hts(in,j) 
                     !print*, "hi3"
                 end if
             end do
         end do
   !   end do
   !end do
   !print*, "PkcLoc(1) ", PkcLoc(1)
   !print*, "PkcLoc(2) ", PkcLoc(2)
   Nc = 1.0
   nc2 = 1.0
   PkcLoc = Nc/nc2 * PkcLoc
   !Pkc = abs(Pkc)**2.0_dp
   !print*, Pkc
   !print*, "jk7"
   !print*, Pkc
   
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagSpectralWeight')
   end if

         !Pkc = Pkc + exp(cmplx_i*dot_product(KG,Rts-Rtsp)) * (dot_product(CONJG(Hts(:,N)), Htsp(:,N)))  
   ! steps
   !Uts = Hts
   !Utsp = Htsp 
   !print*, "jk1"
   !print*, Htsp(:,1)
   !print*, "jk2"
   !print*, dot_product(CONJG(Hts(:,1)), Htsp(:,1))
   !print*, "jk3"
   !print*, matmul(cellts,[1,1,1])
   !print*, "jk4"
   !print*, matmul(celltsp,[1,1,1])
   !print*, "jk5"
   !Rts = matmul(cellts,[1,1,1])
   !Rtsp = matmul(celltsp,[1,1,1])
   !print*, dot_product(KG,Rts-Rtsp)
   !print*, "jk6"
   !print*, KG

end subroutine DiagSpectralWeightWeiKu

subroutine DiagSpectralWeightWeiKuInequivalentOld(N,ns,is,PkcLocA,PkcLocB,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLocA(N,3), PkcLocB(N,3)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   !integer :: kk
   real(dp) :: ll, kk
  
   integer :: cellSize

   integer :: at1, at2

   !integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)
   
   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)
  
   real(dp) :: topBottomRatio
 
   logical :: changeExpSign


   cellts = gcell ! Lattice vectors of PC (either top or bottom layer)
   !celltsp = gcell

   Hts = 0.0_dp
   PkcLocA = 0.0_dp
   PkcLocB = 0.0_dp
   if (frac) call AtomsSetCart()
       G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
       G01 = matmul(cellts,[0,1,0])
       G11 = matmul(cellts,[1,1,0])
       do i=1,N 
           ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
           kk = (Rat(1,i) - ll * G01(1)) / G10(1)
           if (kk.gt.0) then
               kk = floor(kk)
           else
               kk = ceiling(kk)
           end if
           if (ll.gt.0) then
               ll = floor(ll)
           else
               ll = ceiling(ll)
           end if
           !print*, Rat(1,i), Rat(2,i), kk, ll
           !RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
           RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       end do
         do i=1,N
            Hts(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
         if (info/=0) then
            call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
         end if
         ! Note that I'm using the Wei Ku expression for the spectral function
         ! (this is NOT the NISHI one)
!         do j=1,N ! These are the eigenvectors with band index J
!             do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
!                 if (layerIndex(in).eq.1) then
!                     PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!                 else if (layerIndex(in).eq.2) then
!                     PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!                 else if (layerIndex(in).eq.3) then
!                     PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * topBottomRatio ! to account for experimental weight between top and bottom ARPES
!                 end if
!             end do
!         end do
         !call MIO_InputParameter('changeExpSign',changeExpSign,.true.)
         !changeExpSign = .true.
         !if (changeExpSign) then
             do j=1,N
                 do in=1,N
                   if (Species(in).eq.1) then
                     if (layerIndex(in).eq.1) then
                         PkcLocA(j,1) = PkcLocA(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocA(j,2) = PkcLocA(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocA(j,3) = PkcLocA(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.2) then
                     if (layerIndex(in).eq.1) then
                         PkcLocB(j,1) = PkcLocB(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocB(j,2) = PkcLocB(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocB(j,3) = PkcLocB(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     end if
                   end if
                 end do
                 !print*, "here3", abs(PkcLoc(j,1))
                 !print*, "here4", abs(PkcLoc(j,2))
             end do
         !else
         !    do j=1,N
         !        do in=1,N
         !            if (layerIndex(in).eq.1) then
         !                PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            else if (layerIndex(in).eq.2) then
         !                PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            else if (layerIndex(in).eq.3) then
         !                PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            end if
         !        end do
         !    end do
         !end if
   Nc = 1.0
   nc2 = 1.0
   PkcLocA = Nc/nc2 * PkcLocA
   PkcLocB = Nc/nc2 * PkcLocB

end subroutine DiagSpectralWeightWeiKuInequivalentOld

subroutine DiagSpectralWeightWeiKuInequivalent(N,ns,is,PkcLocA,PkcLocB,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLocA(N,4), PkcLocB(N,4)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   !integer :: kk
   real(dp) :: ll, kk
  
   integer :: cellSize

   integer :: at1, at2

   !integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)
   
   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)
  
   real(dp) :: topBottomRatio
 
   logical :: changeExpSign


   cellts = gcell ! Lattice vectors of PC (either top or bottom layer)
   !celltsp = gcell

   Hts = 0.0_dp
   PkcLocA = 0.0_dp
   PkcLocB = 0.0_dp
   if (frac) call AtomsSetCart()
       G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
       G01 = matmul(cellts,[0,1,0])
       G11 = matmul(cellts,[1,1,0])
       do i=1,N 
           ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
           kk = (Rat(1,i) - ll * G01(1)) / G10(1)
           if (kk.gt.0) then
               kk = floor(kk)
           else
               kk = ceiling(kk)
           end if
           if (ll.gt.0) then
               ll = floor(ll)
           else
               ll = ceiling(ll)
           end if
           !print*, Rat(1,i), Rat(2,i), kk, ll
           !RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
           RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       end do
         do i=1,N
            Hts(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
         if (info/=0) then
            call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
         end if
         ! Note that I'm using the Wei Ku expression for the spectral function
         ! (this is NOT the NISHI one)
!         do j=1,N ! These are the eigenvectors with band index J
!             do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
!                 if (layerIndex(in).eq.1) then
!                     PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!                 else if (layerIndex(in).eq.2) then
!                     PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!                 else if (layerIndex(in).eq.3) then
!                     PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * topBottomRatio ! to account for experimental weight between top and bottom ARPES
!                 end if
!             end do
!         end do
         !call MIO_InputParameter('changeExpSign',changeExpSign,.true.)
         !changeExpSign = .true.
         !if (changeExpSign) then
             do j=1,N
                 do in=1,N
                   !if (Species(in).eq.1) then
                   !  if (layerIndex(in).eq.1) then
                   !      PkcLocA(j,1) = PkcLocA(j,1) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.2) then
                   !      PkcLocA(j,2) = PkcLocA(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.3) then
                   !      PkcLocA(j,3) = PkcLocA(j,3) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.4) then
                   !      PkcLocA(j,4) = PkcLocA(j,4) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  end if
                   !else if (Species(in).eq.2) then
                   !  if (layerIndex(in).eq.1) then
                   !      PkcLocB(j,1) = PkcLocB(j,1) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.2) then
                   !      PkcLocB(j,2) = PkcLocB(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.3) then
                   !      PkcLocB(j,3) = PkcLocB(j,3) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.4) then
                   !      PkcLocB(j,4) = PkcLocB(j,4) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  end if
                   !end if
                   if (Species(in).eq.1) then
                     if (layerIndex(in).eq.1) then
                         PkcLocA(j,1) = PkcLocA(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocA(j,2) = PkcLocA(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocA(j,3) = PkcLocA(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocA(j,4) = PkcLocA(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.2) then
                     if (layerIndex(in).eq.1) then
                         PkcLocB(j,1) = PkcLocB(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocB(j,2) = PkcLocB(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocB(j,3) = PkcLocB(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocB(j,4) = PkcLocB(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   end if
                 end do
                 !print*, "here3", abs(PkcLoc(j,1))
                 !print*, "here4", abs(PkcLoc(j,2))
             end do
         !else
         !    do j=1,N
         !        do in=1,N
         !            if (layerIndex(in).eq.1) then
         !                PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            else if (layerIndex(in).eq.2) then
         !                PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            else if (layerIndex(in).eq.3) then
         !                PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            end if
         !        end do
         !    end do
         !end if
   Nc = 1.0
   nc2 = 1.0
   PkcLocA = Nc/nc2 * PkcLocA
   PkcLocB = Nc/nc2 * PkcLocB

end subroutine DiagSpectralWeightWeiKuInequivalent

subroutine DiagSpectralWeightWeiKuInequivalentMoreOrbitals(N,ns,is,PkcLocA,PkcLocB,PkcLocC,PkcLocD,PkcLocE,PkcLocF,PkcLocG,PkcLocH,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLocA(N,4), PkcLocB(N,4)
   complex(dp), intent(out) :: PkcLocC(N,4), PkcLocD(N,4)
   complex(dp), intent(out) :: PkcLocE(N,4), PkcLocF(N,4)
   complex(dp), intent(out) :: PkcLocG(N,4), PkcLocH(N,4)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   !integer :: kk
   real(dp) :: ll, kk
  
   integer :: cellSize

   integer :: at1, at2

   !integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)
   
   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)
  
   real(dp) :: topBottomRatio
 
   logical :: changeExpSign


   cellts = gcell ! Lattice vectors of PC (either top or bottom layer)
   !celltsp = gcell

   Hts = 0.0_dp
   PkcLocA = 0.0_dp
   PkcLocB = 0.0_dp
   PkcLocC = 0.0_dp
   PkcLocD = 0.0_dp
   PkcLocE = 0.0_dp
   PkcLocF = 0.0_dp
   PkcLocG = 0.0_dp
   PkcLocH = 0.0_dp
   if (frac) call AtomsSetCart()
       G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
       G01 = matmul(cellts,[0,1,0])
       G11 = matmul(cellts,[1,1,0])
       do i=1,N 
           ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
           kk = (Rat(1,i) - ll * G01(1)) / G10(1)
           if (kk.gt.0) then
               kk = floor(kk)
           else
               kk = ceiling(kk)
           end if
           if (ll.gt.0) then
               ll = floor(ll)
           else
               ll = ceiling(ll)
           end if
           !print*, Rat(1,i), Rat(2,i), kk, ll
           !RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
           RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       end do
         do i=1,N
            Hts(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
         if (info/=0) then
            call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
         end if
         ! Note that I'm using the Wei Ku expression for the spectral function
         ! (this is NOT the NISHI one)
!         do j=1,N ! These are the eigenvectors with band index J
!             do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
!                 if (layerIndex(in).eq.1) then
!                     PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!                 else if (layerIndex(in).eq.2) then
!                     PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!                 else if (layerIndex(in).eq.3) then
!                     PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * topBottomRatio ! to account for experimental weight between top and bottom ARPES
!                 end if
!             end do
!         end do
         !call MIO_InputParameter('changeExpSign',changeExpSign,.true.)
         !changeExpSign = .true.
         !if (changeExpSign) then
             do j=1,N
                 do in=1,N
                   !if (Species(in).eq.1) then
                   !  if (layerIndex(in).eq.1) then
                   !      PkcLocA(j,1) = PkcLocA(j,1) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.2) then
                   !      PkcLocA(j,2) = PkcLocA(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.3) then
                   !      PkcLocA(j,3) = PkcLocA(j,3) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.4) then
                   !      PkcLocA(j,4) = PkcLocA(j,4) + exp(-cmplx_i*dot_product(KG, Rat(:,1)-Rat(:,in))) * Hts(in,j)
                   !  end if
                   !else if (Species(in).eq.2) then
                   !  if (layerIndex(in).eq.1) then
                   !      PkcLocB(j,1) = PkcLocB(j,1) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.2) then
                   !      PkcLocB(j,2) = PkcLocB(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.3) then
                   !      PkcLocB(j,3) = PkcLocB(j,3) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  else if (layerIndex(in).eq.4) then
                   !      PkcLocB(j,4) = PkcLocB(j,4) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
                   !  end if
                   !end if
                   if (Species(in).eq.1) then
                     if (layerIndex(in).eq.1) then
                         PkcLocA(j,1) = PkcLocA(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocA(j,2) = PkcLocA(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocA(j,3) = PkcLocA(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocA(j,4) = PkcLocA(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.2) then
                     if (layerIndex(in).eq.1) then
                         PkcLocB(j,1) = PkcLocB(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocB(j,2) = PkcLocB(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocB(j,3) = PkcLocB(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocB(j,4) = PkcLocB(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.3) then
                         PkcLocC(j,2) = PkcLocC(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.4) then
                         PkcLocD(j,2) = PkcLocD(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.5) then
                         PkcLocE(j,2) = PkcLocE(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.6) then
                         PkcLocF(j,2) = PkcLocF(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.7) then
                         PkcLocG(j,2) = PkcLocG(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.8) then
                         PkcLocH(j,2) = PkcLocH(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   end if
                 end do
                 !print*, "here3", abs(PkcLoc(j,1))
                 !print*, "here4", abs(PkcLoc(j,2))
             end do
         !else
         !    do j=1,N
         !        do in=1,N
         !            if (layerIndex(in).eq.1) then
         !                PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            else if (layerIndex(in).eq.2) then
         !                PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            else if (layerIndex(in).eq.3) then
         !                PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
         !            end if
         !        end do
         !    end do
         !end if
   Nc = 1.0
   nc2 = 1.0
   PkcLocA = Nc/nc2 * PkcLocA
   PkcLocB = Nc/nc2 * PkcLocB

end subroutine DiagSpectralWeightWeiKuInequivalentMoreOrbitals

subroutine DiagSpectralWeightWeiKuInequivalentLee(N,ns,is,PkcLoc1,PkcLoc2,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh
   use math

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLoc1(N,3)
   complex(dp), intent(out) :: PkcLoc2(N,3)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: TmVec(N,3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy
   real(dp) :: RatG(3,2)

   real(dp) :: Nc, nc2

   integer :: nTS
   !integer :: kk
   real(dp) :: ll, kk
  
   integer :: cellSize

   integer :: at1, at2

   !integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)
   
   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)
  
   real(dp) :: topBottomRatio
 
   logical :: changeExpSign
   
   integer :: ix, iy, ncell(3,9)
   real(dp) :: rMinRp(3), v(3), dmin


   cellts = cell ! Lattice vectors of PC (either top or bottom layer)
   celltsp = gcell

   Hts = 0.0_dp
   PkcLoc1 = 0.0_dp
   PkcLoc2 = 0.0_dp
   if (frac) call AtomsSetCart()
   G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
   G01 = matmul(cellts,[0,1,0])
   G11 = matmul(cellts,[1,1,0])
   G10p = matmul(celltsp,[1,0,0]) ! here, its in real space!!!
   G01p = matmul(celltsp,[0,1,0])
   G11p = matmul(celltsp,[1,1,0])
   do i=1,N 
       ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
       kk = (Rat(1,i) - ll * G01(1)) / G10(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       !print*, Rat(1,i), Rat(2,i), kk, ll
       !RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire

       ll = (Rat(1,i)*G10p(2)/G10p(1) - Rat(2,i)) / (G01p(1)*G10p(2)/G10p(1) - G01p(2)) 
       kk = (Rat(1,i) - ll * G01p(1)) / G10p(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       RtspVec(i,:) = matmul(celltsp,[int(kk),int(ll),0]) ! Graphene
   end do
   do i=1,N
      Hts(i,i) = H0(i)
      do j=1,Nneigh(i)
         in = NList(j,i)
         Tm = matmul(cell,neighCell(:,j,i))
         Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
      end do
   end do
   call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
   end if
   ! Note that I'm using the Wei Ku expression for the spectral function
   ! (this is NOT the NISHI one)
!   do j=1,N ! These are the eigenvectors with band index J
!       do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
!           if (layerIndex(in).eq.1) then
!               PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!           else if (layerIndex(in).eq.2) then
!               PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
!           else if (layerIndex(in).eq.3) then
!               PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * topBottomRatio ! to account for experimental weight between top and bottom ARPES
!           end if
!       end do
!   end do
   !call MIO_InputParameter('changeExpSign',changeExpSign,.true.)
   !changeExpSign = .true.
   !if (changeExpSign) then
   !RatG(:,1) = Rat(:,543)
   !RatG(:,2) = Rat(:,544)
   RatG(:,1) = Rat(:,1)
   RatG(:,2) = Rat(:,2)
   !print*, Rat(:,10)
   do j=1,N
       do in=1,N
           !dmin = 10000.0_dp
           !do ix=-1,1; do iy=-1,1
           !   !if (i==j .and. ix==0 .and. iy==0) cycle
           !   ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
           !   if (Species(in).eq.1) then
           !      v = Rat(:,1) - Rat(:,in) ! - matmul(cell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition. 
           !   else
           !      v = Rat(:,2) - Rat(:,in) !- matmul(cell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition. 
           !   end if
           !   if (norm(v) < dmin) then
           !       rMinRp = v
           !       dmin = norm(v)
           !   end if
           !end do; end do
           
           if (Species(in).eq.1) then
              rMinRp = Rat(:,1)-Rat(:,in)
           else
              rMinRp = Rat(:,2)-Rat(:,in)
           end if
           !rMinRp = -TmVec(in,:)
           !rMinRp = RtsVec - RtspVec(in,:)
           if (layerIndex(in).eq.1) then
               !if (mod(in,2)>0) then
               !   PkcLoc1(j,1) = PkcLoc1(j,1) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
               !else
               !   PkcLoc2(j,1) = PkcLoc2(j,1) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
               !end if
               !PkcLoc(j,1) = PkcLoc(j,1) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
               if (Species(in).eq.1) then
                  PkcLoc1(j,1) = PkcLoc1(j,1) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
                  !PkcLocB(j,2) = PkcLocB(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,2)-Rat(:,in))) * Hts(in,j)
               else
                  PkcLoc2(j,1) = PkcLoc2(j,1) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               end if
           else if (layerIndex(in).eq.2) then
               !if (mod(in,2)>0) then
               !   PkcLoc1(j,2) = PkcLoc1(j,2) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
               !else
               !   PkcLoc2(j,2) = PkcLoc2(j,2) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
               !end if
               !PkcLoc1(j,2) = PkcLoc1(j,2) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               !PkcLoc(j,2) = PkcLoc(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * conjg(Hts(in,j))
               !PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
               if (Species(in).eq.1) then
                  PkcLoc1(j,2) = PkcLoc1(j,2) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               else
                  PkcLoc2(j,2) = PkcLoc2(j,2) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               end if
           else if (layerIndex(in).eq.3) then
               !if (mod(in,2)>0) then
               !   PkcLoc1(j,3) = PkcLoc1(j,3) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
               !else
               !   PkcLoc2(j,3) = PkcLoc2(j,3) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
               !end if
               !PkcLoc1(j,3) = PkcLoc1(j,3) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               !PkcLoc(j,3) = PkcLoc(j,3) + exp(-cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * conjg(Hts(in,j))
               !PkcLoc(j,3) = PkcLoc(j,3) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
               if (Species(in).eq.1) then
                  PkcLoc1(j,3) = PkcLoc1(j,3) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               else
                  PkcLoc2(j,3) = PkcLoc2(j,3) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               end if
           end if
           !if (layerIndex(in).eq.1) then
           !    if (mod(in,2)>0) then
           !       !PkcLoc1(j,1) = PkcLoc1(j,1) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !       PkcLoc1(j,1) = PkcLoc1(j,1) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !    else
           !       !PkcLoc2(j,1) = PkcLoc2(j,1) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !       PkcLoc2(j,1) = PkcLoc2(j,1) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !    end if
           !    !PkcLoc(j,1) = PkcLoc(j,1) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !else if (layerIndex(in).eq.2) then
           !    if (mod(in,2)>0) then
           !       !PkcLoc1(j,2) = PkcLoc1(j,2) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !       PkcLoc1(j,2) = PkcLoc1(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !    else
           !       !PkcLoc2(j,2) = PkcLoc2(j,2) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !       PkcLoc2(j,2) = PkcLoc2(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !    end if
           !    !PkcLoc(j,2) = PkcLoc(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * conjg(Hts(in,j))
           !    !PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !else if (layerIndex(in).eq.3) then
           !    if (mod(in,2)>0) then
           !       PkcLoc1(j,3) = PkcLoc1(j,3) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !       !PkcLoc1(j,3) = PkcLoc1(j,3) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !    else
           !       PkcLoc2(j,3) = PkcLoc2(j,3) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !       !PkcLoc2(j,3) = PkcLoc2(j,3) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !    end if
           !    !PkcLoc(j,3) = PkcLoc(j,3) + exp(-cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * conjg(Hts(in,j))
           !    !PkcLoc(j,3) = PkcLoc(j,3) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !end if
       end do
       !print*, "here3", abs(PkcLoc(j,1))
       !print*, "here4", abs(PkcLoc(j,2))
   end do
   !else
   !    do j=1,N
   !        do in=1,N
   !            if (layerIndex(in).eq.1) then
   !                PkcLoc(j,1) = PkcLoc(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
   !            else if (layerIndex(in).eq.2) then
   !                PkcLoc(j,2) = PkcLoc(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
   !            else if (layerIndex(in).eq.3) then
   !                PkcLoc(j,3) = PkcLoc(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
   !            end if
   !        end do
   !    end do
   !end if
   Nc = 1.0
   nc2 = 1.0
   PkcLoc1 = Nc/nc2 * PkcLoc1
   PkcLoc2 = Nc/nc2 * PkcLoc2

end subroutine DiagSpectralWeightWeiKuInequivalentLee

subroutine DiagSpectralWeightWeiKuInequivalentNishi(N,ns,is,PkcLoc,ELoc1, ELoc2,KptsLoc,KG,cell,gcell1,gcell2,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   !complex(dp), intent(out) :: PkcLoc(N)
   !complex(dp), intent(out) :: Pkc(N,2)
   complex(dp), intent(out) :: PkcLoc(N,2)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc1(N), ELoc2(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell1(3,3)
   real(dp) :: gcell2(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info, inn
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   !integer :: kk
   real(dp) :: ll, kk
  
   integer :: cellSize

   integer :: at1, at2

   !integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)
   
   complex(dp) :: ZWorkLoc1(lwork)
   complex(dp) :: ZWorkLoc2(lwork)
   real(dp) :: DWorkLoc1(3*N-2)
   real(dp) :: DWorkLoc2(3*N-2)


   cellts = gcell1
   celltsp = gcell2

   Hts = 0.0_dp
   Htsp = 0.0_dp
   PkcLoc = 0.0_dp
   if (frac) call AtomsSetCart()
   G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
   G01 = matmul(cellts,[0,1,0])
   G11 = matmul(cellts,[1,1,0])
   G10p = matmul(celltsp,[1,0,0]) ! here, its in real space!!!
   G01p = matmul(celltsp,[0,1,0])
   G11p = matmul(celltsp,[1,1,0])
   do i=1,N 
       ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
       kk = (Rat(1,i) - ll * G01(1)) / G10(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       !print*, Rat(1,i), Rat(2,i), kk, ll
       !RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire

       ll = (Rat(1,i)*G10p(2)/G10p(1) - Rat(2,i)) / (G01p(1)*G10p(2)/G10p(1) - G01p(2)) 
       kk = (Rat(1,i) - ll * G01p(1)) / G10p(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       RtspVec(i,:) = matmul(celltsp,[int(kk),int(ll),0]) ! Moire
   end do
   do i=1,N
      Hts(i,i) = H0(i)
      Htsp(i,i) = H0(i)
      do j=1,Nneigh(i)
         in = NList(j,i)
         Tm = matmul(cell,neighCell(:,j,i)) ! Corresponds to L in Nishi equation
         Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm-RtsVec(i,:))) ! RtsVec is the multiple of ts
         Htsp(in,i) = Htsp(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm-RtspVec(i,:)))
      end do
   end do
   call ZHEEV('V','L',N,Hts,N,ELoc1,ZWorkLoc1,lwork,DWorkLoc1,info) 
   call ZHEEV('V','L',N,Htsp,N,ELoc2,ZWorkLoc2,lwork,DWorkLoc2,info) 
   if (info/=0) then
      call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
   end if
   ! This one is NISHI
   do j=1,N ! These are the eigenvectors with band index J
       do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
           do inn=1,N
               !if (Species(in).eq.1 .and. Species(inn).eq.1 ) then
                  if (layerIndex(in).eq.1 .and. layerIndex(inn).eq.1) then
                      PkcLoc(j,1) = PkcLoc(j,1) + exp(cmplx_i*dot_product(KG, (RtsVec(in,:)-RtsVec(inn,:)))) * conjg(Hts(in,j)) * Hts(inn,j)! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
                      !print*, "hi1"
                  else if (layerIndex(in).eq.2 .and. layerIndex(inn).eq.2) then
                      PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtspVec(inn,:)))) * conjg(Htsp(in,j)) * Htsp(inn,j)! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
                      !PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, RtsVec(in,:))) * Hts(in,j) 
                      !print*, "hi3"
                  end if
               !else if (Species(in).eq.2 .and. Species(inn).eq.2) then
               !   if (layerIndex(in).eq.1 .and. layerIndex(inn).eq.1) then
               !       PkcLoc(j,1) = PkcLoc(j,1) + exp(cmplx_i*dot_product(KG, (RtsVec(in,:)-RtsVec(inn,:)))) * conjg(Hts(in,j)) * Hts(inn,j)! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
               !       !print*, "hi1"
               !   else if (layerIndex(in).eq.2 .and. layerIndex(inn).eq.2) then
               !       PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtspVec(inn,:)))) * conjg(Htsp(in,j)) * Htsp(inn,j)! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
               !       !PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, RtsVec(in,:))) * Hts(in,j) 
               !       !print*, "hi3"
               !   end if
               !end if
           end do
       end do
   end do
   Nc = 1.0
   nc2 = 1.0
   PkcLoc = Nc/nc2 * PkcLoc

end subroutine DiagSpectralWeightWeiKuInequivalentNishi

function convolve(x, h, Epts)
    implicit none
    
    integer :: Epts
    !x is the signal array
    !h is the noise/impulse array
    !real, dimension(:), allocatable :: convolve, y
    !real, dimension(:) :: x, h
    real(dp), allocatable :: convolve(:), y(:)
    real(dp) :: x(Epts), h(Epts)
    !complex(dp), allocatable :: convolve, y
    !complex(dp), dimension(:) :: x, h
    integer :: kernelsize, datasize
    integer :: i,j,k
    
    datasize = size(x)  ! Ake
    kernelsize = size(h) ! gaussian


    
    allocate(y(datasize))
    allocate(convolve(datasize))
   
    !last part
    do i=kernelsize,datasize
        y(i) = 0.0
        j=i
        do k=1,kernelsize
            y(i) = y(i) + x(j)*h(k)
            j = j-1
        end do
    end do
    
    !first part
    do i=1,kernelsize
        y(i) = 0.0
        j=i
        k=1
        do while (j > 0)
            y(i) = y(i) + x(j)*h(k)
            j = j-1
            k = k+1
        end do
    end do
    
    convolve = y
           
end function convolve

pure function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(dp), intent(in) :: A(3,3)   !! Matrix
    real(dp)             :: B(3,3)   !! Inverse matrix
    real(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function

subroutine av (n, v, w)
      integer  ::         n, j
      complex (dp) :: v(n), w(n), rho,  dd, dl, du, s, h, h2
      complex(dp), parameter ::   one = (1.0D+0, 0.0D+0) , two = (2.0D+0, 0.0D+0) 
      common            /convct/ rho

      h = one / dcmplx (n+1)
      h2 = h*h
      s = rho / two
      dd = two / h2
      dl = -one/h2 - s/h
      du = -one/h2 + s/h

      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(n) =  dl*v(n-1) + dd*v(n)
      return
      end


end module diag
