! generator of slip rate functions for the
! ruiz integral k-squared (rik) source model
! (i.e. modified model by ruiz et al., 2011, see gallovic, 2016)
! with random rupture front
! coded by frantisek gallovic (2016)
! version 2.0
!-------------------------------------------------------------------------------

module ruptveloc
    implicit none
    integer nlayers
    real, allocatable:: hvr(:),vr(:)   ! layer-top location (from bottom), rupture velocities inside layers
    real hl,hw                         ! location of hypocenter
end module

module crustaldat
    implicit none
    integer ndepth
    real,allocatable,dimension(:):: depth,vp,vs,rho
    real hypodepth,dip
end module

program kkstf
    use ruptveloc
    use crustaldat
    use mpi
    implicit none
    real,parameter:: pi=3.1415926535
    integer nl,nw
    real mufix
    !pdf for subsource position
    integer,parameter:: pdfnl=2000,pdfnw=1000
    real,allocatable:: pdf2d(:,:),cpdf2d(:,:)
    real pdfdl,pdfdw,pdfgaussl,pdfgaussw,pdfgausss
    real,allocatable:: ruptimegen(:)
    character*256 filename,inputfile
    integer ml(2),pdfoption,filenl,filenw
    !subsource parameters:
    integer sroption,submax,submin
    real,allocatable,dimension(:):: subposl,subposw,subsize,subslip,submoment,subrisetime,submvr,subruptime,subnucll,subnuclw
    integer,allocatable,dimension(:):: subno
    integer subtot,idum1,idum2
    !source parameters:
    real lf,wf,l,w,sml,smw,lwratio,m0,l0,vrsubfact,aparam,dt
    real,allocatable,dimension(:):: srl,srw,srelem,srmu,srslip,srmoment,srstressdrop,sr,stf
    ! [modify]
    real,allocatable,dimension(:):: outsrmoment
    integer nsr,nt
    !others
    real ran2,dum,duml,dumw,dumphi,dumr,totmoment,time,meanvr,ruptime,ruptimesr,momentcorr
    integer i,j,k,m,hits
    integer ierr,rank,nprocs

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

    read(5,*) inputfile
    open(101,file=trim(inputfile))
    write(6,*) "reading"
    read(101,*)
    read(101,*)lf, wf
    read(101,*)
    read(101,*)l, w, sml,smw
    if(sml+l>lf.or.smw+w>wf)then
        write(*,*)'strong-motion area exceeds the fault size!'
        stop
    endif
    read(101,*)
    read(101,*)m0
    read(101,*)
    read(101,*)sroption
    if(sroption==1)then
        read(101,*)nl,nw,mufix
        nsr=nl*nw
        if(mufix==0.)call read_crustaldat()
    elseif(sroption==2)then
        read(101,*)nl,nw,nsr
    else
        write(*,*)'wrong sroption!'
        stop
    endif
    read(101,*)
    read(101,*)hl,hw
    read(101,*)
    read(101,*)hypodepth,dip
    read(101,*)
    read(101,*)l0,vrsubfact,aparam
    read(101,*)
    read(101,*)submin,submax
    read(101,*)
    read(101,*)pdfoption
    if(pdfoption==2)read(101,*)pdfgaussl,pdfgaussw,pdfgausss
    if(pdfoption==3)read(101,*)filenl,filenw,filename
    read(101,*)
    read(101,*)idum1,idum2
    read(101,*)
    read(101,*)dt,nt
    read(101,*)
    read(101,*)nlayers
    if(nlayers>0)then
        allocate(hvr(nlayers),vr(nlayers))
        do i=1,nlayers
        read(101,*)hvr(i),vr(i)
        enddo
    else
        call read_crustaldat()
        read(101,*)dum   !constant vr/vs
        nlayers=ndepth
        allocate(hvr(nlayers),vr(nlayers))
        do i=1,nlayers
        hvr(i)=(hypodepth-depth(nlayers-i+1))/sin(dip/180.d0*pi)+hw
        vr(i)=vs(nlayers-i+1)*dum
        enddo
    endif
    do i=2,nlayers
    if(vr(i)>vr(i-1))then
        write(*,*)'error! vr should decrease upwards!'
    endif
    enddo
    if(hvr(nlayers)<wf)then
        write(*,*)'error! definition of vr does not cover the whole fault!'
        stop
    endif
    close(101)

    open(232,file='nucleationpoint.dat')
    write(232,*)hl,hw
    close(232)
    lwratio=l/w
    l0=l0*w
    
    ! distribute grid point over processors
    call splitgrid(rank,NL,NW)
    !preparing pdf for distribution of subsources
    write(*,*)'preparing pdf for subsource distribution...'
    allocate(pdf2d(pdfnl,pdfnw),cpdf2d(pdfnl,pdfnw))
    call fillpdf(pdfnl,pdfnw,lf,wf,l,w,sml,smw,pdf2d,cpdf2d,pdfoption,filename,filenl,filenw,pdfgaussl,pdfgaussw,pdfgausss)
    pdfdl=lf/real(pdfnl)
    pdfdw=wf/real(pdfnw)
    write(*,*)'... done.'

    !calculate and read random variations of rupture time
    call randomruptvel(lf,wf,nl,nw,hl,hw,idum2)
    allocate(ruptimegen(nw*nl))
    open(232,file='ruptimegen.txt')
    k=0
    do j=1,nw
    do i=1,nl
    k=k+1
    read(232,*)dum,dum,ruptimegen(k)
    enddo
    read(232,*)
    enddo
    close(232)

    !reading location of the slip rate points
    allocate(srl(nsr),srw(nsr),srelem(nsr),srmu(nsr),srslip(nsr),srmoment(nsr),srstressdrop(nsr))
    if(sroption==1)then
        srelem=lf*wf/real(nl*nw)
        k=0
        do j=1,nw
        do i=1,nl
        k=k+1
        srl(k)=lf/real(nl)*(real(i)-.5)
        srw(k)=wf/real(nw)*(real(j)-.5)
        enddo
        enddo
        if(mufix>0.)then
            srmu=mufix
        else
            do k=1,nsr
            dum=(hypodepth+(hw-srw(k))*sin(dip/180.d0*pi))
            if(dum>depth(ndepth))then
                srmu(k)=rho(ndepth)*vs(ndepth)**2*1.d9
            else
                do j=1,ndepth
                if(dum<depth(j))exit
                enddo
                srmu(k)=rho(j-1)*vs(j-1)**2*1.d9
            endif  
            enddo
        endif
    else
        open(101,file='srloc.dat')
        do i=1,nsr
        read(101,*)srl(i),srw(i),srelem(i),srmu(i)
        enddo
    endif
    srelem=srelem*1.e6

    !static subsource parameters
    write(*,*)'preparing distribution and parameters of subsources...'
    allocate(subno(submax))
    subno=0
    do i=submin,submax
    subno(i)=int(float(2*i-1)*lwratio)
    enddo
    subtot=sum(subno)
    allocate(subposl(subtot),subposw(subtot),subsize(subtot),subslip(subtot),&
        submoment(subtot),subrisetime(subtot))
    allocate(subnucll(subtot),subnuclw(subtot),submvr(subtot),subruptime(subtot))
    k=0
    srslip=0.
    srmoment=0.
    srstressdrop=0.
    subslip=0.
    submoment=0.
    do i=submin,submax
    do j=1,subno(i)
    k=k+1
    subsize(k)=w/float(i)/2.   !radius

    ! locating the subsource according to the pdf
    do
    !          if(k==1.and.pdfoption.ne.1)then   !put the largest subsource at the position of the pdf maximum
    !            ml=maxloc(pdf2d(:,:))
    !          else
    ml=minloc(abs(cpdf2d(:,:)-ran2(idum1)))
    !          endif
    subposl(k)=(real(ml(1))-.5)*pdfdl

    if(submin==1.and.i==1)then
        subposw(k)=w/2.
        if(subposl(k)-subsize(k)>=0..and.subposl(k)+subsize(k)<=lf)exit
    else
        subposw(k)=(real(ml(2))-.5)*pdfdw
        if(subposl(k)-subsize(k)>=0..and.subposw(k)-subsize(k)>=0..and.&
            subposl(k)+subsize(k)<=lf.and.subposw(k)+subsize(k)<=wf)exit
    endif
    enddo

    do m=1,nsr
    dum=subsize(k)**2-(srl(m)-subposl(k))**2-(srw(m)-subposw(k))**2
    if(dum>0.)then
        submoment(k)=submoment(k)+sqrt(dum)*srelem(m)*srmu(m)
        subslip(k)=subslip(k)+sqrt(dum)
        srmoment(m)=srmoment(m)+sqrt(dum)*srelem(m)*srmu(m)
        srslip(m)=srslip(m)+sqrt(dum)
        srstressdrop(m)=srstressdrop(m)-srmu(m)/24.*7.*pi/1000. !(oprava na to, ze se pracuje se skluzem v km)
    endif
    enddo
    enddo
    enddo
    totmoment=sum(srmoment(:))
    momentcorr=m0/totmoment
    write(*,*)totmoment*momentcorr
    subslip=subslip*momentcorr
    submoment=submoment*momentcorr
    srslip=srslip*momentcorr
    srmoment=srmoment*momentcorr
    srstressdrop=srstressdrop*momentcorr

    open(201,file='slipdistribution.dat')
    do i=1,nsr
    write(201,'(10e13.5)')srl(i),srw(i),srslip(i),srmoment(i),ruptimegen(i),srstressdrop(i)
    enddo
    close(201)
    open(201,file='slipdistribution.gnuplot.dat')
    do i=1,nw
    write(201,'(1000e13.5)')srslip((i-1)*nl+1:i*nl)
    enddo
    write(201,*);write(201,*)
    do i=1,nw
    write(201,'(1000e13.5)')srstressdrop((i-1)*nl+1:i*nl)
    enddo
    close(201)

    !slip rates on subsources
    do k=1,subtot
    submvr(k)=meanvr(subposl(k),subposw(k),subsize(k))
    if(2.*subsize(k)>=l0)then    !subsources start from the hypocenter
        subrisetime(k)=aparam*l0/submvr(k)
    else                         !subsources start from a random point
        dumphi=ran2(idum1)*2.*pi;dumr=sqrt(ran2(idum1))*subsize(k)
        duml=dumr*cos(dumphi);dumw=dumr*sin(dumphi)
        subnucll(k)=duml+subposl(k)
        subnuclw(k)=dumw+subposw(k)
        subruptime(k)=ruptimegen(int(subnuclw(k)/w*float(nw-1))*nl+int(subnucll(k)/l*float(nl-1))+1)
        submvr(k)=submvr(k)*vrsubfact
        subrisetime(k)=aparam*2.*subsize(k)/submvr(k)
    endif
    enddo
    write(*,*)'... done.'

    open(201,file='subsources.dat')
    do k=1,subtot
    write(201,'(10e13.5)')subposl(k),subposw(k),subsize(k),submoment(k),subruptime(k),subrisetime(subtot),submvr(subtot)
    enddo
    close(201)

    !evaluating slip rates
    write(*,*)'preparing and saving slip rates...'
    allocate(sr(nt),stf(nt),outsrmoment(nt))
    open(201,file='sr.dat')
    ![modify]
    open(261,file='momentrate.dat')
    totmoment=0.
    stf=0.
    do i=1,nsr
        sr=0.
        outsrmoment = 0.0
        ruptimesr=ruptimegen(i)    !comment to go back to the version without rupt. vel. perturbations
!$omp parallel do private(j,k,time,ruptime,dum) default(shared)
        do j=1,nt
            time=dt*(j-1)
            do k=1,subtot
                dum=subsize(k)**2-(srl(i)-subposl(k))**2-(srw(i)-subposw(k))**2
                if(dum>0.)then
                    if(2.*subsize(k)>=l0)then    !subsources start from the hypocenter
                        ruptime=ruptimesr
                    else
                        ruptime=subruptime(k)+sqrt((srl(i)-subnucll(k))**2 &
                            +(srw(i)-subnuclw(k))**2)/submvr(k)
                        !ruptime=ruptimesr    !warning, uncomment if you want all subsources to start from the hypocenter
                    endif
                    if(time>ruptime.and.time<ruptime+subrisetime(k)*5.)then
                        sr(j)=sr(j)+sqrt(dum)*momentcorr*(time-ruptime)*&
                            exp(-(time-ruptime)/subrisetime(k)*pi)/(subrisetime(k)/pi)**2
                    endif
                endif
            enddo
        enddo
!$omp end parallel do
        stf(:)=stf(:)+sr(:)*srelem(i)*srmu(i)
        totmoment=totmoment+sum(sr(:))*srelem(i)*srmu(i)*dt
        outsrmoment = sr(:)*srelem(i)*srmu(i)
    
        do j=1,nt
            write(201,*) dt*(j-1),sr(j)
            write(261,*) dt*(j-1),outsrmoment(j)
        enddo
        write(201,*)
        write(201,*)
        write(261,*)
        write(261,*)
    enddo
    close(201)
    close(261)
    write(*,*)totmoment
    open(201,file='stf.dat')
    do i=1,nt
    write(201,*)dt*(i-1),stf(i)
    enddo
    write(*,*)'... done.'
    call mpi_finalize (ierr)

end program


subroutine fillpdf(pdfnl,pdfnw,lf,wf,l,w,sml,smw,pdf2d,cpdf2d,pdfoption,filename,filenl,filenw,pdfgaussl,pdfgaussw,pdfgausss)  ! creates pdf and cumulative pdf for subsource distribution
    implicit none
    integer pdfnl,pdfnw,pdfoption
    real pdf2d(pdfnl,pdfnw),cpdf2d(pdfnl,pdfnw)
    real lf,wf,l,w,sml,smw,pdfgaussl,pdfgaussw,pdfgausss
    character*256 filename
    integer i,j,k,filenl,filenw
    real cumul,pdfdl,pdfdw,slipdl,slipdw
    real,allocatable:: slip(:,:)
    integer pifrom,pito,pjfrom,pjto
    open(229,file='strongmotionarea.dat')
    write(229,*)sml,smw;write(229,*)sml+l,smw;write(229,*)sml+l,smw+w;write(229,*)sml,smw+w;write(229,*)sml,smw
    close(229)
    pdfdl=lf/real(pdfnl)
    pdfdw=wf/real(pdfnw)
    pifrom=int(sml/pdfdl)+1
    pito=int((sml+l)/pdfdl+0.999)
    pjfrom=int(smw/pdfdw)+1
    pjto=int((smw+w)/pdfdw+0.999)
    pdf2d(:,:)=0.
    select case(pdfoption)
    case(1)
        write(*,*)'uniform spatial pdf for subsources'
        pdf2d(pifrom:pito,pjfrom:pjto)=1.
    case(2)
        write(*,*)'gaussian pdf for subsources'
        do j=pjfrom,pjto
        do i=pifrom,pito
        pdf2d(i,j)=exp(-.5*((((real(i)-.5)*pdfdl-pdfgaussl)**2+((real(j)-.5)*pdfdl-pdfgaussw)**2)/pdfgausss**2)**2)
        enddo
        enddo
    case(3)
        write(*,*)'reading spatial pdf from',trim(filename)
        allocate(slip(filenl,filenw))
        write(*,*)'filenl',filenl
        write(*,*)'filenw',filenw
        slipdl=lf/real(filenl)
        slipdw=wf/real(filenw)
        open(329,file=trim(filename))
        do j=1,filenw
        read(329,*)(slip(i,j),i=1,filenl)
        enddo
        close(329)
        do j=pjfrom,pjto
        do i=pifrom,pito
        pdf2d(i,j)=slip(int(pdfdl/slipdl*(float(i)-0.5))+1,int(pdfdw/slipdw*(float(j)-0.5))+1)
        enddo
        enddo
        deallocate(slip)
    case default
        write(*,*)'wrong pdfoption!'
        stop
    end select
    !normalize and calculate cumulative distribution
    k=0
    cumul=0
    do j=1,pdfnw
    do i=1,pdfnl
    k=k+1
    cumul=cumul+pdf2d(i,j)
    cpdf2d(i,j)=cumul
    enddo
    enddo
    pdf2d=pdf2d/cumul
    cpdf2d=cpdf2d/cumul
end subroutine


function meanvr(x,y,r)   !calculate mean rupture velocity (just slowness mean over subsource depth extent)
    use ruptveloc
    implicit none
    real meanvr,x,y,r
    integer itop,ibottom,i
    itop=1
    ibottom=1
    do i=1,nlayers
    if(hvr(i)<y+r) itop=i+1
    if(hvr(i)<y-r) ibottom=i+1
    enddo
    if(itop==ibottom)then    ! stf point in the layer with the hypocenter
        meanvr=vr(itop)
    else
        meanvr=(y+r-hvr(itop-1))/vr(itop)+(hvr(ibottom)-y+r)/vr(ibottom)
        do i=ibottom+1,itop-1
        meanvr=meanvr+(hvr(i)-hvr(i-1))/vr(i)
        enddo
        meanvr=2.*r/meanvr
    endif
    end


    ! using crustal.dat
    subroutine read_crustaldat()
        use crustaldat
        implicit none
        integer i
        if(allocated(depth))return
        open(10,file='crustal.dat',action='read',status='old')
        write(*,*)'  (using mu values from file crustal.dat)'
        read(10,*)
        read(10,*)
        read(10,*)ndepth
        allocate(depth(ndepth),vp(ndepth),vs(ndepth),rho(ndepth))
        read(10,*)
        read(10,*)
        do i=1,ndepth
        read(10,*)depth(i),vp(i),vs(i),rho(i)
        enddo
        close(10)
        end


        subroutine randomruptvel(l,w,nx,ny,epicx,epicy,idum)
            ! random rupture velocity k^-2 generator for
            ! ruiz integral k-squared (rik) source model
            ! coded by frantisek gallovic (2014)
            !-------------------------------------------------------------------------------
            ! coordinate system on the fault:
            ! the origin is in the left top corner of the fault while the strike direction is to the right,
            ! the x axes is positive to the strike direction and
            ! the y axes is positive in the down-dip direction.
            use ruptveloc
            implicit none
            real*8,parameter:: pi=3.1415926535,vars=0.25d0  !maximum variations in percents
            real*8,parameter:: upperk=1.d0
            complex*16,allocatable:: speq1(:),ac(:,:),speqd(:)
            real*8,allocatable:: a(:,:),aa(:,:),c(:,:),d(:,:),d1(:,:)
            integer i,j,k,nxx,nyy,m,n,fnn,fnm,nx,ny
            real l,w,epicx,epicy
            real*8 dkx,dky,kx,ky
            real*8 dxout,dyout,dd
            real*8 cx,cy,ms,dum,nyql,nyqw,nlw,krad,corner,kcx,kcy
            real ran2
            integer ndepth1,idum,pdfoption
            integer*4 :: nnx, nnz, iostat, ipoint
            integer*4, external :: time_2d
            real*4 xg, zg, eps_init
            real*4, allocatable:: hs(:), t_rupt(:)
            real*8, allocatable:: xintrpl(:),tintrpl(:),yintrpl(:),xspln(:),yspln(:)

            write(*,*)'preparing rupture time variations with rupture velocity sigma ',vars

            !schvalne dd, protoze se diskretizuje na ridsi siti a pak se interpoluje
            m=2**int(log(dble(nx))/log(2.d0)+1.d0)
            n=2**int(log(dble(ny))/log(2.d0)+1.d0)
            dxout=l/dble(nx)
            dyout=w/dble(ny)
            dd=min(l/dble(nx),w/dble(ny))
            11  if(dd*m<l)then
                m=m*2
                goto 11
            endif
            12  if(dd*n<w)then
                n=n*2
                goto 12
            endif
            fnn=n/2+1;fnm=m/2+1

            allocate(speq1(n),ac(m/2,n),speqd(n),a(m,n),aa(m,n),c(m,n),d(nx,ny),d1(m,ny))
            allocate(hs(m*n),t_rupt(m*n))
            allocate(xintrpl(m),xspln(m),yintrpl(n),tintrpl(n),yspln(n))

            dkx=1./dd/real(m);dky=1./dd/real(n)
            kcx=upperk/l        !corner wave-number for along-strike direction
            kcy=upperk/w        !corner wave-number for down-dip direction
            nyql=dkx;nyqw=dky
            nlw=sqrt(nyql**2+nyqw**2)

            !preparing white noise spectrum

            do i=1,m
            do j=1,n
            a(i,j)=ran2(idum)-.5d0
            enddo
            enddo
            call rlft3(a,speq1,m,n,1,1)
            !      speq1=exp(cmplx(0.,atan2(imag(speq1),real(speq1))))
            do i=1,m/2
            ac(i,:)=cmplx(a(2*i-1,:),a(2*i,:))
            enddo

            !adding k^-2 by using the white noise spectrum

            do j=1,n
            if(j<=n/2+1)then
                ky=dky*real(j-1)
            else
                ky=-dky*real(n-j+1)
            endif
            do i=1,m/2+1
            kx=dkx*real(i-1)
            krad=sqrt(kx**2+ky**2)
            if(i<m/2+1.)then
                if(krad>=nlw)then
                    ac(i,j)=ac(i,j)/sqrt(1.+((kx/kcx)**2+(ky/kcy)**2)**1)
                endif
            elseif(krad>=nlw)then
                speq1(j)=speq1(j)/sqrt(1.+((kx/kcx)**2+(ky/kcy)**2)**1)
            endif
            enddo
            enddo

            !back fourier transform

            call rlft3(ac,speq1,m,n,1,-1)
            do i=1,m/2
            a(2*i-1,:)=real(ac(i,:))/m/n*2.
            a(2*i,:)=imag(ac(i,:))/m/n*2.
            enddo

            !adding the mean rupture velocity
            dum=sqrt(sum(a(:,:)**2)/m/n)
            a=a/dum*vars
            do j=1,n
            dum=(dd*(j-1)+dd/2.)
            if(dum>hvr(nlayers))then
                a(:,j)=vr(nlayers)*(1.+a(:,j))
            else
                do k=1,nlayers
                if(dum<hvr(k))then
                    a(:,j)=vr(k)*(1.+a(:,j))
                    exit
                endif
                enddo
            endif
            enddo

            !writing 2d rupture velocity distribution
            open(102,file='ruptvelgen.txt')
            do j=1,int(w/dd)
            do i=1,int(l/dd)
            write(102,'(3e14.6)') real(i-1)*dd+.5*dd,real(j-1)*dd+.5*dd,a(i,j)
            enddo
            write(102,*)
            enddo
            close(102)

            !translating rupture velocities to rupture times
            do j = 1,n
            do i = 1,m
            ipoint = i + (j-1) * m
            hs(ipoint) = dd/a(i,j)
            enddo
            enddo
            eps_init = 0.
            nnx=m
            nnz=n
            xg=epicx/dd+1.
            zg=epicy/dd+1.

            iostat = time_2d(hs, t_rupt, nnx, nnz, xg, zg, eps_init, 0)

            !preinterpolovani do vystupni diskretizace:
            do j=1,n
            yintrpl(j)=dd*(j-1)+dd/2.d0
            enddo
            do i=1,m
            do j=1,n
            ipoint = i + (j-1) * m
            tintrpl(j)=t_rupt(ipoint)
            a(i,j)=t_rupt(ipoint)
            enddo
            call spline(yintrpl,tintrpl,n,1.d30,1.d30,yspln)
            do j=1,ny
            dum=(j-1)*dyout+dyout/2.d0
            call splint(yintrpl,tintrpl,yspln,n,dum,d1(i,j))
            enddo
            enddo
            do i=1,m
            xintrpl(i)=dd*(i-1)+dd/2.d0  !schvalne dy, protoze se diskretizuje na ridsi siti a pak se interpoluje
            enddo
            do j=1,ny
            call spline(xintrpl,d1(:,j),m,1.d30,1.d30,xspln)
            do i=1,nx
            dum=(i-1)*dxout+dxout/2.d0
            call splint(xintrpl,d1(:,j),xspln,m,dum,d(i,j))
            enddo
            enddo
            open(101,file='ruptimegen.txt')
            do j=1,ny
            do i=1,nx
            write(101,'(3e14.6)') real(i-1)*dxout+.5*dxout,real(j-1)*dyout+.5*dyout,d(i,j)
            enddo
            write(101,*)
            enddo
            close(101)

            !forward fourier transform

            call rlft3(a,speq1,m,n,1,1)
            do i=1,m/2
            ac(i,:)=cmplx(a(2*i-1,:),a(2*i,:))
            enddo
            ac=ac/real(m*n/2);speq1=speq1/real(m*n/2)

            !writing amplitude spectrum along y:
            !    open(106,file='specy.txt')
            !    do i=1,n/2+1
            !      write(106,*)(i-1)*dky,abs(ac(1,i))
            !    enddo

            !writing amplitude spectrum along x:
            !    open(105,file='specx.txt')
            !    do i=1,m/2
            !      write(105,*)(i-1)*dkx,abs(ac(i,1))
            !    enddo
            !      write(105,*)(m/2)*dkx,abs(speq1(1))

            end    


            ! numerical recipes    

            subroutine rlft3(data,speq,nn1,nn2,nn3,isign)
                implicit none
                integer isign,nn1,nn2,nn3
                complex*16 data(nn1/2,nn2,nn3),speq(nn2,nn3)
                integer i1,i2,i3,j1,j2,j3,nn(3)
                double precision theta,wi,wpi,wpr,wr,wtemp
                complex*16 c1,c2,h1,h2,w
                c1=dcmplx(0.5d0,0.0d0)
                c2=dcmplx(0.0d0,-0.5d0*isign)
                theta=6.28318530717959d0/dble(isign*nn1)
                wpr=-2.0d0*sin(0.5d0*theta)**2
                wpi=sin(theta)
                nn(1)=nn1/2
                nn(2)=nn2
                nn(3)=nn3
                if(isign.eq.1)then
                    call fourn(data,nn,3,isign)
                    do 12 i3=1,nn3
                    do 11 i2=1,nn2
                    speq(i2,i3)=data(1,i2,i3)
                    11        continue
                    12      continue
                endif
                do 15 i3=1,nn3
                j3=1
                if (i3.ne.1) j3=nn3-i3+2
                wr=1.0d0
                wi=0.0d0
                do 14 i1=1,nn1/4+1
                j1=nn1/2-i1+2
                do 13 i2=1,nn2
                j2=1
                if (i2.ne.1) j2=nn2-i2+2
                if(i1.eq.1)then
                    h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
                    h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
                    data(1,i2,i3)=h1+h2
                    speq(j2,j3)=conjg(h1-h2)
                else
                    h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
                    h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
                    data(i1,i2,i3)=h1+w*h2
                    data(j1,j2,j3)=conjg(h1-w*h2)
                endif
                13        continue
                wtemp=wr
                wr=wr*wpr-wi*wpi+wr
                wi=wi*wpr+wtemp*wpi+wi
                w=dcmplx(dble(wr),dble(wi))
                14      continue
                15    continue
                if(isign.eq.-1)then
                    call fourn(data,nn,3,isign)
                endif
                return
                end


                subroutine fourn(data,nn,ndim,isign)
                    implicit none
                    integer isign,ndim,nn(ndim)
                    double precision data(*)
                    integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
                    double precision tempi,tempr
                    double precision theta,wi,wpi,wpr,wr,wtemp
                    ntot=1
                    do 11 idim=1,ndim
                    ntot=ntot*nn(idim)
                    11    continue
                    nprev=1
                    do 18 idim=1,ndim
                    n=nn(idim)
                    nrem=ntot/(n*nprev)
                    ip1=2*nprev
                    ip2=ip1*n
                    ip3=ip2*nrem
                    i2rev=1
                    do 14 i2=1,ip2,ip1
                    if(i2.lt.i2rev)then
                        do 13 i1=i2,i2+ip1-2,2
                        do 12 i3=i1,ip3,ip2
                        i3rev=i2rev+i3-i2
                        tempr=data(i3)
                        tempi=data(i3+1)
                        data(i3)=data(i3rev)
                        data(i3+1)=data(i3rev+1)
                        data(i3rev)=tempr
                        data(i3rev+1)=tempi
                        12            continue
                        13          continue
                    endif
                    ibit=ip2/2
                    1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
                        i2rev=i2rev-ibit
                        ibit=ibit/2
                        goto 1
                    endif
                    i2rev=i2rev+ibit
                    14      continue
                    ifp1=ip1
                    2       if(ifp1.lt.ip2)then
                        ifp2=2*ifp1
                        theta=isign*6.28318530717959d0/(ifp2/ip1)
                        wpr=-2.d0*sin(0.5d0*theta)**2
                        wpi=sin(theta)
                        wr=1.d0
                        wi=0.d0
                        do 17 i3=1,ifp1,ip1
                        do 16 i1=i3,i3+ip1-2,2
                        do 15 i2=i1,ip3,ifp2
                        k1=i2
                        k2=k1+ifp1
                        tempr=dble(wr)*data(k2)-dble(wi)*data(k2+1)
                        tempi=dble(wr)*data(k2+1)+dble(wi)*data(k2)
                        data(k2)=data(k1)-tempr
                        data(k2+1)=data(k1+1)-tempi
                        data(k1)=data(k1)+tempr
                        data(k1+1)=data(k1+1)+tempi
                        15            continue
                        16          continue
                        wtemp=wr
                        wr=wr*wpr-wi*wpi+wr
                        wi=wi*wpr+wtemp*wpi+wi
                        17        continue
                        ifp1=ifp2
                        goto 2
                    endif
                    nprev=n*nprev
                    18    continue
                    return
                    end


                    subroutine spline(x,y,n,yp1,ypn,y2)
                        integer n,nmax
                        double precision yp1,ypn,x(n),y(n),y2(n)
                        parameter (nmax=50000)
                        integer i,k
                        double precision p,qn,sig,un,u(nmax)
                        if (yp1.gt..99d30) then
                            y2(1)=0.d0
                            u(1)=0.d0
                        else
                            y2(1)=-0.5d0
                            u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
                        endif
                        do 11 i=2,n-1
                        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
                        p=sig*y2(i-1)+2.d0
                        y2(i)=(sig-1.d0)/p
                        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
                        11    continue
                        if (ypn.gt..99d30) then
                            qn=0.d0
                            un=0.d0
                        else
                            qn=0.5d0
                            un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
                        endif
                        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
                        do 12 k=n-1,1,-1
                        y2(k)=y2(k)*y2(k+1)+u(k)
                        12    continue
                        return
                        end


                        subroutine splint(xa,ya,y2a,n,x,y)
                            integer n
                            double precision x,y,xa(n),y2a(n),ya(n)
                            integer k,khi,klo
                            double precision a,b,h
                            klo=1
                            khi=n
                            1     if (khi-klo.gt.1) then
                                k=(khi+klo)/2
                                if(xa(k).gt.x)then
                                    khi=k
                                else
                                    klo=k
                                endif
                                goto 1
                            endif
                            h=xa(khi)-xa(klo)
                            if (h.eq.0.d0) pause 'bad xa input in splint'
                            a=(xa(khi)-x)/h
                            b=(x-xa(klo))/h
                            y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
                            return
                            end


                            function ran2(idum)
                                integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
                                double precision am,eps,rnmx
                                real ran2
                                parameter (im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1,ia1=40014,ia2=40692,iq1=53668,&
                                    iq2=52774,ir1=12211,ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=3.d-16,rnmx=1.d0-eps)
                                integer idum2,j,k,iv(ntab),iy
                                save iv,iy,idum2
                                data idum2/123456789/, iv/ntab*0/, iy/0/
                                if (idum.le.0) then
                                    idum=max(-idum,1)
                                    idum2=idum
                                    do 11 j=ntab+8,1,-1
                                    k=idum/iq1
                                    idum=ia1*(idum-k*iq1)-k*ir1
                                    if (idum.lt.0) idum=idum+im1
                                    if (j.le.ntab) iv(j)=idum
                                    11      continue
                                    iy=iv(1)
                                endif
                                k=idum/iq1
                                idum=ia1*(idum-k*iq1)-k*ir1
                                if (idum.lt.0) idum=idum+im1
                                k=idum2/iq2
                                idum2=ia2*(idum2-k*iq2)-k*ir2
                                if (idum2.lt.0) idum2=idum2+im2
                                j=1+iy/ndiv
                                iy=iv(j)-idum2
                                iv(j)=idum
                                if(iy.lt.1)iy=iy+imm1
                                ran2=min(am*iy,rnmx)
                                return
                            end function
