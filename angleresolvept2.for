c	This program computes depth profiles from angle resolved XPS data.
c	sig(i,j)   = observed signal at theta(i) for species j
c			NOTE: j=1 must be underlying substrate species
c	sigma(i,j) = standard deviation in sig(i,j) (assumed to be 0.1*sig(i,j))
c	theta(i)   = take-off angle (wrt surface)
c	sth(i)	   = relative detection efficiency of electrons emerging
c		     from the surface as a function of theta DEGREES (NOT
c		     RADIANS), determined for this SSL system from C 
c		     surface impurities on Si.
c	l(k,j)	   = inelastic mean free path for slab k
c	c(k,j)     = concentration of species A(j) at depth k*dz for
c		     slab i
c	in(i,j)	   = intensity at theta(i) of species j summed over all slabs
c	dz	   = slab height
c	eal	   = electron attenuation length in angstroms for the
c		     underlying substrate at the electron kinetic energy
c		     of the substrate material (22 for Si)
c
c
c	input file: concvsth.dat
c		file structure:
c			nspec     iatom(1)     iatom(2)     ...
c			be(1)     be(2)        ...
c			theta(1)  sig(1,1)     sig(1,2)     ...
c			theta(2)  sig(2,1)     sig(2,2)     ...
c			...
      program angleresolvept   

c
      double precision a,b1,b2,b3,b4,b5,b6,b7
      double precision c,theta,sth,cold,ent,entold
      real*8 in(20,10),insub(20),insun(20),inold(20,10)
      real na,msi,nsi,nsio2
      real l(200)
      real deg(20)
c removed l(200) it was unused
      common /blk1/ in,sig(20,10),sigma(20,10),insub,inold,insun
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg
      common /blk3/ iatom(10)
      common /blk4/ atten(120,200), be(10), xray, senexp
      common /blk5/ c0(200,10), cold(10)
      dimension d0(20), alp(20), sumc(10), chis10(10)
c
      data pi,ealf,eals,xray,niter,densi,msi,na
     * /3.14159,35.0,22.0,1486.7,100,2.35,28.0,6.023e23/
      data nsi,nsio2 /5.00e22,5.00e22/
c	nsio2=6.92e22
      eal=ealf
      write(*,*)'enter npt and alpha0'
      read(*,*)npt, alpha0
c
c	Set input = 0 for input from keyboard, 1 for SSL values or 2
c	for Lehigh values with no keyboard input
c
      input=0
c
      if(input.eq.1)then
      imach=1
      icolang=15
      ires=4
      iforce=0
      endif
c
      if(input.eq.2)then
      imach=2
      iforce=0
      go to 99
      endif
c
      if(input.eq.0)then
      write(*,*)'Enter 1 for SSL, 2 for Lehigh, or 3 for Phi machine'
      read(*,*)imach
      if(imach.eq.1)then
      icolang=15
      write(*,*)'Enter 1 for high resolution, 4 for low resolution'
      read(*,*)ires
      endif
      if(imach.eq.2)then
      write(*,*)'sth not fixed for >15 and < 15 deg yet'
      go to 99
      write(*,*)'What is +- collection angle?'
      read(*,*)icolang
      endif
      if(imach.eq.3)then
      icolang=7
      endif
c	   write(*,*)'Enter 1 to force species 1 to zero at surface, any'
c	   write(*,*)'number to let it float'
c	   read(*,*)iforce
      iforce=0
      endif
c	ilast=int(4.0*eal)
      ilast=200
      if(ilast.gt.200)then
      write(*,*)'ilast too large for dimension allocations'
      go to 99
      endif
      dz=1.0
c	alpha0=0.4
c	alpha0=0.1
      senexp=0.7
c
      open(10,file='concvsth.dat')
      rewind 10
      open(11,file='concvsd.dat')
      rewind 11
      open(12,file='sigvsth.dat')
      rewind 12
      open(13,file='signorm.dat')
      rewind 13
      open(14,file='converge.dat')
      rewind 14
      open(15,file='areal.dat')
      rewind 15
c

      read(10,*)nspec,(iatom(j),j=1,nspec)
      read(10,*)(be(j),j=1,nspec)
      do i=1,20
      
c	   deg(i) is read in degrees and then converted into theta(i)
c	   in radians

      read(10,*,end=10)deg(i),(sig(i,j),j=1,nspec)
c
c	Set the sigmas of the intensities to 3% of the intensity
c
      do j=1,nspec
      sigma(i,j)=0.03*sig(i,j)
      enddo
c
c	   NOTE: First compute sth(i) from theta(i) in degrees, and THEN
c	   convert theta into radians
c
      if(deg(i) .gt. 90 .and. deg(i) .le. 360)deg(i)=180.0-deg(i)
c errors maybe???
      if(imach .eq. 1)then
c another error???
c		SSL machine
c
c		Low res (#4):
c
      if(ires.eq.4)then
      A=0.38421
      B1=0.00472
      B2=2.33479E-5
      B3=0.0
      B4=0.0
      B5=0.0
      B6=0.0
      B7=0.0
      endif
c
c		From high res (#1) scans
c
      if(ires.eq.1)then
      A=0.71814
      B1=0.00333
      B2=-2.0294E-6
      B3=0.0
      B4=0.0
      B5=0.0
      B6=0.0
      B7=0.0
      endif
      endif
      if(imach.eq.2)then
      if(theta(i).lt.15.0)then
c
c		   Lehigh machine (toa=0-15)
c		   ----------------------------------------------------
      A=0.6763
      B1=0.09341036
      B2=0.002702
      B3=0.0
      B4=0.0
      B5=0.0
      B6=0.0
      B7=0.0
c		   ----------------------------------------------------
c
      endif
      if(theta(i).ge.15.0)then
c
c		   Lehigh machine (toa=15-90)
c		   ---------------------------------------------------
      A=5.95091
      B1=-0.33508
      B2=0.01011
      B3=-1.63926E-4
      B4=1.45911E-6
      B5=-6.69174E-9
      B6=1.23921E-11
      B7=-4.88415E-23
c		   ---------------------------------------------------
      endif
      endif
      if(imach.eq.3)then
      A=0.4368
      B1=0.01252
      B2=-6.95401E-5
      B3=-6.81597E-20
      endif
      if(imach .gt. 3 .or. imach .lt. 1)go to 99
      theta(i)=deg(i)*pi/180.0
      nang=i
      enddo
c
10    do i=1,90+icolang
c	   NOTE: Compute sth(i) from angle in degrees
c
      if(i.le.90)degi=float(i)
      if(i.gt.90)degi=float(i-2*(i-90))
c	degi=deg(i)
      sth(i)=A+B1*degi+B2*degi**2+B3*degi**3
     * +B4*degi**4+B5*degi**5+B6*degi**6+B7*degi**7
      write(*,*)degi,sth(i)
      enddo
c
c	Normalize signals and sigmas at each angle theta(i) to intensity
c	and of the bulk component (j=1) at theta(i), i.e. sig(i,1)=1.0
c	for all i.  This throws away some information, but removes any error
c	due to the sample not being at the right height at any angle, or to
c	shadowing of some of the signal, due to the clips on the sample
c	holder, etc.  Maybe this should be an option in future work.
c
      do i=1,nang
      do j=2,nspec
      sig(i,j)=sig(i,j)/sig(i,1)
      sigma(i,j)=sigma(i,j)/sig(i,1)
      enddo
      sigma(i,1)=sigma(i,1)/sig(i,1)
      sig(i,1)=1.0
      write(13,11)theta(i)*180.0/pi,(sig(i,j),j=1,nspec)
      enddo
c
c
c	First assume a uniform layer with no Si (j=1).  Compute its
c	thickness for the different take-off angles.
c
c	Compute the thickness of the layer assuming that no substrate
c	material is in the layer.
c
      do i=1,nang
      sumsig=0.0
      do k=2,nspec
      sumsig=sumsig+iatom(k)*(sig(i,k)/sig(i,1))*
     * ((xray-be(1))/(xray-be(j)))**senexp
      enddo
      d0(i)=ealf* sin(theta(i))*
     * log(1.0+(eals/ealf)*sumsig*(nsi/nsio2))
      write(*,*)theta(i)*180.0/pi,d0(i)
      enddo
c
c	Calculate the proportionality constant a
c
      do i=1,nang
      alp(i)=0.0
      do j=2,nspec
      alp(i)=alp(i)+iatom(j)*sig(i,j)
      enddo
      alp(i)=1.0/alp(i)
      enddo
c
c	Calculate the number of total number of atoms in all species
c
      sumat=0.0
      do j=1,nspec
      sumat=sumat+iatom(j)
      enddo
c
c	Set up initial depth profile, using the values at the last (assumed
c	to be largest) theta(i)
c
      do k=1,ilast
      z(k)=k*dz
      do j=1,nspec
      if(z(k).le.d0(nang))then
      if(j.eq.1)c(k,j)=0.0
      if(j.ne.1)c(k,j)=alp(nang)*sig(nang,j)*nsio2/nsi
      else
      if(j.eq.1)c(k,j)=1.0
      if(j.ne.1)c(k,j)=0.0
      endif
      enddo
      enddo
c   Do 5-point smooth of the initial profile one time
c	npt=5
c	npt=3
      ntimes=1
      do j=1,nspec
      call smooth(npt,ntimes,j,ilast)
      enddo
      do k=1,ilast
      call normal(k,nspec,nsi,nsio2)
      enddo
      call beers(nang,nspec,ilast,eals,ealf,icolang)
c
c	Store this as the initial estimate of the depth profile
c
      do k=1,ilast
      do j=1,nspec
      c0(k,j)=c(k,j)
      if(c0(k,j).lt.1.0e-6)c0(k,j)=1.0e-6
      enddo
      enddo

c	Compute signals at all theta's from this depth profile
c
c
      call exsignal(nang,nspec,ilast,kjunk,icolang)
c
      call chisqr(chi,nang,nspec)
      call entrop(ent,ilast,nspec,kjunk,entjunk)
      chilst=chi-100.0*alpha0*ent
      entold=ent
c
c	MAIN ROUTINE BEGINS HERE:
c
      do iter=1,niter
      if(iter.le.100)then
c		alpha=(100-iter)*alpha0
      alpha=20.0*alpha0
      else
      alpha=20.0*alpha0
      endif
      if(iter.le.50)then
      delcon=0.01
c		delcon=0.02
      else
      delcon=0.01
c		delcon=0.02
      endif
      do k=1,ilast
      do j=1,nspec
      do j2=1,nspec
      cold(j2)=c(k,j2)
      enddo
      c(k,j)=c(k,j)+delcon
      do j2=1,nspec
      if(j2.ne.j)c(k,j2)=c(k,j2)
     * -delcon*iatom(j)/((nspec-1)*iatom(j2))
      if(c(k,j2).lt.0.0)c(k,j2)=0.0
      enddo
      call normal(k,nspec,nsi,nsio2)
      call sigslb(nang,nspec,ilast,k,icolang)
c	call signal(nang,nspec,ilast,kjunk,icolang)
      call chisqr(chi,nang,nspec)
      call entslb(ent,ilast,nspec,k,entold)
      if(chi-alpha*ent.lt.chilst)then
      chilst=chi-alpha*ent
      entold=ent
      do i2=1,nang
      insub(i2)=insun(i2)
      do j2=1,nspec
      inold(i2,j2)=in(i2,j2)
      enddo
      enddo
      else
      do j2=1,nspec
      c(k,j2)=cold(j2)
      enddo
      if(c(k,j)-delcon.ge.0.0)then
      c(k,j)=c(k,j)-delcon
      do j2=1,nspec
      if(j2.ne.j)c(k,j2)=c(k,j2)+delcon
     * *iatom(j)/((nspec-1)*iatom(j2))
      enddo
      call normal(k,nspec,nsi,nsio2)
      call sigslb(nang,nspec,ilast,k,icolang)
c	call signal(nang,nspec,ilast,kjunk,icolang)
      call chisqr(chi,nang,nspec)
      call entslb(ent,ilast,nspec,k,entold)
      if(chi-alpha*ent.lt.chilst)then
      chilst=chi-alpha*ent
      entold=ent
      do i2=1,nang
      insub(i2)=insun(i2)
      do j2=1,nspec
      inold(i2,j2)=in(i2,j2)
      enddo
      enddo
      else
      do j2=1,nspec
      c(k,j2)=cold(j2)
      enddo
      endif
      endif
      endif
      enddo
      enddo
c	   Subtract off zeros for all but species 1 to force profiles
c	   to go to zero at max depth
c
      c1surf=c(1,1)
      do k=1,ilast
      if(iforce.eq.1 .and. k.le.5)then
      c(k,1)=0.0
c		   c(k,1)=c(k,1)-c1surf
      if(c(k,1).lt.0.0)c(k,1)=0.0
      endif
      do j=2,nspec
      c(k,j)=c(k,j)-c(ilast,j)
      if(c(k,j).lt.0.0)c(k,j)=0.0
      enddo
      call normal(k,nspec,nsi,nsio2)
      enddo
c
c	   Do 5-pt smooth
c
c	   npt=5
c	   npt=3
      ntimes=1
      do j=1,nspec
      call smooth(npt,ntimes,j,ilast)
      enddo
      do k=1,ilast
      call normal(k,nspec,nsi,nsio2)
      enddo
      call exsignal(nang,nspec,ilast,kjunk,icolang)
      call chisqr(chi,nang,nspec)
      call entrop(ent,ilast,nspec,kjunk,entjunk)
      entold=ent
      chilst=real(chi-alpha*ent)
      write(*,*)'iter,chisqr,S= ',iter,chi,ent,sumchi/10.0
      write(14,55)iter,chi,ent,sumchi/10.0
      if(iter.ge.1)then
      do ichi=2,10
c		   chis10(ichi)=chis10(ichi-1)
      chis10(12-ichi)=chis10(11-ichi)
      enddo
      chis10(1)=chilst
      sumchi=0.0
      do itot=1,10
      sumchi=sumchi+chis10(itot)
      enddo
c	write(*,333)iter,k,(chis10(itot),itot=1,10),sumchi
      if(iter.gt.15)then
      if(1.01*chilst.gt.sumchi/10.0)then
      write(*,*)'Convergence attained'
      go to 99
      endif
      endif
      endif
      call beers(nang,nspec,ilast,eals,ealf,icolang)
      enddo
c	call normal(k,nspec,nsi,nsio2)
c
c	END MAIN ROUTINE
c
99    do i=1,nang
      write(*,11)theta(i)*180.0/pi,(in(i,j)/in(i,1),j=2,nspec)
      write(12,11)theta(i)*180.0/pi,(in(i,j)/in(i,1),j=2,nspec)
      enddo
      do k=1,ilast
      write(11,22)k,(c(k,j),j=1,nspec)
      enddo
c
c	Express the ariel density (atoms/cm^2) in terms of Si atom density
c	PROBABLY OK BUT THINK ABOUT A LITTLE MORE:
c
      do j=1,nspec
      sumc(j)=0.0
      do k=1,ilast
      sumc(j)=sumc(j)+1.0e-8*dz*c(k,j)*na*densi/msi
      enddo
      enddo
      write(*,*)'Species: ',(j,j=1,nspec)
      write(*,*)(sumc(j),j=1,nspec)
      write(15,*)(be(j),j=1,nspec)
      write(15,*)(sumc(j),j=1,nspec)
c
c	Determine film thickness (where Si reaches 0.5)
c
      thickn=-1.0
      do k=1,ilast
      if(thickn .lt. 0.0 .and. c(k,1) .gt. 0.5)then
      thickn=dz*(k-1)
      endif
      enddo
      write(15,*)thickn
      close unit=10
      close unit=11
      close unit=12
      close unit=13
      close unit=14
      close unit=15
11    format(1f6.1,8f9.3)
22    format(1i4,6f10.5)
55    format(1i6,3f12.6)
c333	format(2i4,10f5.0,1f6.0)
      stop
      end
c
c
c
c	SUBROUTINES
c
c
      subroutine chisqr(chi,nang,nspec)
      real*8 in(20,10), insub(20),insun(20),inold(20,10)
      common /blk1/ in,sig(20,10),sigma(20,10),insub,inold,insun
      chi=0.0
      do i=1,nang
      do j=1,nspec
      chi=chi+((in(i,j)-sig(i,j))/sigma(i,j))**2
      enddo
      enddo
      return
      end
c
c
      subroutine entrop(ent,ilast,nspec,kjunk,entjunk)
c
c	This routine computes the entropy term the first time from the
c	entire depth profile for each species.  It only needs to be called
c	the first time through.
c
c      real l(200)
c removed l(200) it was unused
      double precision theta,sth,c,cold,ent
      real deg(20)
      real l(200)
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg  !      common /blk2/ c(200,10),z(200),dz,sth(120),l,theta(20),deg(20)
      
      common /blk5/ c0(200,10), cold(10)
      ent=0.0
      do k=1,ilast
      do j=1,nspec
      if(c(k,j).lt.1.0e-6)then
      ent=ent-c0(k,j)
      else
      ent=ent+c(k,j)-c0(k,j)-c(k,j)*dlog10(c(k,j)/c0(k,j))
      endif
      enddo
      enddo
      return
      end
c
      subroutine entslb(ent,ilast,nspec,k,entold)
c
c	This routine computes the entropy term by using the previous entropy
c	and subtracting off the term for slab k before it was changed and
c	then adding the term for slab k for this trial change.  This is much
c	faster than computing the entire entropy again, using subroutine
c	entrop.
c
      real l(200)
      real deg(20)
c removed l(200) it was unused
      double precision theta,sth,c,cold,eoldk,enewk,ent,entold
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg                         !common /blk2/ c(200,10),z(200),dz,sth(120),l,theta(20),deg(20)
      common /blk5/ c0(200,10), cold(10)
      eoldk=0.0
      enewk=0.0
      do j=1,nspec
c
c	   First compute entropy term for old values for slab k
c
      if(cold(j).lt.1.0e-6)then
      eoldk=eoldk-c0(k,j)
      else
      eoldk=eoldk+cold(j)-c0(k,j)-cold(j)*dlog10(cold(j)/c0(k,j))
      endif
c
c	   Now compute entropy term for trial new values for slab k
c
      if(c(k,j).lt.1.0e-6)then
      enewk=enewk-c0(k,j)
      else
      enewk=enewk+c(k,j)-c0(k,j)-c(k,j)*dlog10(c(k,j)/c0(k,j))
      endif
      enddo
      ent=entold-eoldk+enewk
      return
      end
c
      subroutine exsignal(nang,nspec,ilast,kjunk,icolang)
      real*8 in(20,10),insub(20),insun(20),inold(20,10)
      double precision theta,sth,c
      real in0(200,10)
      real l(200) 
      real deg(20)
      common /blk1/ in,sig(20,10),sigma(20,10),insub,inold,insun
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg                          !common /blk2/ c(200,10),z(200),dz,sth(120),l,theta(20),deg(20)
      common /blk4/ atten(120,200), be(10), xray, senexp
      do i=1,90+icolang
      do j=1,nspec
      in0(i,j)=0.0
      enddo
      enddo
      do i=1,90+icolang
      do k=1,ilast
      do j=1,nspec
      in0(i,j)=in0(i,j)+dz*sth(i)*c(k,j)*atten(i,k)**
     * ((xray-be(1))/(xray-be(j)))**senexp
c     1		   sqrt((xray-be(1))/(xray-be(j)))
      enddo
      enddo
      enddo
c
      do i=1,nang
      degsp=deg(i)
      itheta=int(degsp)
      do j=1,nspec
      in(i,j)=0.0
      do i2=itheta-icolang,itheta+icolang
c
c		   NOTE: 2 lines below were incorrectly used in some
c		   cases.  The results were unpredictable, since in()
c		   is only dimensioned for 20.  Such a treatment for
c		   angles > 90 is unnecessary anyway.
c
c		   if(i2.gt.90)in(i2,j)=in(i2-2*(i2-90),j)
c		   if(i2.ge.1.and.i2.le.90)in(i,j)=in(i,j)+in0(i2,j)/
c     1		   (2*icolang+1)
      in(i,j)=in(i,j)+in0(i2,j)/(2*icolang+1)
      enddo
      enddo
      insub(i)=in(i,1)
      do j=1,nspec
      in(i,j)=in(i,j)/insub(i)
      inold(i,j)=in(i,j)
      enddo
      enddo
      return
      end
c
c
      subroutine sigslb(nang,nspec,ilast,k,icolang)
      real*8 in(20,10),insub(20),insun(20),inold(20,10),indif0(200,10),
     * indif
      double precision theta,sth,c,cold
      real l(200)
      real deg(20)
c removed l(200) it was unused
      common /blk1/ in,sig(20,10),sigma(20,10),insub,inold,insun
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg
      common /blk4/ atten(120,200), be(10), xray, senexp
      common /blk5/ c0(200,10), cold(10)
      do i=1,90+icolang
      do j=1,nspec
c
c		Compute intensities from old and new slabs and new intensity
c		integrated over slabs
c
c		Assume electron attenuation length is inversely proportional
c		to the electron kinetic energy to the senexp power (e.g. 0.7)
c
c
      indif0(i,j)=dz*sth(i)*(c(k,j)-cold(j))*atten(i,k)**
     * ((xray-be(1))/(xray-be(j)))**senexp
c     1		sqrt((xray-be(1))/(xray-be(j)))
      enddo
      enddo
      do i=1,nang
      degsp=deg(i)
      itheta=int(degsp)
      do j=1,nspec
c		in(i,j)=0.0
      indif=0.0
      do i2=itheta-icolang,itheta+icolang
c		   SEE NOTE IN SUBROUTINE SIGNAL
c		   if(i2.gt.90)in(i2,j)=in(i2-2*(i2-90),j)
c		   if(i2.ge.1.and.i2.le.90)indif=indif+
c     1		   indif0(i2,j)/(2*icolang+1)
      indif=indif+indif0(i2,j)/(2*icolang+1)
      enddo
      in(i,j)=insub(i)*inold(i,j)+indif
      enddo
      insun(i)=in(i,1)
      do j=1,nspec
      in(i,j)=in(i,j)/insun(i)
      enddo
      enddo
      return
      end
c
      subroutine smooth(npt,ntimes,j,ilast)
      double precision theta,sth,c
      dimension cruff(-20:200)
      real l(200)
      real deg(20)
c removed l(200) it was unused
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg              !common /blk2/ c(200,10),z(200),dz,sth(120),l,theta(20),deg(20)
      ist=-((npt-3)/2)
      do i=1,ntimes
      do k=ist,ilast+(npt-1)/2
      if(k .ge. 1 .and. k .le. ilast)then
      cruff(k)= (c(k,j))
      elseif(k .lt. 1)then
      cruff(k)= (c(-k+2,j))
      elseif(k .gt. ilast)then
      cruff(k)= (c(k-2*(k-ilast),j))
      endif
      enddo
      do k=1,ilast
      c(k,j)=0.0
      do i2=1,npt
      c(k,j)=c(k,j)+cruff(k-(npt+1)/2+i2)
      enddo
      c(k,j)=c(k,j)/npt
      enddo
      enddo
      return
      end
c
      subroutine normal(k,nspec,nsi,nsio2)
      double precision theta,sth,c
      real nsi, nsio2
      real l(200) 
      real deg(20)
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg
      common /blk3/ iatom(10)
      sum=0.0
      do j=1,nspec
      if(c(k,j).lt.0.0)c(k,j)=0.0
      sum=real(sum+iatom(j)*c(k,j))
c WHAT is sum???
      enddo
      fsi=real(c(k,1)/sum)
      fsio2=1-fsi
      atomden=fsi+fsio2*nsio2/nsi
      do j=1,nspec
      c(k,j)=atomden*(c(k,j)/sum)
      enddo
      return
      end
c
      subroutine beers(nang,nspec,ilast,eals,ealf,icolang)
c
c
c	Beers law is used to compute atten(i,k) the attenuation factor
c	at angle theta(i) and distance z(k) below the surface for Si.
c
c	In principle, this routine should be called every time signals are to
c	be computed.  However, since the mean-free-path is only a mild
!	function of composition, it will be called only once per cycle
c	through the main loop.  Subroutine normal must be called before beers
c	to get the revised c(k,j) values,so that l(k)'s can be calculated.
c
      double precision theta,sth,c
      real l(200)
      real deg(20)
c removed l(200) it was unused
      common /blk2/ c(200,10), z(200), dz, sth(120), l, theta(20), deg                          !common /blk2/ c(200,10), z(200),dz,sth(120), l, theta(20),deg(20)

      common /blk4/ atten(120,200), be(10), xray, senexp
      data pi /3.14159/
c
c	Compute the electron attenuation length at depth z, assuming that
c	it is a linear combination of the fraction of film material and
c	substrate material.  Energy is for Si(2p) excited by Al K-alpha,
c	i.e. 1486.7-99.4=1387.4
c
      
      do k=1,ilast
            l(k) = c(k,1)*eals + (1.0 - c(k,1))*ealf
      end do
c
      do i=1,90+icolang
      do k=1,ilast
      rad = i * pi / 180.0
!      real rad = i * pi / 180.0
      if(k.eq.1) then
            atten(i,1) = exp(-0.5 * dz / (l(1) * sin(rad)))
      else
            atten(i,k) = atten(i,k-1) * exp(-0.5 * dz / (l(k-1))) ! * sin(rad)))
            
      ! Alternatively, if you meant the commented line below, you can uncomment it
      ! atten(i,k) = atten(i,k-1) * exp(-0.5 * dz / (l(k) * sin(rad)))
      endif
      end do
      end do


      return
      end

