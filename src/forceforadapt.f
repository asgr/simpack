      subroutine orbitfor(poso,velo,accelo,
     & masso,No,gpeo,keo,timestepo,Nstep,soft,dolist,adapt,
     & minadsteps, maxadsteps)

      implicit none

      integer No,Nstep,io,jo,ko,year, dolist(No)
      real*8 poso(No,3), velo(No,3), accelo(No,3), accelnewo(No,3)
      real*8 masso(No), gpeo(No), timestepo, soft, keo(No), totvel(No)
      real*8 totaccelo(No), adapt
      integer minadsteps, maxadsteps
C     Conversion m/s to km/s
      real*8 kms_to_ms
      parameter (kms_to_ms=1.d3) 
C     Converstion pc to m
      real*8 pc_to_m
      parameter (pc_to_m=3.08568d16)
C     Units
      real*8 munit,lunit,vunit,convert
      parameter (munit=1.e10,lunit=1.e6,vunit=1.)
      year=31557600
      convert=(year*1e6)/pc_to_m
C     This converts 977.8139 km/s = kpc/Myr, i.e. convert is 1/977.8139
      convert=convert/1e3
C     This gets us to Mpc/Myr since internally distances are all in Mpc

C     Set accelo to 0

C     Get initial acceleration
     
      do  5 io= 1, No
          
          accelo(io,1)=0
          accelo(io,2)=0
          accelo(io,3)=0

5     continue

      call accelfor(poso,masso,No,accelo,gpeo,soft,dolist)

      do 15 io= 1, No

          accelnewo(io,1)=accelo(io,1)
          accelnewo(io,2)=accelo(io,2)
          accelnewo(io,3)=accelo(io,3)

15    continue

      do 10 jo= 1, Nstep

          do  20 io= 1, No

          poso(io,1)=poso(io,1)+(velo(io,1)*timestepo
     &    +(accelo(io,1)*timestepo**2)/2)*convert
          poso(io,2)=poso(io,2)+(velo(io,2)*timestepo
     &    +(accelo(io,2)*timestepo**2)/2)*convert
          poso(io,3)=poso(io,3)+(velo(io,3)*timestepo
     &    +(accelo(io,3)*timestepo**2)/2)*convert

20        continue

C     Get at r[jo+1] acceleration

          call accelfor(poso,masso,No,accelnewo,gpeo,soft,dolist)
      
          do  30 io= 1, No	

          velo(io,1)=velo(io,1)
     &    +(accelo(io,1)+accelnewo(io,1))*timestepo/2
          velo(io,2)=velo(io,2)
     &    +(accelo(io,2)+accelnewo(io,2))*timestepo/2
          velo(io,3)=velo(io,3)
     &    +(accelo(io,3)+accelnewo(io,3))*timestepo/2

30        continue

          do  40 io= 1, No	

          accelo(io,1)=accelnewo(io,1)
          accelo(io,2)=accelnewo(io,2)
          accelo(io,3)=accelnewo(io,3)

	  if(dolist(io).eq.0) then

          totvel(io)=sqrt(velo(io,1)**2+velo(io,2)**2+velo(io,3)**2)
          totaccelo(io)=sqrt(accelo(io,1)**2+accelo(io,2)**2
     &    +accelo(io,3)**2)
	  keo(io)=0.5*masso(io)*totvel(io)**2

          dolist(io)=floor(adapt*(totvel(io)/totaccelo(io)))

          else

             dolist(io)=dolist(io)-1

          endif

          if(dolist(io).lt.minadsteps) then
             dolist(io)=0
          endif

          if(dolist(io).gt.maxadsteps) then
             dolist(io)=maxadsteps
          endif


40        continue

10    continue

      return

      end

      subroutine accelfor(posa,massa,Na,accela,gpea,soft,dolist)

      implicit none

C     Classify integers
      integer Na,ia,ja,ka,year, dolist(Na)
C     Classify reals
      real*8 r,rT2,g,a,rx,ry,rz,convert,soft,softT2,kea(Na)
C     Setup input/output parameters
      real*8 posa(Na,3),  massa(Na), accela(Na,3), gpea(Na)
C     Gravitational constant, MKS units
      real*8 grav
      parameter (grav=6.67384d-11)
C     Conversion m/s to km/s
      real*8 kms_to_ms
      parameter (kms_to_ms=1.d3) 
C     Conversion Msol to kg
      real*8 msol_to_kg
      parameter (msol_to_kg=1.98892d30) 
C     Converstion pc to m
      real*8 pc_to_m
      parameter (pc_to_m=3.08568d16) 
C     Units
      real*8 munit,lunit,vunit
      parameter (munit=1.e10,lunit=1.e6,vunit=1.)
C     Soft
      softT2=soft**2

C     Evaluate gravitational constant in desired units
      g = grav*msol_to_kg/pc_to_m/kms_to_ms/kms_to_ms
      g = g*munit/(vunit*vunit*lunit)
C     Final units are (Lunit*pc)*(Vunit*km)^2/(Munit*Msun*s^2)
      year=31557600
      convert=1e6*year*(kms_to_ms*vunit)/(lunit*pc_to_m)

      do  10 ia= 1, Na

      if(dolist(ia).eq.0) then
      
      gpea(ia)=0

          do 20 ja= 1, Na

		   if(ia.ne.ja) then

                   rx=posa(ja,1)-posa(ia,1)
                   ry=posa(ja,2)-posa(ia,2)
                   rz=posa(ja,3)-posa(ia,3)
                   rT2=(rx**2+ry**2+rz**2)+softT2
                   a=g*massa(ja)/rT2
C		   a in units ((kms_to_ms*vunit)^2/(Lunit*pc_to_m))/s
                   r=sqrt(rT2)
                   accela(ia,1)=accela(ia,1)+a*rx/r
                   accela(ia,2)=accela(ia,2)+a*ry/r
                   accela(ia,3)=accela(ia,3)+a*rz/r
                   gpea(ia)=gpea(ia)+g*massa(ja)*massa(ia)/r
C		   gpe in units (Munit*Msun)*(kms_to_ms*vunit)^2 This definitely works!
                   endif

20        continue

          accela(ia,1)=accela(ia,1)*convert
          accela(ia,2)=accela(ia,2)*convert
          accela(ia,3)=accela(ia,3)*convert
C	  accel in units (Vunit*km/s)/Myr This definitely works!

      endif

10    continue

      return

      end

      subroutine potdenfor(posa, massa, Na, soft, grida,
     & potgrida, dengrida, bina, Ncella)

      implicit none

C     Classify integers
      integer Na, ia, ja, ka, la, Ncella
C     Classify reals
      real*8 r,g,rx,ry,rz,soft,softT2, bina(3)
C     Setup input/output parameters
      real*8 posa(Na,3), massa(Na)
      real*8 grida(Ncella,3)
      real*8 potgrida(Ncella)
      real*8 dengrida(Ncella)
C     Gravitational constant, MKS units
      real*8 grav
      parameter (grav=6.67384d-11)
C     Conversion m/s to km/s
      real*8 kms_to_ms
      parameter (kms_to_ms=1.d3) 
C     Conversion Msol to kg
      real*8 msol_to_kg
      parameter (msol_to_kg=1.98892d30) 
C     Converstion pc to m
      real*8 pc_to_m
      parameter (pc_to_m=3.08568d16) 
C     Units
      real*8 munit,lunit,vunit
      parameter (munit=1.e10,lunit=1.e6,vunit=1.)
C     Soft
      softT2=soft**2

C     Evaluate gravitational constant in desired units
      g = grav*msol_to_kg/pc_to_m/kms_to_ms/kms_to_ms
      g = g*munit/(vunit*vunit*lunit)
C     Final units are (Lunit*pc)*(Vunit*km)^2/(Munit*Msun*s^2)

      do  10 ia= 1, Ncella

            do 40 la= 1, Na

                   rx=ABS(posa(la,1)-grida(ia,1))
                   ry=ABS(posa(la,2)-grida(ia,2))
                   rz=ABS(posa(la,3)-grida(ia,3))
                   r=sqrt((rx**2+ry**2+rz**2)+softT2)
                   potgrida(ia)=potgrida(ia)-g*massa(la)/r
		   if(rx .le. (bina(1)/2) .and. ry .le. (bina(2)/2)
     & 		     .and. rz .le. (bina(3)/2)) then
                     dengrida(ia)=dengrida(ia)+massa(la)
		   endif
C		   gpe in units (Munit*Msun)*(kms_to_ms*vunit)^2 This definitely works!

40          continue

10    continue

      return

      end


