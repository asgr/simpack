C     Subroutine to generate an equilibrated dark matter halo with a Hernquist 
C     profile and an isotropic velocity distribution.
C
C     Inputs 
C     npart -- number of particles
C     sigma -- velocity dispersion of halo, in km/s
C     conc -- concentration of underlying Hernquist halo
C     iseed -- random number seed
C
C     Outputs
C     sigma -- velocity dispersion of halo, in km/s
C     (x,y,z) -- positions, in Mpc
C     (vx,vy,vz) -- velocities, in km/s

      subroutine genhalofortran(npart,x,y,z,vx,vy,vz,mpart,sigma,conc,
     &     iseed)
  
      implicit none

C     Overdensity criterion      
      real*8 delta_vir
C     Critical density
      real*8 rho_crit
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
C     PI
      real*8 PI
C     Accuracy parameter
      real*8 acc
      parameter (acc=1d-6)

C     Gravitational constant in code units
      real(kind=8) :: g         
C     Characteristic density of halo      
      real(kind=8) :: rho_0     
C     Characteristic density of halo
      real(kind=8) :: pot_0     
C     Virial radius      
      real(kind=8) :: rvir      
C     Scale radius
      real(kind=8) :: rscale    

      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      
C     Tabulated quantities
      integer*4 ntab
      parameter(ntab=5000)
      integer*4 nkeep
      real*8 rtab(ntab),ptab(ntab),mtab(ntab)
      real*8 dtab(ntab),ddtab(ntab),dddtab(ntab)
      real*8 drtab(ntab),ddrtab(ntab),dptab(ntab)
      real*8 etab(ntab),ftab(ntab)
      real*8 xx
      
      common /array/ ptab,dptab,xx,nkeep
      
C     Quantities for looping over energy and radius
      real*8 lemin,lemax,dle,le
      real*8 lrmin,lrmax,dlr,lr
      
      real*8 frmax
      real*8 rmax
      
      real*8 mtot,fsum,fmax,eran,fran,fsamp,mvir,msum,mtotal,sigma

C     Abel integrand      
      real*8 fabel   
      external fabel
      
C     Virial radius, mass, sigma
      real*8 get_rvir,get_mvir,get_sigma
      external get_rvir,get_mvir,get_sigma

      integer*4 i,j,k,l,m,n,nmin,nmax
      
      real*8 conc
      
      integer*4 iseed
      real*8 ran3
      external ran3
      
      real*8 xmin,xmid,xmax
      real*8 radius,rcynd,phi,theta,mofr,vmag,vesc,mpart
C     Particle mass
      real*8 mp
C     Positions
      real*8 x(npart),y(npart),z(npart)
C     Velocities
      real*8 vx(npart),vy(npart),vz(npart)
C     Radii
      real*8 r(npart)
      real*8 rav,r0,v1,v2,v3,vrav,x1,x2,x3

      integer*4 ncount,npart,nvir
      
C     Number of particles per species
      integer*4 np(6), nall(6)
      real*8 massarr(6)
      real*8 time
      parameter(time=0.0)
      real*8 dummy
      integer*4 idummy
      parameter(idummy=0)
      integer*4 NumFiles
      parameter(NumFiles=1)
      integer*4 highword(6)
      character*1 buffer(56)
      
      integer*4 nit
      integer*4 nitmax
      parameter(nitmax=10000)

      real*8 rho_hern,pot_hern,get_mtot_hern,get_rho0_hern
      external rho_hern,pot_hern,get_mtot_hern,get_rho0_hern
      
C     Evaluate gravitational constant in desired units
      
      g = grav*msol_to_kg/pc_to_m/kms_to_ms/kms_to_ms
      g = g*munit/(vunit*vunit*lunit)

      delta_vir=200.
      rho_crit=27.755
      PI=3.14159265359

C     Now evaluate characteristic halo properties

C     if(sigma.lt.0.0) then
C         mvir=sigma
C     else
C        mvir=get_mvir(sigma)   ! Virial mass in 1e10 Msol/h
C     end if
      
      mvir=npart*mpart
      sigma=get_sigma(mvir)

      rvir=get_rvir(mvir)       ! Virial radius in Mpc/h
      rscale=rvir/conc          ! Scale radius in Mpc/h
      
      rho_0=get_rho0_hern(rvir/rscale) ! Characteristic density
      
      mtot=get_mtot_hern(1d30)

C     Now tabulate the density and potential, for inversion
      pot_0 = pot_hern(1.d10)
  
      lrmin=log10(rscale)-15.0
      lrmax=log10(rscale)+15.0
      dlr=(lrmax-lrmin)/real(ntab-1)
      lr=lrmin
      
      n=0
      
      do
         if((lr-lrmax).gt.1.0d-6) exit
         n=n+1
         rtab(n)=10**lr
         dtab(n)=rho_hern(rtab(n))/mtot
         ptab(n)=(-pot_hern(rtab(n))+pot_0)*rscale/(g*mtot)
         lr=lr+dlr
      end do
      
C Compute the 1st derivative of radius with respect to relative gravitational potential
      do n=2,ntab-1
         drtab(n)=0.5*(((dtab(n+1)-dtab(n))/(ptab(n+1)-ptab(n))+
     &        (dtab(n)-dtab(n-1))/(ptab(n)-ptab(n-1))))
      end do
      
      drtab(1)=drtab(2)
      drtab(ntab)=drtab(ntab-1)
      
C Compute the 2nd derivative of density with respect to radius
      do n=2,ntab-1
         dptab(n)=0.5*(((drtab(n+1)-drtab(n))/(ptab(n+1)-ptab(n))+
     &        (drtab(n)-drtab(n-1))/(ptab(n)-ptab(n-1))))
      end do
      
      dptab(1)=dptab(2)
      dptab(ntab)=dptab(ntab-1)
      
! Compute the distribution function
      
      lemin=-10.0
      lemax=2.0
      dle=(lemax-lemin)/real(ntab-1)
      
      le=lemin
      
      nkeep=ntab
      
      n=0
      
      do 
         if(le.gt.lemax) exit
         if(le.ge.dlog10(ptab(1))) exit
         
         xx=10**le

         call qromo(fabel,1.0d-20,xx,fsum)
         n=n+1
         etab(n)=xx
         ftab(n)=fsum/sqrt(8.)/PI/PI
         le=le+dle
      end do
      
C     Now generate the halo by sampling the mass profile and using accept/reject
C     to sample the distribution function for velocities...

      nkeep=n
      
      ncount=0; nvir=0
      
      xmin=1.d-12
      xmax=1.d12
      
      frmax=1.0

      mpart=mvir/real(npart)
      
      npart = 1+floor(frmax*mtot/mpart)
      
      mtotal=0.0
      
      do
         if(ncount.ge.npart) exit

C     Choose a random mass

         mofr=frmax*ran3(iseed)

         xmin=1.d-12
         xmax=1.d12
         
         if((get_mtot_hern(xmin)/mtot-mofr)*(get_mtot_hern(xmax)/mtot-
     &        mofr).gt.0.0) cycle
         
         do
            xmid=0.5*(xmin+xmax)
            if((get_mtot_hern(xmid)/mtot-mofr)*(get_mtot_hern(xmax)/mtot
     &           -mofr).ge.0.0) then
               xmax=xmid
            else
               xmin=xmid
            end if
            
            if(dabs(xmax-xmin).lt.acc) exit
         end do
 
         radius=rscale*xmid
         
C     We know what the radius is -- now we need to estimate the relative potential...
         i=1
         do
            if(rtab(i).gt.radius) exit
            i=i+1
         end do
         
         phi=ptab(i+1)*(radius-rtab(i))/(rtab(i+1)-rtab(i))+ptab(i)
     &        *(rtab(i+1)-radius)/(rtab(i+1)-rtab(i))
         
         vesc=dsqrt(2.d0*phi)
         
         nit=0
         
         do
            if(nit.gt.nitmax) exit
            
            nit=nit+1
            
            v1=(2*ran3(iseed)-1.)*vesc
            v2=(2*ran3(iseed)-1.)*vesc
            v3=(2*ran3(iseed)-1.)*vesc
            
            eran=phi-0.5*(v1*v1+v2*v2+v3*v3)
            
            if(eran.lt.0.0) cycle
            
            nmin=1
            nmax=nkeep
            
            if((etab(1)-eran)*(etab(nmax)-eran).gt.0.0) cycle
            
            do
               n=(nmin+nmax)/2
               if((etab(n)-eran)*(etab(nmax)-eran).ge.0.0) then
                  nmax=n
               else
                  nmin=n
               end if
               
               if(abs(nmax-nmin).le.1) exit
            end do
        
            fsamp=ftab(nmax)*(eran-etab(nmin))/(etab(nmax)-etab(nmin))+
     &           ftab(nmin)*(etab(nmax)-eran)/(etab(nmax)-etab(nmin)) 
            
C     We also need to estimate where the maximum of the distribution function lies...

            nmin=1
            nmax=nkeep
            
            if((etab(1)-phi)*(etab(nmax)-phi).gt.0.0) cycle
            
            do
               n=(nmin+nmax)/2
               if((etab(n)-phi)*(etab(nmax)-phi).ge.0.0) then
                  nmax=n
               else
                  nmin=n
               end if
               
               if(abs(nmax-nmin).le.1) exit
            end do
            
            fmax=ftab(nmax)*(phi-etab(nmin))/(etab(nmax)-etab(nmin))+
     &           ftab(nmin)*(etab(nmax)-phi)/(etab(nmax)-etab(nmin)) 
            
            fran=fmax*ran3(iseed)
            
            if(fsamp.gt.fran) exit
            
         end do
            
 2102    ncount=ncount+1    
         
C     Now assign the selected radii cartesian coordinates.
C     These will be uniformly distributed on the unit sphere.
         
         phi=2*PI*ran3(iseed)
         
         do
            x1=(2.0*ran3(iseed)-1.0)
            rcynd=sqrt(1.0-x1*x1)
            x2=rcynd*cos(phi)
            x3=rcynd*sin(phi)
            if(abs(1.0-(x1*x1+x2*x2+x3*x3)).lt.acc) exit
         end do
         
         x(ncount)=radius*x1
         y(ncount)=radius*x2
         z(ncount)=radius*x3
         
         vx(ncount)=v1*sqrt(g*mtot/rscale)
         vy(ncount)=v2*sqrt(g*mtot/rscale)
         vz(ncount)=v3*sqrt(g*mtot/rscale)
         
         mtotal=mtotal+mpart

      end do
      
      npart=ncount

      return

      end subroutine genhalofortran

C     The Abel integral, needed to compute the distribution function. Can be derived
C     analytically for some halo profiles, but for generality assume that it is computed
C     numerically and tabulated.
      
      function fabel(y)

      implicit none

      real*8 y,p,dp,fabel

      integer*4 n,nmin,nmax

      integer*4 ntab
      parameter(ntab=5000)
      integer*4 nkeep

      real*8 xx
      real*8 ptab(ntab),dptab(ntab)

      common /array/ ptab,dptab,xx,nkeep

      nmin=1
      nmax=nkeep
      
C      if((ptab(1)-y)*(ptab(nmax)-y).gt.0.0) stop 
C     &     '(ptab(1)-y)*(ptab(nmax)-y).gt.0.0'
      
      do
         n=(nmin+nmax)/2
         if((ptab(n)-y)*(ptab(nmax)-y).ge.0.0) then
            nmax=n
         else
            nmin=n
         end if
         
         if(abs(nmax-nmin).le.1) exit
      end do
      
      dp = dptab(nmin)*(ptab(nmax)-y)/(ptab(nmax)-ptab(nmin)) + 
     &     dptab(nmax)*(y-ptab(nmin))/(ptab(nmax)-ptab(nmin))
      
      fabel = dp/dsqrt(xx-y)

      return

      end 

C     Compute virial radius, given virial mass

      function get_rvir(mvir)

      implicit none
      
      real*8 mvir
      real*8 get_rvir

      real*8 PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      
      get_rvir = (3.*mvir/4/PI/rho_crit/delta_vir)**(1./3.)
      
      return
      
      end 

C     Compute virial mass, given velocity dispersion of halo

      function get_mvir(sigma)
      
      implicit none
      
      real*8 sigma
      real*8 get_mvir

      real*8 PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      
      get_mvir = sqrt(3./(32.*PI*g*g*g*delta_vir*rho_crit))*sigma*sigma
     &     *sigma
 
      return
      
      end 

C     Compute velocity dispersion of halo given the virial mass

      function get_sigma(mvir)
      
      implicit none
      
      real*8 mvir
      real*8 get_sigma

      real*8 PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      
      get_sigma = (mvir*sqrt(32.*PI*g*g*g*delta_vir*rho_crit/3.))
     &     **(1./3.)
      
      return
      
      end 

      SUBROUTINE qromo(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=2.e-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do j=1,JMAX
         call midpnt(func,a,b,s(j),j)
         if (j.ge.K) then
            call polint(h(j-KM),s(j-KM),K,dble(0.),ss,dss)
            if (abs(dss).le.EPS*abs(ss)) return
         endif
         s(j+1)=s(j)
         h(j+1)=h(j)/9.
      end do
C     pause 'too many steps in qromo'
      end 

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
         dift=abs(x-xa(i))
         if (dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
 11   continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
         do 12 i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
C           if(den.eq.0.)pause 'failure in polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
 12      continue
         if (2*ns.lt.n-m)then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
 13   continue
      return
      end 

      SUBROUTINE midpnt(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
         s=(b-a)*func(0.5*(a+b))
      else
         it=3**(n-2)
         tnm=it
         del=(b-a)/(3.*tnm)
         ddel=del+del
         x=a+0.5*del
         sum=0.
         do 11 j=1,it
            sum=sum+func(x)
            x=x+ddel
            sum=sum+func(x)
            x=x+del
 11      continue
         s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END 

      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
         iff=1
         mj=MSEED-iabs(idum)
         mj=mod(mj,MBIG)
         ma(55)=mj
         mk=1
         do i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.MZ)mk=mk+MBIG
            mj=ma(ii)
         end do
         do k=1,4
            do i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
            end do
         end do
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END

C     Hernquist profile central density
      function get_rho0_hern(c)
      
      implicit none
      
      real*8 c
      real*8 get_rho0_hern

      real*8 PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      
      get_rho0_hern = delta_vir * (2./3.) * c * (1+c)*(1+c) * rho_crit
      
      return
      
      end
      
C     Hernquist profile total mass
      function get_mtot_hern(c)
      
      implicit none
      
      real*8 c
      real*8 get_mtot_hern

      real*8 PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      
      get_mtot_hern=2.*PI*rho_0*rscale**3 * (c*c)/(1+c)**2
      
      return
      
      end

C     Hernquist profile total mass      
      function rho_hern(r)
      
      implicit none
      
      real*8 r,x
      real*8 rho_hern

      real*8 PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      
      x=r/rscale
      
      rho_hern = rho_0 * (x*(1+x)*(1+x)*(1+x))**(-1)
      
      return
      
      end
      
C     Hernquist profile total mass
      function pot_hern(r)
      
      implicit none
      
      real*8 r,x
      real*8 pot_hern

      real*8 PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale
      common /params/ PI,delta_vir,rho_crit,g,rho_0,pot_0,rscale

      x=r/rscale
      
      pot_hern = -4.*PI*g*rho_0*rscale*rscale*(2.*(1+x))**(-1.)
      
      return
      
      end 
