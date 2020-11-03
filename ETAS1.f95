!   simulatinon of ETAS model
!     power-law decrease of aftershock density with distance p(r)~ 1/(d(m)+r)^{1+mu} in 3D
!   Characteristic triggering distance d(m)=d0*10^(gamma M) 
!   Uniform background rate (in time and space)
!   Corrective term rho^* to account for unobserved events with M0<m<Md
!     Written by Karen Felzer and Yu Gu, spring 2001
!   Modified by Agnes Helmstetter 2002-2011
!     References: 
!   - Felzer, K. R., T. W. Becker, R. E. Abercrombie, G. Ekstro¬ m, and J. R. Rice,
!     Triggering of the 1999 MW 7.1 Hector Mine earthquake by aftershocks of the 1992 MW 7.3 Landers earthquake, 
!     J. Geophys. Res., 107(B9), 2190, doi:10.1029/2001JB000911, 2002. 
!     - Helmstetter, A. and D. Sornette, 
!     Sub-critical and supercritical regimes in epidemic models of earthquake aftershocks, 
!     J. Geophys. Res., 107, doi:10.1029/2001JB001580, 2002 

!   Compilation: g77 ETAS.f -O3 -o ETAS -lm

     !subroutine ETAS(q)

      program ETAS3D
      implicit none
      integer q,MS,lcat,laft2,ii,jj,iseed,i,j,it,NS,nsi,NMAX,lcatold
      parameter(MS = 7000000)       ! max number of aftershocks
      real M0,M1,Md,Mc,Mmax,K,b,p,c,alpha,tmax,MB,d,mu,gamma,drho,zmax
      double precision mubg,R,r1,r2,r3,r4,r5,tbg,mbg,n,tstar,NDIR
      double precision t1(MS),t2(MS),t(MS),ran3
      real mag(MS),mag1(MS),mag2(MS),x(MS),y(MS),z(MS)
      real x1(MS),y1(MS),z1(MS),x2(MS),y2(MS),z2(MS)
      integer ind1(MS),ind2(MS),ims(MS)
      integer g(MS),naft(MS),imsi(MS),nafti(MS)
      character*100 outputfile
      integer rank(MS)
      external ran3
      common/CAT/t,x,y,z,mag,ims,g,naft,imsi,nafti
      common/CATA/t1,x1,y1,z1,mag1,ind1,t2,x2,y2,z2,mag2,ind2
      common/ETAS/M0,M1,Md,Mmax,MB,K,b,p,c,alpha,drho,tmax,d,mu,gamma
      common/ETAS/zmax

 2004  format('% b= ',f4.2,' p= ',f5.3,' alpha= ',f5.3,' K= ',f9.7, &
        & ' c= ',f6.4)
 2005  format('% M0= ',f4.2,' M1= ',f4.2,' Md= ',f4.2,' Mc=',f4.2,' Mmax=',f4.2,' iseed= ',i12,' tmax= ',g8.2,' n=',f8.5)
 2006  format('% d= ',f5.3,' mu= ',f4.2,' gamma= ',f4.2,' mub= ',g8.2,' R= ',g8.2,' Zmax= ',g8.2) 
 2003  format(f24.8,f6.2,1x,3f11.5,6i8)  !  catalog format
 
! ------- read input parameters and writes parameter values in header of output file
      open(2,file='ETAS.par')
      read(2,*) b          ! GR exponent
      read(2,*) p          ! "local" Omori low (direct aftershocks) 
      read(2,*) c     ! characteristic time Omori law
      read(2,*) n     ! branching ratio = Kb/(alpha-b)/c^theta/theta
      read(2,*) alpha ! aft. productivity scaling with m
      read(2,*) d     ! characteristic triggering distance for M=M0
      read(2,*) mu    ! exponent decrease in aft density with distance from mainshock 
      read(2,*) gamma ! increase in triggering distance with mainshock magnitude:  d = md*10^(gamma*M)
      read(2,*) M0     ! min magnitude of triggering events
      read(2,*) M1     ! magnitude of first event
      read(2,*) Md     ! detection magnitude
      read(2,*) Mmax     ! maximum magnitude (tapered GR law between Md and Mmax)
      read(2,*) mubg  ! external background (# of EQs per unit time)
      read(2,*) R     ! spatial window 0<x,y<R
      read(2,*) zmax  ! max depth for background  and aft. 
      read(2,*) tmax  ! max duration of catalog
      read(2,*) MB     ! mag. min. screen output     
      read(2,*) Mc     ! mag. min. output catalog
      read(2,*) outputfile  ! output file for catalog (t,m,x,y,z)
      read(2,*) NMAX  ! max number of events
      read(2,*) iseed ! integer, seed to initialize random number generator
      close(2)
   
       K=n/(b/(b-alpha)/(p-1.)/c**(p-1.)*(1.-10.**((alpha-b)*(Mmax-M0))))
      if (alpha.eq.b) K=n/(b*log(10.)*(Mmax-M0)/(p-1.)/c**(p-1.)) ! if alpha=b the above equation is incorrect
      K=K*10.**(-b*(Md-M0))                              ! aftershock productivity. 
      tstar=c*(n/(1.-n))**(1./(p-1))                    ! characteristic time of aftershock decay, see Helmstetter and Sornette JGR 2002
      NDIR=K*10**(alpha*(M1-M0))/(p-1.)/c**(p-1.) ! number of direct m>Md aftershocks
      drho=K*b/(b-alpha)*10.**(b*(Md-M0))*(1.-10.**((alpha-b)*(Md-M0))) ! coorrection for missing undetected events
      if (n.gt.1) then
           write(*,*) '!! WARNING: branching number n=',n,' When n>1, there may be an infinite number of earthquakes!'
      endif
      open(3,file=outputfile)
      write(*,2004) b,p,alpha,K,c
      write(3,2004) b,p,alpha,K,c
      write(*,2005) M0,M1,Md,Mc,Mmax,iseed,tmax,n
      write(3,2005) M0,M1,Md,Mc,Mmax,iseed,tmax,n     
      write(*,2006) d,mu,gamma,mubg,R,zmax
      write(3,2006) d,mu,gamma,mubg,R,zmax
!---- starts simulation of ETAS catalog     
      it=1   ! number of generations
!     initialisation
      t(1) = 0.      ! time
      mag(1) = M1    ! magnitude
      ims(1) = 0     ! ancetre
      g(1)   = 0     ! generation
      x(1) = R/2.    ! position (x,
      y(1) = R/2.    !      ! y)
      z(1)=zmax/2.   ! z
      
!   'first generation' 
      t1(1) = 0.
      mag1(1) = M1
      ind1(1) = 1
      x1(1) = R/2.
      y1(1) = R/2.
      z1(1) = zmax/2.
      ims(1)=0
      imsi(1)=0
      i=2
      write(*,'(a,a)') '%             time     mag      x       y        z         rank     ims     imsi generation'
!--------- tectonic loading
      if (mubg.gt.0.) then
      write(*,*) ' *** Tectonic loading (list of events with mag>',MB,')'
        tbg=0.
        do while(tbg.lt.tmax) 
           r1=ran3(iseed)
           tbg=tbg-log(r1)/mubg
           r2=ran3(iseed)
           mbg = M0 - log10(r2)/b
           r3=ran3(iseed)*R
           r4=ran3(iseed)*R
           r5=ran3(iseed)*zmax        
           t(i)=tbg
           mag(i)=mbg        
           ims(i)=0          ! mainshock index
           imsi(i)=0
           g(i)=0
           x(i)=r3
           y(i)=r4
           z(i)=r5
           x1(i)=r3
           y1(i)=r4
           z1(i)=r5
           mag1(i)=mbg
           t1(i)=tbg
           ind1(i)=i
            if (mag(i).gt.MB) write(*,2003) t(i),mag(i),x(i),y(i),z(i),i,ims(i),imsi(i),g(i)
           i=i+1
        enddo
         lcat=i        
          write(*,*) ' *** N=',lcat-1,' it=0'
      endif   
      lcat=i     
      lcatold=i
      laft2=i


! -------- Time magnitude, and aftershock generation        
      write(*,*) ' *** Direct and indirect aftershocks (list of events with mag>',MB,')'
10     call Getbigaftershocks(lcat,laft2,iseed,it)
      do i=1, laft2-1
         x1(i)=x2(i); y1(i)=y2(i); z1(i)=z2(i);
         t1(i)=t2(i); mag1(i)=mag2(i); ind1(i)=ind2(i);
      enddo
      write(*,*) ' *** N=',lcat-1,' dN',lcat-lcatold,' it=',it     
      lcatold=lcat
      if(lcat.gt.NMAX)  write(*,*) 'ERROR: Too many events : N=',lcat,' NMAX=',NMAX
      if(laft2.eq.1.or.lcat.gt.NMAX)  goto 20   ! no aftershock 
          it=it+1
      goto 10       ! iteration, next generation
20     continue
     lcat=lcat-1   ! number of EQS 
!---------  reorder catalog by increasing time and write output catalog
      write(*,*) 'END N=',lcat,' N(M>Mc)=', &
      & lcat*10.**(-b*(Mc-M0)),' it=',it
      if (lcat.ge.Nmax-1) then
          write(*,*) 'WARNING: too many events generated N=',lcat,' Nmax=',Nmax
          write(*,*) '     try decreasing background rate mub or branching ratio n'
          write(*,*) '     or increase  Nmax (must be <7000000)'
          return
      endif
      write(*,*) 'Reorder catalog by increasing time'
      call tricat(lcat,rank) ! reorder catalog by increasing time
      write(*,*) 'Write output catalog: ',outputfile
         write(3,'(a,a)')' %          time     mag     x        y       z     rank     imd     imi  generation naftd nafti'
        do i=1,lcat        
           if (mag(i).ge.Mc.and.x(i).le.R.and.y(i).le.R.and.x(i).ge.0..and.y(i).ge.0.) then
               write(3,2003) t(i),mag(i),x(i),y(i),z(i),rank(i),ims(i),imsi(i),g(i),naft(i),nafti(i)
           endif
        enddo   
     close(3)	
     end

!--------------------------------------------------------------------------- 
      subroutine Getbigaftershocks(lcat,laft2,iseed,it)
!       generate next generation aftershocks
      implicit none
      integer MS,lcat,iseed,laft2,idum,i,j,it,ii
      parameter(MS = 7000000)
      real M0,M1,Md,Mmax,MB,K,b,p,c,alpha,tmax,d,mu,gamma,n
      real a2,rho,n1,n2,n3,zmax,M,rr,pi,lambda1,lambda2,lambda,rh,drho     
      double precision t1(MS),t2(MS),t(MS)
      double precision ran3,timeaft,timerem,TT,r
      real mag(MS),mag1(MS),mag2(MS),x(MS),y(MS),z(MS)
      real x1(MS),y1(MS),z1(MS),x2(MS),y2(MS),z2(MS)
      integer ind1(MS),ind2(MS),g(MS),ims(MS),naft(MS)
      INTEGER imsi(MS),nafti(MS)
      external ran3
      common/CAT/t,x,y,z,mag,ims,g,naft,imsi,nafti
      common/CATA/t1,x1,y1,z1,mag1,ind1,t2,x2,y2,z2,mag2,ind2
      common/ETAS/M0,M1,Md,Mmax,MB,K,b,p,c,alpha,drho,tmax,d,mu
      common/ETAS/gamma,zmax
2003     format(f24.8,f6.2,1x,3f11.5,6i8)  !  catalog format
      idum = laft2-1
      laft2 = 1
      pi=3.14159265
      do i=1, idum
      if (mag1(i).gt.Md) then
        timeaft = 0.0
        timerem = tmax - t1(i)
        rho = K*10**(alpha*(mag1(i)-M0)) + drho
        do while(timeaft.lt.timerem.and.lcat.lt.MS)
          r=ran3(iseed)
          if(p.eq.1.0) then
            TT = exp(log(c+timeaft) - (1/rho)*log(r)) - c - timeaft
            timeaft = timeaft + TT
          else 
             n3 = exp((rho/(1.0 - p))*((timeaft+c)**(1.0-p)))
             if(p.gt.1.0.and.r.lt.n3) then
                    timeaft = timerem + 1.0
               else
                    n1 = (timeaft + c)**(1.0 - p)
                    n2 = ((1.0 - p)*log(r)/rho)
                    TT = (n1 - n2)**(1.0/(1.0-p)) - c - timeaft
                    timeaft = timeaft + TT   
               end if
          end if
          if(timeaft.lt.timerem) then
             M=Mmax+1.
             do while (M.gt.Mmax)
                M = Md - log10(ran3(iseed))/b
             enddo
             ! aft-mainshock distance position rr (isotrop. linear  power-law pdf of distance):
             rr=d*((1.-ran3(iseed))**(-1./mu)-1.)*10**(gamma*mag1(i))
             ! aft. depth  0<z<zmax, angle lambda between horizontal and ms-aft
             lambda1=-asin(min(1.,z1(i)/rr)) ! min angle 
             lambda2=asin(min(1.,(zmax-z1(i))/rr)) ! max angle 
             lambda=lambda1+ran3(iseed)*(lambda2-lambda1) ! uniform angle between lambda1 and lambda2
             z(lcat) = z1(i) + rr*sin(lambda)   ! aft. depth
             rh=rr*cos(lambda) ! horizontal distance
             r=ran3(iseed)
             x(lcat) = x1(i) + rh*cos(2.*pi*r)  
             y(lcat) = y1(i) + rh*sin(2.*pi*r)
             t(lcat) = t1(i) + timeaft
             mag(lcat) = M
             ims(lcat) = ind1(i)    ! index of 'direct mainshock' (0 if background EQ)
             if (g(ind1(i)).eq.0) then ! direct mainshock is a background EQ
                imsi(lcat) = ind1(i) 
             else  ! eq i is a secondary aft.
                imsi(lcat)=imsi(ind1(i)) ! index of 'indirect mainshock' (first generation EQ)
             endif
             g(lcat) = it        ! generation  (0 if background EQ)
             naft(ims(lcat))  =naft(ims(lcat))+1   ! number of direct (only) aftershocks
             nafti(ims(lcat)) =nafti(ims(lcat))+1  ! number of direct ...
             if (ims(lcat).ne.imsi(lcat)) then
                nafti(imsi(lcat))=nafti(imsi(lcat))+1 ! ... and indirect aftershocks
             endif
             t2(laft2)=t(lcat); mag2(laft2)=mag(lcat)
             x2(laft2)=x(lcat); y2(laft2)=y(lcat); z2(laft2)=z(lcat);
             ind2(laft2)=lcat
             if (M.gt.MB) then
!                ii=ind1(i) ! 'mainshock'
!                write(*,2003) t(ii),mag(ii),x(ii),y(ii),z(ii),ii,ims(ii),imsi(ii),g(ii),naft(ii),nafti(ii)
                ii=lcat    ! 'aftershock'
                write(*,2003) t(ii),mag(ii),x(ii),y(ii),z(ii),ii,ims(ii),imsi(ii),g(ii)
             endif  
             lcat = lcat + 1
             laft2 = laft2 + 1
             if (lcat.gt.MS) write(*,*) 'PROBLEM: too many events N=',lcat
          endif
        enddo ! next "aftershock"
        endif ! if m>md
      enddo    ! next "mainshock"
      return     
      end
!-----------------------------------------------------------------------------
      FUNCTION RAN3(IDUM)
      INTEGER IDUM     
      integer MBIG,MSEED,MZ
      real*8 RAN3,FAC     
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      integer I,IFF,II,INEXT,INEXTP,K
      integer MJ,MK,MA(55)
      SAVE IFF,INEXT,INEXTP,MA     
      DATA IFF /0/     
 1    IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
      IFF=1
      MJ=MSEED-IABS(IDUM)
      MJ=MOD(MJ,MBIG)
      MA(55)=MJ
      MK=1
      DO I=1,54
       II=MOD(21*I,55)
       MA(II)=MK
       MK=MJ-MK
       IF(MK.LT.MZ) MK=MK+MBIG
       MJ=MA(II)
      enddo  
      DO K=1,4
       DO I=1,55
         MA(I)=MA(I)-MA(1+MOD(I+30,55))
         IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
       enddo
      enddo  
      INEXT=0
      INEXTP=31
      IDUM=1         !
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC     
      if ((RAN3.le.0).or.(RAN3.ge.1)) then
      goto 1
      endif
      RETURN
      END      
!--------------------------------------------------------------------------- 
      subroutine tricat(l,rank) ! shell method, num. rec. p323
!       sort the catalog by increasing time t(1:l)
      implicit none
      integer i,j,inc,l,r,MS,gg,imss,naftt,imsii,naftii
      parameter (MS=7000000)
      double precision t(MS),tt
      real x(MS),y(MS),z(MS),mag(MS),xx,yy,zz,mm
      integer g(MS),ims(MS),imsi(MS),rank(MS),naft(MS),nafti(MS)
      common/CAT/t,x,y,z,mag,ims,g,naft,imsi,nafti
      do i=1,l
         rank(i)=i
      enddo
      inc=1
 1      inc=3*inc+1
      if(inc.le.l) goto 1
 2      continue
      inc=inc/3
      do 11 i=inc+1,l
         tt=t(i); xx=x(i); yy=y(i); zz=z(i); mm=mag(i) 
         gg=g(i); imss=ims(i); naftt=naft(i); r=rank(i)
         imsii=imsi(i); naftii=nafti(i)
         j=i
 3         if(t(j-inc).gt.tt) then
            t(j)=t(j-inc); x(j)=x(j-inc); y(j)=y(j-inc); z(j)=z(j-inc) 
            mag(j)=mag(j-inc); g(j)=g(j-inc); ims(j)=ims(j-inc) 
            naft(j)=naft(j-inc); rank(j)=rank(j-inc)
            imsi(j)=imsi(j-inc); nafti(j)=nafti(j-inc) 
            j=j-inc
            if(j.le.inc) goto 4
            goto 3
         endif
 4         t(j)=tt; x(j)=xx; y(j)=yy; z(j)=zz; mag(j)=mm
         g(j)=gg; ims(j)=imss; naft(j)=naftt; rank(j)=r
         imsi(j)=imsii; nafti(j)=naftii;
 11      continue
      if(inc.gt.1) goto 2
      return
      !END
      !END subroutine ETAS1
      END
