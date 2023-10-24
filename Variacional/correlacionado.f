c programa para helio variacional 
c psi=exp(-alfa*(r_1+r_2))(1+b*exp(-beta_i*r_{12}))
c
      program helio
      implicit real *8 (a-h,o-z)
      parameter (nmf=40)
      parameter (nsmax=20)
      common /constantes / pi
      common /nfac/fa(0:nmf)
      common /ion /znuc
      dimension plibre(nsmax),dplibre(nsmax)
c
      pi=4.d0*datan(1.d0)
      call factorial
c
      write(*,*) 'carga nuclear'
      read(*,*) znuc
      write(*,*) 'funcion prueba'
      write(*,*) 'exp(-alfa*(r_1+r_2))*(1+b*exp(-beta*r_{12})'
      write(*,*) 'parametro alfa (Z-5/16)'
      read(*,*) alfa
      write(*,*) 'parametros b y beta (0.2 y 0.2)'
      read(*,*) b,beta
      write(*,*) 'carga nuclear',znuc
      call energia(alfa,b,beta,etot)
      write(*,*) 'parametro alfa inicial',alfa
      write(*,*) 'parametros b y beta inicial',b,beta
      write(*,*) 'la energia inicial es', etot
c
c simplex
c
      ndim=3
      plibre(1)=alfa
      plibre(2)=b
      plibre(3)=beta
      dplibre(1)=0.2d0
      dplibre(2)=0.4d0
      dplibre(3)=0.4d0
      ftol=1.d-6
      call simplex(plibre,dplibre,fmin,ftol,iter,ndim)
      write(*,*)'minimo',fmin
      write(*,*) (plibre(ic),ic=1,3)
c
      end
c
      subroutine func(x,fx,ndim)
      implicit real *8 (a-h,o-z)
      dimension x(ndim)
      alfa=abs(x(1))
      b=x(2)
      beta=abs(x(3))
      call energia(alfa,b,beta,etot)
      fx=etot
      return
      end
      subroutine energia(alfa,b,beta,etot)
      implicit real *8 (a-h,o-z)
      parameter (nint=20,nmf=40)
      parameter (nmcor=2,nmpot=2)
      common /constantes / pi
      common /nfac/fa(0:nmf)
      common /ion /znuc
      common /integ/xi(-1:nint,-1:nint,0:nint)
      dimension bcor(nmcor),betacor(nmcor)
      dimension dnt(-1:nmpot,-1:nmpot,-1:nmpot,nmcor,nmcor)
c
      a=2.d0*alfa
      ncor=2
      bcor(1)=1.d0
      bcor(2)=b
      betacor(1)=0.d0
      betacor(2)=beta
c
      a=2.d0*alfa
      do 10 icor=1,ncor
        do 20 jcor=icor,ncor
          c=betacor(icor)+betacor(jcor)
          call isrk(a,a,c)
          do 30 i12=-1,nmpot
            do 40 i2=-1,nmpot
              do 50 i1=-1,nmpot
                dnt(i1,i2,i12,jcor,icor)=xi(i1+1,i2+1,i12+1)
   50         continue
   40       continue
   30     continue
   20   continue
   10 continue
      do 60 icor=1,ncor
        do 70 jcor=icor+1,ncor
          do 80 i12=-1,nmpot
            do 90 i2=-1,nmpot
              do 100 i1=-1,nmpot
                dnt(i1,i2,i12,icor,jcor)=dnt(i1,i2,i12,jcor,icor)
  100         continue
   90       continue
   80     continue
   70   continue
   60 continue
c
      alfa2=alfa**2
      dnor=0.d0
      datr=0.d0
      drep=0.d0
      dkin=0.d0
      do 110 icor=1,ncor
        do 120 jcor=1,ncor
          cij=bcor(icor)*bcor(jcor)
          betaj=betacor(jcor)
          betaj2=betaj**2
          albe=alfa*betaj
          dnor=dnor+cij*dnt(0,0,0,jcor,icor)
          datr=datr+cij*(dnt(-1,0,0,jcor,icor)+dnt(0,-1,0,jcor,icor))
          drep=drep+cij*dnt(0,0,-1,jcor,icor)
          dl1=alfa2*dnt(0,0,0,jcor,icor)
     >       -2.d0*alfa*dnt(-1,0,0,jcor,icor)
          dl2=alfa2*dnt(0,0,0,jcor,icor)
     >       -2.d0*alfa*dnt(0,-1,0,jcor,icor)
          dl12=betaj2*dnt(0,0,0,jcor,icor)
     >        -2.d0*betaj*dnt(0,0,-1,jcor,icor)
          dg112=albe*(dnt(1,0,-1,jcor,icor)+dnt(-1,0,1,jcor,icor)
     >               -dnt(-1,2,-1,jcor,icor))
          dg212=albe*(dnt(0,1,-1,jcor,icor)+dnt(0,-1,1,jcor,icor)
     >               -dnt(2,-1,-1,jcor,icor))
          dkin=dkin+cij*(dl1+dl2+2.d0*dl12+dg112+dg212)
  120   continue
  110 continue
      etot=(-0.5d0*dkin-znuc*datr+drep)/dnor
      return
      end
c
c calculo de las integrales
c \int d\vec{r}_1 d\vec{r}_2 exp(-a*r_1)exp(-b*r_2)exp(-c*r_{12})*
c                            r_1*(n_1-1) r_2*(n_2-1) r_{12}*(n_{12}-1)
c  Sack, Roothaan y Kolos, J. Math. Phys., 8, 1093 (1967)
c
      subroutine isrk(a,b,c)
      implicit real *8 (a-h,o-z)
      parameter (nint=20,nmf=40)
      dimension xa(0:nint,0:nint,0:nint),xb(-1:nint,0:nint,-1:nint)
      common /constantes / pi
      common /nfac/fa(0:nmf)
      common /integ/xi(-1:nint,-1:nint,0:nint)
c
      amb=1.d0/(a+b)
      amc=1.d0/(a+c)
      bmc=1.d0/(b+c)
c
      xgam=1.d0
      do 10 m=0,nint
        xgam=xgam*bmc   
        xn=xgam
        do 20 n=0,nint
          ns=n+m
          xa(0,m,n)=fa(ns)*xn*16.d0*pi*pi
          xn=xn*bmc
          do 30 l=1,nint
            xa(l,m,n)=0.d0
   30     continue
          xb(-1,n,m)=0.d0
          xb(m,n,-1)=0.d0
          xi(-1,m,n)=0.d0
          xi(m,-1,n)=0.d0
   20   continue
   10 continue
      
      do 40 l=0,nint
        do 50 n=0,nint
          do 60 m=0,nint
            xb(l,m,n)=(l*xb(l-1,m,n)+n*xb(l,m,n-1)+xa(l,m,n))*amc
   60 continue
   50 continue
   40 continue

      do 70 l=0,nint
        do 80 n=0,nint
          do 90 m=0,nint
            xi(l,m,n)=(l*xi(l-1,m,n)+m*xi(l,m-1,n)+xb(l,m,n))*amb
   90     continue
   80   continue
   70 continue
      return
      end
c
      subroutine factorial
      implicit real *8 (a-h,o-z)
      parameter (nmf=40)
      common /nfac/ fa(0:nmf)
      fa(0)=1.d0
      do 10 n=1,nmf
      fa(n)=n*fa(n-1)
 10   continue
      return
      end
c
      subroutine simplex(x,dx,fmin,ftol,iter,ndim)
      implicit real *8 (a-h,o-z)
      dimension p(ndim+1,ndim),y(ndim+1)
      dimension x(ndim),aux(ndim),dx(ndim)
      mp=ndim+1
      np=ndim
      do 10 j=1,ndim+1
      do 20 i=1,ndim
      if(i.eq.j) then
      p(j,i)=x(i)+dx(i)
      else
      p(j,i)=x(i)
      endif
      aux(i)=p(j,i)
   20 continue
      call func(aux,y(j),ndim)
   10 continue
      call amoeba(p,y,mp,np,ndim,ftol,iter)
      fmin=y(1)
      do 30 i=1,ndim
      x(i)=p(1,i)
   30 continue
      iter=iter+ndim
      return
      end
      subroutine amoeba(p,y,mp,np,ndim,ftol,iter)
      implicit real *8 (a-h,o-z)
      parameter (nmax=50,itmax=5000)
      parameter (eps=1.d-10)
      dimension p(mp,np),y(mp)
      dimension psum(nmax)
      iter=0
    1 do 10 i=1,ndim
      psum(i)=0.d0
      do 20 j=1,ndim+1
      psum(i)=psum(i)+p(j,i)
   20 continue
   10 continue
    2 ilo=1
      if(y(1).gt.y(2)) then
      ihi=1
      inhi=2
      else
      ihi=2
      inhi=1
      endif
      do 30 i=1,ndim+1
      if(y(i).le.y(ilo)) ilo=i
      if(y(i).gt.y(ihi)) then
      inhi=ihi
      ihi=i
      elseif(y(i).gt.y(inhi)) then
      if(i.ne.ihi) inhi=i
      endif
   30 continue
c     write(*,*)y(ihi),y(inhi),y(ilo)
      rtol=2.d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+eps)
      if(rtol.lt.ftol) then
      swap=y(1)
      y(1)=y(ilo)
      y(ilo)=swap
      do 40 i=1,ndim
      swap=p(1,i)
      p(1,i)=p(ilo,i)
      p(ilo,i)=swap
   40 continue
      return
      endif
      if(iter.gt.itmax) then
      write(*,*)'itmax exceeded in amoeba'
      stop
      endif
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,ihi,-1.d0)
      if(ytry.le.y(ilo)) then
       ytry=amotry(p,y,psum,mp,np,ndim,ihi,2.d0)
      elseif(ytry.ge.y(inhi)) then
       ysave=y(ihi)
       ytry=amotry(p,y,psum,mp,np,ndim,ihi,0.5d0)
       if(ytry.ge.ysave) then
        do 50 i=1,ndim+1
        if(i.ne.ilo) then
         do 60 j=1,ndim
         psum(j)=0.5d0*(p(i,j)+p(ilo,j))
         p(i,j)=psum(j)
   60    continue
         call func(psum,y(i),ndim)
        endif
   50   continue
        iter=iter+ndim
        goto 1
       endif
      else
       iter=iter-1
      endif
      goto 2
      end
      double precision function amotry(p,y,psum,mp,np,ndim,ihi,fac)
      implicit real *8 (a-h,o-z)
      parameter (nmax=50)
      dimension p(mp,np),y(mp),psum(np)
      dimension ptry(nmax)
      fac1=(1.d0-fac)/ndim
      fac2=fac1-fac
      do 10 j=1,ndim
      ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
   10 continue
      call func(ptry,ytry,ndim)
      if(ytry.lt.y(ihi)) then
       y(ihi)=ytry
       do 20 j=1,ndim
       psum(j)=psum(j)-p(ihi,j)+ptry(j)
       p(ihi,j)=ptry(j)
   20  continue
      endif
      amotry=ytry
      return
      end
