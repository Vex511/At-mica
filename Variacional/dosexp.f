c
c programa para calcular variacionalmente el estado fundamental del
c atomo de dos electrones con una funcion prueba
c psi=A*exp(-alfa*r)+B*exp(-beta*r)
c
      program hevargs
      implicit real *8 (a-h,o-z)
      parameter (nsmax=20)
      dimension plibre(nsmax),dplibre(nsmax)
      common /ion /znuc
c     a=2.60505d0
c     b=2.08144d0
c     alfa=1.41d0
c     beta=2.61d0
c     z=2
      write(*,*) 'carga nuclear'
      read(*,*) znuc
      write(*,*) 'funcion prueba'
      write(*,*) 'Psi(r_1,r_2)=u(r_1)*u(r_2)'
      write(*,*) 'u(r)=A*exp(-alfa*r)+B*exp(-beta*r)'
      write(*,*) 'parametros alfa y beta (1.4 y 2.6)'
      read(*,*) alfa,beta
      write(*,*) 'parametros A y B (2.6 y 2.1)'
      read(*,*) a,b
      call energia(a,b,alfa,beta,ener,dnor)
      write(*,*)'la energia inicial es',ener
c
c simplex
c
      ndim=4
      plibre(1)=alfa
      plibre(2)=beta
      plibre(3)=a
      plibre(4)=b
      dplibre(1)=0.2d0
      dplibre(2)=0.2d0
      dplibre(3)=0.2d0
      dplibre(4)=0.2d0
      ftol=1.d-6
      call simplex(plibre,dplibre,fmin,ftol,iter,ndim)
      alfa=plibre(1)
      beta=plibre(2)
      a=plibre(3)
      b=plibre(4)
      call energia(a,b,alfa,beta,ener,dnor)
      a=a/dnor**0.25d0
      b=b/dnor**0.25d0
      write(*,*)'alfa,beta',alfa,beta
      write(*,*)'a,b',a,b
      write(*,*) a**2/(2.0*alfa)**2+b**2/(2.0*beta)**2
     >                     +2.0*a*b/(alfa+beta)**2
      call energia(a,b,alfa,beta,ener,dnor)
      write(*,*)'ener',ener
      write(*,*)'dnor',dnor
      end
      subroutine func(x,fx,ndim)
      implicit real *8 (a-h,o-z)
      dimension x(ndim)
      alfa=abs(x(1))
      beta=abs(x(2))
      a=x(3)
      b=x(4)
      call energia(a,b,alfa,beta,ener,dnor)
      fx=ener
      return
      end
      subroutine energia(a,b,alfa,beta,ener,dnor)
      implicit real *8 (a-h,o-z)
      parameter (nmax=10)
      common /ion /znuc
      dimension di1aa(0:2*nmax),di1bb(0:2*nmax),di1ab(0:2*nmax)
      dimension dj22a2a(0:nmax,0:nmax),dk22a2a(0:nmax,0:nmax)
      dimension dj22aab(0:nmax,0:nmax),dk22aab(0:nmax,0:nmax)
      dimension dj22a2b(0:nmax,0:nmax),dk22a2b(0:nmax,0:nmax)
      dimension dj22b2b(0:nmax,0:nmax),dk22b2b(0:nmax,0:nmax)
      dimension dj22bab(0:nmax,0:nmax),dk22bab(0:nmax,0:nmax)
      dimension dj22b2a(0:nmax,0:nmax),dk22b2a(0:nmax,0:nmax)
      dimension dj2ab2a(0:nmax,0:nmax),dk2ab2a(0:nmax,0:nmax)
      dimension dj2abab(0:nmax,0:nmax),dk2abab(0:nmax,0:nmax)
      dimension dj2ab2b(0:nmax,0:nmax),dk2ab2b(0:nmax,0:nmax)
      n=4
      alfa2=2.d0*alfa
      beta2=2.d0*beta
      albe=alfa+beta
      a2=a**2
      b2=b**2
      ab=a*b
      call integ1(n,alfa2,di1aa)
      call integ1(n,beta2,di1bb)
      call integ1(n,albe,di1ab)
      dnor=(a2*di1aa(2)+b2*di1bb(2)+2.d0*ab*di1ab(2))**2
      datr=(a2*di1aa(1)+b2*di1bb(1)+2.d0*ab*di1ab(1))
     >    *(a2*di1aa(2)+b2*di1bb(2)+2.d0*ab*di1ab(2))
      dkin=(a2*(alfa**2*di1aa(2)-alfa2*di1aa(1))
     >     +b2*(beta**2*di1bb(2)-beta2*di1bb(1))
     >     +ab*((alfa**2+beta**2)*di1ab(2)-2.d0*albe*di1ab(1)))
     >    *(a2*di1aa(2)+b2*di1bb(2)+2.d0*ab*di1ab(2))
      n=2
      call integ2(n,2.d0*alfa,2.d0*alfa,dj22a2a,dk22a2a)
      call integ2(n,2.d0*alfa,alfa+beta,dj22aab,dk22aab)
      call integ2(n,2.d0*alfa,2.d0*beta,dj22a2b,dk22a2b)
      call integ2(n,2.d0*beta,2.d0*beta,dj22b2b,dk22b2b)
      call integ2(n,2.d0*beta,alfa+beta,dj22bab,dk22bab)
      call integ2(n,2.d0*beta,2.d0*alfa,dj22b2a,dk22b2a)
      call integ2(n,alfa+beta,2.d0*alfa,dj2ab2a,dk2ab2a)
      call integ2(n,alfa+beta,alfa+beta,dj2abab,dk2abab)
      call integ2(n,alfa+beta,2.d0*beta,dj2ab2b,dk2ab2b)
      drep=
     >     a**4*dj22a2a(1,2)
     >    +2.d0*a**3*b*dj22aab(1,2)+a**2*b**2*dj22a2b(1,2)
     >    +2.d0*a**3*b*dj2ab2a(1,2)+4.d0*a**2*b**2*dj2abab(1,2)
     >    +2.d0*a*b**3*dj2ab2b(1,2)
     >    +a**2*b**2*dj22b2a(1,2)+2.d0*a*b**3*dj22bab(1,2)
     >    +b**4*dj22b2b(1,2)
     >    +a**4*dk22a2a(2,1)
     >    +2.d0*a**3*b*dk22aab(2,1)+a**2*b**2*dk22a2b(2,1)
     >    +2.d0*a**3*b*dk2ab2a(2,1)+4.d0*a**2*b**2*dk2abab(2,1)
     >    +2.d0*a*b**3*dk2ab2b(2,1)
     >    +a**2*b**2*dk22b2a(2,1)+2.d0*a*b**3*dk22bab(2,1)
     >    +b**4*dk22b2b(2,1)
      datr=2.d0*(-znuc)*datr
      dkin=2.d0*(-1.d0/2.d0)*dkin
      ener=(datr+dkin+drep)/dnor
      return
      end
      subroutine integ2(n,alfa,beta,dj2,dk2)
      implicit real *8 (a-h,o-z)
      parameter (nmax=10)
      dimension dj2(0:nmax,0:nmax),dk2(0:nmax,0:nmax)
      dimension di1(0:2*nmax)
      if(alfa.le.0.d0.or.beta.le.0.d0) then
        write(*,*)'error en integ1 alfa o beta cero o negativo'
        write(*,*)'alfa',alfa
        write(*,*)'beta',beta
        stop
      endif
      albe=alfa+beta
      call integ1(n,albe,di1)
      dj2(0,0)=1.d0/(alfa*albe)
      dk2(0,0)=1.d0/(beta*albe)
      do in=1,n
          dj2(in,0)=(in*dj2(in-1,0)+di1(in))/alfa
          dk2(in,0)=(in*dk2(in-1,0)-di1(in))/alfa
        enddo
      do im=1,n
        dj2(0,im)=(im*dj2(0,im-1)-di1(0+im))/beta
        dk2(0,im)=(im*dk2(0,im-1)+di1(0+im))/beta
        do in=1,n
          dj2(in,im)=(in*dj2(in-1,im)+di1(in+im))/alfa
          dk2(in,im)=(in*dk2(in-1,im)-di1(in+im))/alfa
        enddo
      enddo
      return
      end
      subroutine integ1(n,alfa,di1)
      implicit real *8 (a-h,o-z)
      parameter (nmax=10)
      dimension di1(0:nmax)
      if(alfa.le.0.d0) then
        write(*,*)'error en integ1 alfa cero o negativo'
        write(*,*)'alfa',alfa
        stop
      endif
      di1(0)=1.d0/alfa
      do i=1,n
        di1(i)=i*di1(i-1)/alfa
      enddo
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
