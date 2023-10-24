c
c programa para calcular la correccion perturbativa al estado 1s2s del helio
c
      program he1s2s
      implicit real *8 (a-h,o-z)
      write(*,*)'dime z'
      read(*,*)z
      enp=-z**2/2*(1.d0+1.d0/2.d0**2)
      call perturbacion(z,ejs,eks,ejp,ekp)
      write(*,*)'energia sin perturbar',enp
      write(*,*)'1s2s energia singlete',enp+(ejs+eks)
      write(*,*)'1s2s energia triplete',enp+(ejs-eks)
      write(*,*)'1s2p energia singlete',enp+(ejp+ekp)
      write(*,*)'1s2p energia triplete',enp+(ejp-ekp)
      write(*,*)'J_21',ejp
      write(*,*)'J_20',ejs
      write(*,*)'K_21',ekp
      write(*,*)'K_20',eks
      end
      subroutine perturbacion(z,ejs,eks,ejp,ekp)
      implicit real *8 (a-h,o-z)
      parameter (nmax=10)
      dimension dj2(0:nmax,0:nmax),dk2(0:nmax,0:nmax)
      z6=z**6
      z8=z**8
      n=4
      alfa=2.d0*z
      beta=z
      call integ2(n,alfa,beta,dj2,dk2)
      ejs=dj2(1,2)-z*dj2(1,3)+z**2/4.d0*dj2(1,4)
     >  +dk2(2,1)-z*dk2(2,2)+z**2/4.d0*dk2(2,3)
      ejp=dj2(1,4)+dk2(2,3)
      ejs=2.d0*z6*ejs
      ejp=z8/6.d0*ejp
      alfa=3.d0*z/2.d0
      beta=alfa
      call integ2(n,alfa,beta,dj2,dk2)
      eks=dj2(1,2)-z/2.d0*dj2(1,3)-z/2.d0*dj2(2,2)+z**2/4.d0*dj2(2,3)
     >  +dk2(2,1)-z/2.d0*dk2(2,2)-z/2.d0*dk2(3,1)+z**2/4.d0*dk2(3,2)
      ekp=dj2(1,4)+dk2(4,1)
      eks=2.d0*z6*eks
      ekp=z8/18.d0*ekp
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
