      program excithe
      implicit real *8 (a-h,o-z)
      common /sistema / znuc,sesp,na,la,nda,nb,lb,ndb
      common /ligadura /optdos
      data  htoev/ 27.2113846d0/
      logical optdos
c     write(*,*)'Valor de la carga nuclear'
c     read(*,*)znuc
      znuc=2.d0
      write(*,*)'Valores de n l'
      read(*,*)ncuan,lcuan
      write(*,*)'Valor del espin total (0 o 1)'
      read(*,*)iespin
      if(ncuan.eq.1.and.lcuan.eq.0.and.iespin.ne.0) then
        write(*,*)'no existe un triplete de espin con n=0 y L=0'
        stop
      endif
      write(*,*)'0 los dos alfas iguales para el estado fundamental'
      read(*,*)iexpon
      optdos=.true.
      if(ncuan.eq.1.and.lcuan.eq.0.and.iexpon.eq.0)optdos=.false.

         
      if(iespin.eq.0) then
         sesp=1.d0
      else
        sesp=-1.d0
      endif
      na=1
      la=0
      nb=ncuan
      lb=lcuan
      nda=na-la
      ndb=nb-lb
      call perturbativo(eindp,epert)
      call variacional(ev1,ev2)
      write(*,*)
      write(*,*)'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      write(*,*)'**********************************************'
      write(*,*)'|                  Sistema                   |'
      write(*,*)'**********************************************'
      write(*,*)'carga nuclear',znuc
      write(*,*)'n,l          ',ncuan,lcuan
      if(iespin.eq.0) then
        write(*,*)'Calculo para el singlete de espin'
      else
        write(*,*)'Calculo para el triplete de espin'
      endif
      write(*,*)'energia particula indep    ',eindp
      write(*,*)'energia perturbativa       ',epert
!     write(*,*)'energia variacional inicial',ev1
      write(*,*)'energia variacional optima ',ev2
      write(*,*)'--------------Resultados en eV--------------'
      write(*,*)'energia particula indep    ',eindp*htoev
      write(*,*)'energia perturbativa       ',epert*htoev
!     write(*,*)'energia variacional inicial',ev1*htoev
      write(*,*)'energia variacional optima ',ev2*htoev
      write(*,*)'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

      end
      subroutine variacional(ev1,ev2)
      implicit real *8 (a-h,o-z)
      parameter (nimax=20)
      parameter (nprop=20)
      parameter (nsmax=20)
      dimension plibre(nsmax),dplibre(nsmax)
      dimension otrose(nprop)
      common /sistema / znuc,sesp,na,la,nda,nb,lb,ndb
      common /coeflin/ cfrad(nimax,2)
      common /ligadura /optdos
      logical optdos
      data  htoev/ 27.2113846d0/
      alfa=znuc/na
      alfb=znuc/nb
      call energia(alfa,alfb,ener,otrose)
      eunc=otrose(1)
      edosc=otrose(2)
      e2dir=otrose(3)
      e2exc=otrose(4)
      itcon=otrose(5)
      enea=otrose(6)
      eneb=otrose(7)
      sumab=otrose(8)
      difer=abs((enea-eneb)/(enea+eneb))
      pdir=abs(e2dir)/(abs(e2dir)+abs(e2exc))*100.d0
      pexc=abs(e2exc)/(abs(e2dir)+abs(e2exc))*100.d0
      write(*,*)'**********************************************'
      write(*,*)'            Calculo Variacional               '
      write(*,*)'   Energia con los parametros iniciales       '
      write(*,*)'**********************************************'
      write(*,*)'-----Exponentes de la funcion de onda------'
      write(*,*)'alfa inicial                      ',alfa
      write(*,*)'alfb inicial                      ',alfb
      write(*,*)'Z alfa efectiva inicial (alfa)    ',alfa
      write(*,*)'Z alfb efectiva inicial (alfb*n)  ',alfb*nb
      write(*,*)'-----------Resultados en hartree-----------'
      write(*,*)'energia a un cuerpo   (h)  <h1>   ',eunc
      write(*,*)'energia a dos cuerpos   (h) <h12> ',edosc
      write(*,*)'componente directa     (h)  Jnl   ',e2dir
      write(*,*)'componente intercambio (h) Knl    ',e2exc
      write(*,*)'energia total       (h)           ',ener
      write(*,*)'--------------Resultados en eV--------------'
      write(*,*)'energia a un cuerpo   (eV)  <h1>   ',eunc*htoev
      write(*,*)'energia a dos cuerpos   (eV) <h12> ',edosc*htoev
      write(*,*)'energia total       (eV)           ',ener*htoev
      write(*,*)'-----------------Porcentajes-----------------'
      write(*,*)'100*(eprima/ecero)',100.d0*abs(edosc/eunc)
      write(*,*)'porcentaje directo     ',pdir
      write(*,*)'porcentaje intercambio ',pexc
      write(*,*)'---------Orbitales monoparticulares----------'
      write(*,*)'Exponentes             ',alfa,alfb
      write(*,*)'Coefientes primero     ',(cfrad(ic,1),ic=1,nda)
      write(*,*)'Coefientes segundo     ',(cfrad(ic,2),ic=1,ndb)
      write(*,*)'ortogonalidad orbitales',sumab
      write(*,*)'---------Datos tecnicos del calculo----------'
      write(*,*)'iteraciones hasta autoconsistencia',itcon
      write(*,*)'autovalor monoparticular Ea       ',enea
      write(*,*)'autovalor monoparticular Eb       ',eneb
      write(*,*)'diferencia relativa entre ellos   ',difer
      write(*,*)
      ev1=ener
      write(20,*) na,la,nda
      write(20,*) alfa
      do ic=1,nda
        write(20,*) cfrad(ic,1)
      enddo
      write(20,*) nb,lb,ndb
      write(20,*) alfb
      do ic=1,ndb
        write(20,*) cfrad(ic,2)
      enddo
      write(20,*) sesp
      write(20,*) znuc
      write(20,*) ener
c
c simplex
c
      if(optdos) then
        nsimplex=2
        plibre(1)=alfa
        plibre(2)=alfb
      else
        nsimplex=1
        plibre(1)=alfa
      endif
      do ic=1,nsimplex
        dplibre(ic)=plibre(ic)/10.d0
      enddo
      ftol=1.d-8
      call simplex(plibre,dplibre,fmin,ftol,iter,nsimplex)
      write(*,*)'optimizacion finalizada'
      write(*,*)'tolerancia optimizacion',ftol
      if(optdos) then
        alfa=plibre(1)
        alfb=plibre(2)
      else
        alfa=plibre(1)
        alfb=alfa
      endif
      write(*,*)'parametros optimos',alfa,alfb
      call energia(alfa,alfb,ener,otrose)
      eunc=otrose(1)
      edosc=otrose(2)
      e2dir=otrose(3)
      e2exc=otrose(4)
      itcon=otrose(5)
      enea=otrose(6)
      eneb=otrose(7)
      sumab=otrose(8)
      difer=abs((enea-eneb)/(0.5d0*(enea+eneb)))
      pdir=abs(e2dir)/(abs(e2dir)+abs(e2exc))*100.d0
      pexc=abs(e2exc)/(abs(e2dir)+abs(e2exc))*100.d0
      write(*,*)'**********************************************'
      write(*,*)'            Calculo Variacional               '
      write(*,*)'   Energia con los parametros optimos         '
      write(*,*)'**********************************************'
      write(*,*)'-----Exponentes de la funcion de onda------'
      write(*,*)'alfa optimo                       ',alfa
      write(*,*)'alfb optimo                       ',alfb
      write(*,*)'Z alfa efectiva optima (alfa)     ',alfa
      write(*,*)'Z alfb efectiva optima (alfb*n)   ',alfb*nb
      write(*,*)'-----------Resultados en hartree-----------'
      write(*,*)'energia a un cuerpo   (h)  <h1>   ',eunc
      write(*,*)'energia a dos cuerpos   (h) <h12> ',edosc
      write(*,*)'componente directa     (h)  Jnl   ',e2dir
      write(*,*)'componente intercambio (h) Knl    ',e2exc
      write(*,*)'energia total       (h)           ',ener
      write(*,*)'--------------Resultados en eV--------------'
      write(*,*)'energia a un cuerpo   (eV)  <h1>   ',eunc*htoev
      write(*,*)'energia a dos cuerpos   (eV) <h12> ',edosc*htoev
      write(*,*)'energia total       (eV)           ',ener*htoev
      write(*,*)'-----------------Porcentajes-----------------'
      write(*,*)'100*(eprima/ecero)',100.d0*abs(edosc/eunc)
      write(*,*)'porcentaje directo     ',pdir
      write(*,*)'porcentaje intercambio ',pexc
      write(*,*)'---------Orbitales monoparticulares----------'
      write(*,*)'Exponentes             ',alfa,alfb
      write(*,*)'Coefientes primero     ',(cfrad(ic,1),ic=1,nda)
      write(*,*)'Coefientes segundo     ',(cfrad(ic,2),ic=1,ndb)
      write(*,*)'ortogonalidad orbitales',sumab
      write(*,*)'---------Datos tecnicos del calculo----------'
      write(*,*)'iteraciones del simplex           ',iter
      write(*,*)'iteraciones hasta autoconsistencia',itcon
      write(*,*)'autovalor monoparticular Ea       ',enea
      write(*,*)'autovalor monoparticular Eb       ',eneb
      write(*,*)'diferencia relativa entre ellos   ',difer
      write(*,*)
      ev2=ener
      write(21,*) na,la,nda
      write(21,*) alfa
      do ic=1,nda
        write(21,*) cfrad(ic,1)
      enddo
      write(21,*) nb,lb,ndb
      write(21,*) alfb
      do ic=1,ndb
        write(21,*) cfrad(ic,2)
      enddo
      write(21,*) sesp
      write(21,*) znuc
      write(21,*) ener
      return
      end
      subroutine func(x,fx,ndim)
      implicit real *8 (a-h,o-z)
      parameter (nprop=20)
      dimension otrose(nprop)
      dimension x(ndim)
      common /ligadura /optdos
      logical optdos
      do ic=1,ndim
       x(ic)=abs(x(ic))
      enddo
      if(optdos) then
        alfa=x(1)
        alfb=x(2)
      else
        alfa=x(1)
        alfb=alfa
      endif
      call energia(alfa,alfb,ener,otrose)
      fx=ener
      return
      end
      subroutine energia(alfa,alfb,ener,otrose)
      implicit real *8 (a-h,o-z)
      parameter (nprop=20)
      parameter (nimax=20)
      dimension otrose(nprop)
      common /sistema / znuc,sesp,na,la,nda,nb,lb,ndb
      common /muncuerpo /dhaa(nimax,nimax),dnaa(nimax,nimax),
     >                   dhbb(nimax,nimax),dnbb(nimax,nimax),
     >                   dhab(nimax,nimax),dnab(nimax,nimax)
      common /mdoscuerpos /dh2d(nimax,nimax,nimax,nimax),
     >                     dh2e(nimax,nimax,nimax,nimax)
      dimension dhtot(nimax,nimax,nimax,nimax),
     >          dntot(nimax,nimax,nimax,nimax)
      dimension dha1(nimax,nimax),dna1(nimax,nimax),
     >          dhb1(nimax,nimax),dnb1(nimax,nimax)
      common /coeflin/ cfrad(nimax,2)
      data nmiter /50/
      data precision /1.d-10/
      call matrhn(znuc,alfa,alfb,na,nb,la,lb,nda,ndb)
      call matrh2(alfa,alfb,na,nb,la,lb,nda,ndb)
      cexc=sesp/(2.d0*lb+1.d0)
      if(sesp.lt.0.d0.and.lb.eq.0) then
        nautb=ndb-1
      else
        nautb=ndb
      endif
      do i4=1,ndb
        do i3=1,ndb
          do i2=1,nda
            do i1=1,nda 
              dn1=2.d0*(dnaa(i1,i2)*dnbb(i3,i4)
     >           +dnab(i1,i3)*dnab(i2,i4)*sesp)
              dh1=dhaa(i1,i2)*dnbb(i3,i4)+dhbb(i3,i4)*dnaa(i1,i2)
     >           +dhab(i1,i3)*dnab(i2,i4)*2.d0*sesp
              dh2=dh2d(i1,i2,i3,i4)+cexc*dh2e(i1,i2,i3,i4)
              dhtot(i1,i2,i3,i4)=2.d0*dh1+2.d0*dh2
              dntot(i1,i2,i3,i4)=dn1
            enddo
          enddo
        enddo
      enddo
      call diagon(dhaa,dnaa,nda,nda,enea,cfrad(1,1))
      call diagon(dhbb,dnbb,ndb,ndb,eneb,cfrad(1,2))
      difer=2.d0*precision
      itcon=0
      do iter=1,nmiter
        if(difer.gt.precision) then
          do j=1,nda
            do i=1,nda
              sume=0.d0
              sumn=0.d0
              do i4=1,ndb
                do i3=1,ndb
                  coef=cfrad(i3,2)*cfrad(i4,2)
                  sume=sume+coef*(dhtot(i,j,i3,i4)+dhtot(j,i,i3,i4))
                  sumn=sumn+coef*(dntot(i,j,i3,i4)+dntot(j,i,i3,i4))
                enddo
              enddo
              dha1(i,j)=sume
              dna1(i,j)=sumn
            enddo
          enddo
          do j=1,ndb
            do i=1,ndb
              sume=0.d0
              sumn=0.d0
              do i2=1,nda
                do i1=1,nda
                  coef=cfrad(i1,1)*cfrad(i2,1)
                  sume=sume+coef*(dhtot(i1,i2,i,j)+dhtot(i1,i2,j,i))
                  sumn=sumn+coef*(dntot(i1,i2,i,j)+dntot(i1,i2,j,i))
                enddo
              enddo
              dhb1(i,j)=sume
              dnb1(i,j)=sumn
            enddo
          enddo
          call diagon(dha1,dna1,nda,nda,enea,cfrad(1,1))
          call diagon(dhb1,dnb1,ndb,nautb,eneb,cfrad(1,2))
          difer=abs((enea-eneb)/(0.5d0*(enea+eneb)))
          itcon=itcon+1
        endif 
      enddo
      dnor=0.d0
      eunc=0.d0
      e2dir=0.d0
      e2exc=0.d0
      do i4=1,ndb
        do i3=1,ndb
          coefb=cfrad(i3,2)*cfrad(i4,2)
          do i2=1,nda
            do i1=1,nda
              coef=cfrad(i1,1)*cfrad(i2,1)*coefb
              dn1=2.d0*(dnaa(i1,i2)*dnbb(i3,i4)
     >           +dnab(i1,i3)*dnab(i2,i4)*sesp)
              dh1=dhaa(i1,i2)*dnbb(i3,i4)+dhbb(i3,i4)*dnaa(i1,i2)
     >           +dhab(i1,i3)*dnab(i2,i4)*2.d0*sesp
              dnor=dnor+coef*dn1
              eunc=eunc+coef*2.d0*dh1
              e2dir=e2dir+coef*2.d0*dh2d(i1,i2,i3,i4)
              e2exc=e2exc+coef*2.d0*dh2e(i1,i2,i3,i4)
            enddo
          enddo
        enddo
      enddo
      sumab=0.d0
      do i3=1,ndb
        do i1=1,nda
          sumab=sumab+cfrad(i1,1)*dnab(i1,i3)*cfrad(i3,2)
         enddo
      enddo
      eunc=eunc/dnor
      e2dir=e2dir/dnor
      e2exc=e2exc/dnor
      e2exc=sesp/(2.d0*lb+1.d0)*e2exc
      edosc=e2dir+e2exc
      ener=eunc+edosc
      otrose(1)=eunc
      otrose(2)=edosc
      otrose(3)=e2dir
      otrose(4)=e2exc
      otrose(5)=itcon
      otrose(6)=enea
      otrose(7)=eneb
      otrose(8)=sumab
      return
      end
      subroutine perturbativo(eunc,epert)
      implicit real *8 (a-h,o-z)
      parameter (nimax=20)
      common /sistema / znuc,sesp,na,la,nda,nb,lb,ndb
      common /muncuerpo /dhaa(nimax,nimax),dnaa(nimax,nimax),
     >                   dhbb(nimax,nimax),dnbb(nimax,nimax),
     >                   dhab(nimax,nimax),dnab(nimax,nimax)
      common /mdoscuerpos /dh2d(nimax,nimax,nimax,nimax),
     >                     dh2e(nimax,nimax,nimax,nimax)
      dimension cfrad(nimax,2)
      data  htoev/ 27.2113846d0/
      alfa=znuc/na
      alfb=znuc/nb
      call matrhn(znuc,alfa,alfb,na,nb,la,lb,nda,ndb)
      call matrh2(alfa,alfb,na,nb,la,lb,nda,ndb)
      call diagon(dhaa,dnaa,nda,nda,enea,cfrad(1,1))
      call diagon(dhbb,dnbb,ndb,ndb,eneb,cfrad(1,2))
      dnor=0.d0
      eunc=0.d0
      e2dir=0.d0
      e2exc=0.d0
      do i4=1,ndb
        do i3=1,ndb
          coefb=cfrad(i3,2)*cfrad(i4,2)
          do i2=1,nda
            do i1=1,nda
              coef=cfrad(i1,1)*cfrad(i2,1)*coefb
              dn1=2.d0*(dnaa(i1,i2)*dnbb(i3,i4)
     >           +dnab(i1,i3)*dnab(i2,i4)*sesp)
              dh1=dhaa(i1,i2)*dnbb(i3,i4)+dhbb(i3,i4)*dnaa(i1,i2)
     >           +dhab(i1,i3)*dnab(i2,i4)*2.d0*sesp
              dnor=dnor+coef*dn1
              eunc=eunc+coef*2.d0*dh1
              e2dir=e2dir+2.d0*coef*dh2d(i1,i2,i3,i4)
              e2exc=e2exc+2.d0*coef*dh2e(i1,i2,i3,i4)
            enddo
          enddo
        enddo
      enddo
      eunc=eunc/dnor
      e2dir=e2dir/dnor
      e2exc=e2exc/dnor
      e2exc=sesp/(2.d0*lb+1.d0)*e2exc
      edosc=e2dir+e2exc
      ener=eunc+edosc
      eaut=enea+eneb
      ener=eunc+edosc
      pdir=abs(e2dir)/(abs(e2dir)+abs(e2exc))*100.d0
      pexc=abs(e2exc)/(abs(e2dir)+abs(e2exc))*100.d0
      sumab=0.d0
      do i3=1,ndb
        do i1=1,nda
          sumab=sumab+cfrad(i1,1)*dnab(i1,i3)*cfrad(i3,2)
         enddo
      enddo
      write(*,*)'**********************************************'
      write(*,*)'            Calculo Peturbativo               '
      write(*,*)'**********************************************'
      write(*,*)'-----Exponentes de la funcion de onda------'
      write(*,*)'alfa                              ',alfa
      write(*,*)'alfb                              ',alfb
      write(*,*)'Z alfa efectiva        (alfa)     ',alfa
      write(*,*)'Z alfb efectiva        (alfb*n)   ',alfb*nb
      write(*,*)'-----------Resultados en hartree-----------'
      write(*,*)'energia sin perturbar (h)  valesp ',eunc
      write(*,*)'energia sin perturbar (h)  autoval',eaut
      write(*,*)'correccion perturbativa (h) eprima',edosc
      write(*,*)'componente directa     (h)  Jnl   ',e2dir
      write(*,*)'componente intercambio (h) Knl    ',e2exc
      write(*,*)'energia total       (h)           ',ener
      write(*,*)'--------------Resultados en eV--------------'
      write(*,*)'energia sin perturbar (eV)  ecero  ',eunc*htoev
      write(*,*)'correccion perturbativa (eV) eprima',edosc*htoev
      write(*,*)'energia total       (eV)           ',ener*htoev
      write(*,*)'---------Orbitales monoparticulares----------'
      write(*,*)'Exponentes             ',alfa,alfb
      write(*,*)'Coefientes primero     ',(cfrad(ic,1),ic=1,nda)
      write(*,*)'Coefientes segundo     ',(cfrad(ic,2),ic=1,ndb)
      write(*,*)'ortogonalidad orbitales',sumab
      write(*,*)'-----------------Porcentajes-----------------'
      write(*,*)'100*(eprima/ecero)',100.d0*abs(edosc/eunc)
      write(*,*)'porcentaje directo     ',pdir
      write(*,*)'porcentaje intercambio ',pexc
      write(*,*)
      epert=ener
      return
      end
      subroutine diagon(dene,dnor,nbase,naut,ener,clin)
      implicit real *8 (a-h,o-z)
      parameter (nimax=20)
      parameter (ndim=80)
      dimension dene(nimax,nimax),dnor(nimax,nimax),clin(nimax)
      common /auxdiag/ VENER(NDIM,NDIM),OV(NDIM),EIGENV(NDIM)
     >                ,EVEC(NDIM,NDIM)
      do jc=1,nbase
        vener(jc,jc)=dene(jc,jc)
        ov(jc)=dnor(jc,jc)
        if(jc.gt.1) then
          do ic=1,jc-1
            vener(ic,jc)=dene(ic,jc)
            vener(jc,ic)=dnor(ic,jc)
          enddo
        endif
      enddo
      if(nbase.eq.1) then
        ener=dene(1,1)/dnor(1,1)
        clin(1)=1.d0/sqrt(dnor(1,1))
      else 
        call optimizd(nbase,naut)
        call eigvec(nbase,nbase,naut)
        ener=eigenv(naut)
        sumn=0.d0
        do ic=1,nbase
          sumn=sumn+evec(naut,ic)*dnor(ic,ic)*evec(naut,ic)
          do jc=ic+1,nbase
            sumn=sumn+2.d0*evec(naut,jc)*dnor(jc,ic)*evec(naut,ic)
          enddo
        enddo
        sumn=sqrt(sumn)
        do ic=1,nbase
          clin(ic)=evec(naut,ic)/sumn
        enddo
        sum=0.d0
      endif
      return
      end
      subroutine matrhn(znuc,alfa,alfb,na,nb,la,lb,nda,ndb)
      implicit real *8 (a-h,o-z)
      parameter (nimax=20,limax=10)
      dimension dint1(0:nimax+2*limax)
      common /muncuerpo /dhaa(nimax,nimax),dnaa(nimax,nimax),
     >                   dhbb(nimax,nimax),dnbb(nimax,nimax),
     >                   dhab(nimax,nimax),dnab(nimax,nimax)
c AB
      if(la.eq.lb) then
        asum=alfa+alfb
        lsum=la+lb
        a2=alfb**2
        dosal=2.d0*alfb
        call integ1(asum,dint1)
        do jc=1,ndb
          coef=(jc-1)*(2.d0*lb+jc)
          do ic=1,nda
            ijs=lsum+ic+jc
            dnab(ic,jc)=dint1(ijs)
            dhab(ic,jc)=-0.5d0*(a2*dint1(ijs)-dosal*(lb+jc)*dint1(ijs-1)
     >                         +coef*dint1(ijs-2))-znuc*dint1(ijs-1)
          enddo
        enddo
      endif
c AA
      asum=alfa+alfa
      lsum=la+la
      a2=alfa**2
      dosal=2.d0*alfa
      call integ1(asum,dint1)
      do jc=1,nda
        coef=(jc-1)*(2.d0*la+jc)
        do ic=1,nda
          ijs=lsum+ic+jc
          dnaa(ic,jc)=dint1(ijs)
          dhaa(ic,jc)=-0.5d0*(a2*dint1(ijs)-dosal*(la+jc)*dint1(ijs-1)
     >                      +coef*dint1(ijs-2))-znuc*dint1(ijs-1)
        enddo
      enddo
c BB
      asum=alfb+alfb
      lsum=lb+lb
      a2=alfb**2
      dosal=2.d0*alfb
      call integ1(asum,dint1)
      do jc=1,ndb
        coef=(jc-1)*(2.d0*lb+jc)
        do ic=1,ndb
          ijs=lsum+ic+jc
          dnbb(ic,jc)=dint1(ijs)
          dhbb(ic,jc)=-0.5d0*(a2*dint1(ijs)-dosal*(lb+jc)*dint1(ijs-1)
     >                       +coef*dint1(ijs-2))-znuc*dint1(ijs-1)
        enddo
      enddo
      return
      end
      subroutine matrh2(alfa,alfb,na,nb,la,lb,nda,ndb)
      implicit real *8 (a-h,o-z)
      parameter (nimax=20,limax=10)
      dimension d2dir(nimax+limax,nimax+limax)
      dimension d2exc(nimax+limax,nimax+limax)
      common /mdoscuerpos /dh2d(nimax,nimax,nimax,nimax),
     >                     dh2e(nimax,nimax,nimax,nimax)
      adir=alfa+alfa
      bdir=alfb+alfb
      ldir=0
      l1dir=la+la
      l2dir=lb+lb
      call integ2(adir,bdir,ldir,d2dir)
      aexc=alfa+alfb
      bexc=alfa+alfb
      lexc=lb
      l1exc=la+lb
      l2exc=la+lb
      call integ2(aexc,bexc,lexc,d2exc)
      do i4=1,ndb
        do i3=1,ndb
          do i2=1,nda
            do i1=1,nda
              dh2d(i1,i2,i3,i4)=d2dir(l1dir+i1+i2,l2dir+i3+i4)
              dh2e(i1,i2,i3,i4)=d2exc(l1exc+i1+i3,l2exc+i2+i4)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine integ1(a,dint1)
      implicit real *8 (a-h,o-z)
      parameter (nimax=20,limax=10)
      dimension dint1(0:nimax+2*limax)
      if(a.le.0.d0) then
        write(*,*)'en integ1 problema con el exponente',a
        write(*,*)'deberia ser positivo'
        stop
      endif
      dint1(0)=1.d0/a
      do ic=1,nimax+2*limax
        dint1(ic)=ic/a*dint1(ic-1)
      enddo
      return
      end
      subroutine integ2(a,b,l,dint2)
      implicit real *8 (a-h,o-z)
      parameter (nimax=20,limax=10)
      dimension dint2(nimax+limax,nimax+limax)
      dimension b2(0:nimax,0:nimax),c2(0:nimax,0:nimax)
      dimension dint1(0:2*nimax+2*limax)
      if(a.le.0.d0.or.b.le.0.d0) then
        write(*,*)'en integ2 problema con los exponentes',a,b
        write(*,*)'deberian ser positivos'
        stop
      endif
      if(l.lt.0) then
        write(*,*)'en integ2 problema con l',l
        write(*,*)'deberia ser positivo o cero'
        stop
      endif
      amb=a+b
      dint1(0)=1.d0/amb
      do ic=1,2*nimax+2*limax
        dint1(ic)=ic/amb*dint1(ic-1)
      enddo
      b2(0,0)=1.d0/(a*amb)
      c2(0,0)=1.d0/(b*amb)
      do ic=1,nimax
        b2(ic,0)=(ic*b2(ic-1,0)+dint1(ic))/a
        c2(ic,0)=(ic*c2(ic-1,0)-dint1(ic))/a
      enddo
      do ic=0,nimax
        do jc=1,nimax
          b2(ic,jc)=(jc*b2(ic,jc-1)-dint1(ic+jc))/b
          c2(ic,jc)=(jc*c2(ic,jc-1)+dint1(ic+jc))/b
        enddo
      enddo
      do jc=l+1,nimax
        do ic=l+1,nimax
          dint2(ic,jc)=b2(ic-l-1,jc+l)+c2(ic+l,jc-l-1)
        enddo
      enddo
      return
      end
C
C			RUTINA OPTIMIZ.FOR
C		      GENERAL DIAGONALIZATION
C
	SUBROUTINE OPTIMIZD(NX1,NEIGEN)
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (NDIM=80)
	DIMENSION DELMT(NDIM),UDELMT(NDIM),VENX(NDIM,NDIM),OVX(NDIM),
     >	          EIGX(NDIM)
C
      common /auxdiag/ VENER(NDIM,NDIM),OV(NDIM),EIGENV(NDIM)
     >                ,EVEC(NDIM,NDIM)
C
	DO 10 I1=1,NX1
	OVX(I1)=OV(I1)
	DO 10 I2=1,NX1
	VENX(I1,I2)=VENER(I1,I2)
 10	CONTINUE
C
	CALL CHOLESKYD(veNX,OVX,NDIM,NX1)
C
	CALL HOUSED(VENX,DELMT,UDELMT,NDIM,NX1)
	CALL TRIDIAD(DELMT,UDELMT,NX1,NEIGEN,EIGX)
C
	DO 11 I1=1,neigen
	EIGENV(I1)=EIGX(I1)
 11	CONTINUE
C
	RETURN	
	END
C
C	TRANSFORMA EL PROBLEMA V.P.  H.X=E.O.X EN H'.X=E.X
C	H Y O SON MATRICES SIMETRICAS QUE ENTRAN ASI:
C	H SE GUARDA EN H(I,J) PARA I MENOR O IGUAL A J (TRIANG. SUPERIOR)
C	O SE GUARDA EN H(I,J) PARA I MAYOR QUE J
C	O(I,I) ESTA EN OV(I)
C
C	A LA SALIDA EL TRIANGULO SUPERIOR CONTIENE LA H TRANSFORMADA
C	LA INFORMACION DEL TRIANGULO INFERIOR SE PIERDE
C
C	FASE UNO: DESCOMPOSICION CHOLESKY DEL OveRLAP = L*L(TRANS)
C
	subroutine CHOLESKYD (H,OV,NM,N)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION H(NM,NM),OV(NM)
C
	OV(1)=DSQRT(OV(1))
	DO 140 K=2,N
	T=OV(K)
	DO 115 I=1,K-1
	S=H(K,I)
	IF (I.EQ.1) GOTO 110
	DO 105 J=1,I-1
105	S=S-H(I,J)*H(K,J)
110	H(K,I)=S/OV(I)
115	CONTINUE
	DO 120 J=1,K-1
	Z=H(K,J)
120	T=T-Z*Z
	OV(K)=DSQRT(T)
140	CONTINUE
C
C	INveRSION MATRIZ L
C
	DO 250 K=N,2,-1
	OV(K)=1./OV(K)
	H(K,K-1)=-OV(K)*H(K,K-1)/OV(K-1)
	IF (K.EQ.2) GOTO 250
	DO 240 I=K-2,1,-1
	S=OV(K)*H(K,I)
	DO 230 J=I+1,K-1
230	S=S+H(K,J)*H(J,I)
240	H(K,I)=-S/OV(I)
250	CONTINUE
	OV(1)=1.d0/OV(1)
C
C	MULTIPLICA H POR L INveRSA Y TRANSPUESTA
C
	DO 330 J=N,2,-1
	Z=OV(J)
	DO 320 I=J,1,-1
	S=Z*H(I,J)
	DO 310 K=1,J-1
	Q=H(I,K)
	IF(I.GT.K) Q=H(K,I)
310	S=S+Q*H(J,K)
320	H(I,J)=S
330	CONTINUE
	H(1,1)=OV(1)*H(1,1)
C
C	MULTIPLICA L INveRSA POR H ANTERIOR
C
	DO 430 J=N,2,-1
	DO 420 I=J,1,-1
	S=OV(I)*H(I,J)
	IF (I.EQ.1) GOTO 420
	DO 410 K=1,I-1
	Q=H(K,J)
	IF (K.GT.J) Q=H(J,K)
410	S=S+Q*H(I,K)
420	H(I,J)=S
430	CONTINUE
440	H(1,1)=OV(1)*H(1,1)
	RETURN
	END
C	TRIDIAGONALIZA LA MATRIZ SIMETRICA A
C	SALIDA  ...  DIAGONAL  ...  D (1,... N)
C	        UPPERDIAGONAL  ...  E (1,...,N-1)
C	ALGORITMO DE HOUSEHOLDER
C
C	SOLO USA LA PARTE SUPERIOR DE A
C
	subroutine HOUSED (A,D,E,NM,N)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(NM,NM),D(NM),E(NM)
	DOUBLE PRECISION K2
C
	DO 2330 J=1,N-2
	JP=J+1
	S=0.d0
	DO 2130 I=JP,N
	S1=A(J,I)
2130	S=S+S1*S1
	S1=A(J,JP)
	S=DSQRT(S)
	IF(S1.LT.0.d0) S=-S
	E(JP)=S1+S
	DO 2160 I=J+2,N
2160	E(I)=A(J,I)
	K2=S*(S+S1)
	A(J,JP)=-S
	S1=0.d0
	DO 2260 I=JP,N
	S=0.d0
	DO 2210 K=JP,I
2210	S=S+A(K,I)*E(K)
	IF (I.EQ.N) GOTO 2240
	DO 2230 K=I+1,N
2230	S=S+A(I,K)*E(K)
2240	if(k2.ne.0.d0) S=S/K2
	S1=S1+E(I)*S
2260	D(I)=S
	if(k2.ne.0.d0) S1=-0.5E0*S1/K2
	DO 2280 I=JP,N
2280	D(I)=D(I)+E(I)*S1
	DO 2320 I=JP,N
	S=-E(I)
	S1=-D(I)
	DO 2320 K=I,N
	A(I,K)=A(I,K)+S*D(K)+S1*E(K)
2320	CONTINUE
2330	CONTINUE
	DO 2345 I=1,N-1
	D(I)=A(I,I)
2345	E(I)=A(I,I+1)
	D(N)=A(N,N)
	RETURN
	END
C
C	DETERMINA VALORES PROPIOS MATRIZ TRIDIAGONAL DEFINIDA POR (D,E)
C	
C	W ... CONTIENE LOS VALORES PROPIOS
C	N ... DIMENSION MATRIZ
C	NV .. NUMERO VALORES PROPIOS DESEADOS
C	      ORDENADOS DE MENOR A MAYOR
C
C	OJO, E SE DESTRUYE
C
	subroutine TRIDIAD (D,E,N,NV,W)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(N),E(N),W(NV)
c
	NB=52 
	S=0.d0
	DO 1140 I=1,N
1140	S=S+D(I)
	S=S/N
	DO 1150 I=1,N
1150	D(I)=D(I)-S
	XL1=DABS(D(1))+DABS(E(1))
	XL=DABS(D(N))+DABS(E(N-1))
	IF (XL1.GT.XL) XL=XL1
	DO 1200 I=2,N-1
	XL1=DABS(D(I))+DABS(E(I))+DABS(E(I-1))
	IF(XL1.GT.XL) XL=XL1
1200	CONTINUE
	DO 1210 I=1,N-1
	X=E(I)
1210	E(I)=X*X
	XL1=2.d0*XL/N
	X0=-XL
C	ALGORITMO BISECCION
	DO 1390 I=1,NV
1240	X=X0+XL1
	CALL STURMD (D,E,N,X,NC)
	IF(NC.GE.I) GOTO 1260
	X0=X
	GOTO 1240
1260	IF(NC.EQ.I) GOTO 1310
	X2=X
1280	X=0.5E0*(X0+X2)
	CALL STURMD (D,E,N,X,NC)
	IF(NC.GE.I) GOTO 1260
	X0=X
	GOTO 1280
1310	X2=X
	X1=X
	DO 1360 J=1,NB
	X=0.5E0*(X0+X2)
	CALL STURMD (D,E,N,X,NC)
	IF(I.NE.NC) GOTO 1350
	X2=X
	GOTO 1360
1350	X0=X
1360	CONTINUE
	W(I)=0.5E0*(X0+X2)+S
	X0=X1
1390	CONTINUE
	RETURN
	END
	SUBROUTINE STURMD (D,E,N,X,NC)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(N),E(N)
	P=D(1)-X
	J=1
 	NC=0
	GOTO 1000
900	P=D(J)-X-E(J-1)/P
1000	IF(P.LE.0.) NC=NC+1
	J=J+1
	IF(J.GT.N) RETURN
	IF(P.NE.0.) GOTO 900
	J=J+1
	IF(J.GT.N) RETURN
	P=D(J)-X
	GOTO 1000
	END
C
C			RUTINA EIGVEC.FOR
C
C	=================================================================
C	Rutina para la determinacion de los autovectores de un problema
C	generalizado de autovalores. 
C
C			A.x = Eig B.x
C
C	Como input se requiere las matrices simetricas A, B asi como los
C	autovalores. 
C
C	NEIGV	 : Numero de autovalores.
C	NDIM	 : Dimension fisica de las matrices.
C	N	 : Dimension del problema concreto.
C	XMAT(I,J): Matriz que contiene A en la diagonal y los elementos
C		   no diagonales en los elementos con I<J. Los elementos
C		   no diagonales de B se acomulan en I>J.
C	OV(I)	 : Contiene los elementos diagonales de B.
C	WW(I)	 : Contiene los autovalores.
C	EVEC(I,J): Contiene los autovectores. I es el indice correspondiente
C		   al autovector y J sus componentes.
C	DD(I)	 : A la salida de DLUBKSB contiene el corres. autovector
C	=================================================================
C
	SUBROUTINE EIGVEC(NDIM,N,NEIGV)
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (NMAX=80)
	DIMENSION VPOT(NMAX,NMAX),DD(NMAX),INDEX(NMAX)
        dimension itest(nmax),XMATX(NMAX,NMAX)
C
      common /auxdiag/ XMAT(NMAX,NMAX),OV(NMAX),WW(NMAX)
     >                ,EVEC(NMAX,NMAX)
C
C     Determinamos si el vector N-simo tiene cero toda la fila 
C     excepto diagonal.
C     Ademas vemos si la matriz es diagonal
       IPIV=0
       ntest=0
       do 2001 IEIG=1,N
          itest(ieig)=0
          VTEST=0.d0
          DO 2002 I=1,N
              IF(I.GT.IEIG) VTEST=VTEST+XMAT(IEIG,I)
              IF(I.LT.IEIG) VTEST=VTEST+XMAT(I,IEIG)
 2002     CONTINUE
          IF(VTEST.EQ.0.d0) THEN
             IPIV=IPIV+1
             itest(ieig)=1
             if(ieig.eq.N) ntest=1
          END IF
 2001   CONTINUE
C
C    Caso 1.
C     NTEST.EQ.0: El ultimo elemento no es del tipo especial.
        IF(NTEST.EQ.0) THEN
	DO 4000 IEIG=1,NEIGV
C
	DO 3020 I=1,N-1 
	DO 3020 J=I+1,N
	VPOT(I,J)=XMAT(I,J)-WW(IEIG)*XMAT(J,I)
 3020	VPOT(J,I)=VPOT(I,J)
        iporras=0
	DO 3030 I=1,N
	VPOT(I,I)=XMAT(I,I)-WW(IEIG)*OV(I)
        if(dabs(vpot(i,i)).lt.1.e-14) then
        iporras=1
        jlin1=i
        end if
 3030   continue
C
C		;Como el sistema tiene determinante nulo se elimina la
C		 ultima fila, y la ultima columna juega el papel de 
C		 vector en el sistema C.x=D. Naturalmente la dimension del
C		 sistema disminuye en una unidad.
	DO 3040 I=1,N
 3040	DD(I)=-VPOT(I,N)
C
        if(iporras.eq.0) then
	DETER=1.D0
	CALL DLUDCMP(VPOT,N-1,NMAX,INDEX,DETER)
	CALL DLUBKSB(VPOT,N-1,NMAX,INDEX,DD)
C
	DO 3070 I=1,N-1
 3070	EveC(IEIG,I)=DD(I)
	EveC(IEIG,N)=1.D0
        else
	DO 3071 I=1,N
	EveC(IEIG,I)=0.d0
 3071   CONTINUE
        EveC(IEIG,JLIN1)=1.D0
        end if
C
 4000	CONTINUE
        RETURN
        end if
C
C    Caso 2.
C    IPIV.EQ.N: La matriz es diagonal
        if(ipiv.eq.N) then
	DO 1001 IEIG=1,NEIGV
        energy=WW(IEIG)
        do 1002 kk=1,N
        ediag=xmat(kk,kk)
        EveC(IEIG,kk)=0.D0
        if(energy.eq.ediag) EveC(IEIG,kk)=1.D0
 1002   continue
 1001   continue
        return
        end if
C
C    Caso 3.
C    IPIV.ne.N.and.itest.eq.1: El ultimo elemento es del tipo especial.
        if(ntest.eq.1) then
C     Localizo la linea mas proxima a la ultima que no es del tipo especial
        jlin=0
        do 1101 IEIG=1,N
          if(itest(ieig).eq.0) jlin=ieig
 1101   continue
C     Construyo matriz auxiliar XMATX donde cambio la linea N por la jlin.
        do 1102 kk=1,N
        do 1102 ll=1,N
 1102   XMATX(kk,ll)=XMAT(kk,ll)
        do 1103 kk=1,N
        XMATX(N,kk)=XMAT(JLIN,kk)
        XMATX(JLIN,kk)=XMAT(N,kk)
 1103   CONTINUE
        do 1104 kk=1,N
        XMATAUX=XMATX(kk,JLIN)
        XMATX(kk,JLIN)=XMATX(kk,N)
        XMATX(kk,N)=XMATAUX
 1104 CONTINUE
C
	DO 1201 IEIG=1,NEIGV
C
	DO 1220 I=1,N-1 
	DO 1220 J=I+1,N
	VPOT(I,J)=XMATX(I,J)-WW(IEIG)*XMATX(J,I)
 1220	VPOT(J,I)=VPOT(I,J)
        iporras=0
	DO 1230 I=1,N
	VPOT(I,I)=XMATX(I,I)-WW(IEIG)*OV(I)
        if(dabs(vpot(i,i)).lt.1.e-14) then
        iporras=1
        jlin1=i
        end if
 1230   continue
C
C		;Como el sistema tiene determinante nulo se elimina la
C		 ultima fila, y la ultima columna juega el papel de 
C		 vector en el sistema C.x=D. Naturalmente la dimension del
C		 sistema disminuye en una unidad.
	DO 1240 I=1,N
 1240	DD(I)=-VPOT(I,N)
C
        if(iporras.eq.0) then
	DETER=1.D0
	CALL DLUDCMP(VPOT,N-1,NMAX,INDEX,DETER)
	CALL DLUBKSB(VPOT,N-1,NMAX,INDEX,DD)
C
	DO 1270 I=1,N-1
        IF(I.EQ.JLIN) GOTO 1270
	EveC(IEIG,I)=DD(I)
 1270   CONTINUE
        EveC(IEIG,N)=DD(JLIN)
        EveC(IEIG,JLIN)=1.D0
        else
	DO 1271 I=1,N
	EveC(IEIG,I)=0.d0
 1271   CONTINUE
        EveC(IEIG,JLIN1)=1.D0
        evecaux=EveC(IEIG,N)
        EveC(IEIG,N)=EveC(IEIG,jlin)
        EveC(IEIG,JLIN)=evecaux
        end if
C
 1201	CONTINUE
        return
        end if
C
	RETURN
	END
C
C	==============================================================
C	Given a NxN matrix A, with PHYSICAL dimension NP
C	the routine replaces it by its Lower-Upper LU
C	decomposition, with a row-wise ordering given by
C	INDX (dimensioned at N) and sign exchange D
C	Used with DLUBKSB to solve linear systems or
C	inverting a matrix
C
C	Parameters
C	A	matrix dimensioned at (NP,NP), logical NxN
c	INDX	integer array dimensioned to N at least
C	D	Real*8 number indicating the sign of permutations
C
C	Limitations: N max is 250
C
C	A is destroyed
C	==============================================================
C
	SUBROUTINE DLUDCMP(A,N,NP,INDX,D)
	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX=250,TINY=1.0D-20)
	DIMENSION A(NP,NP),INDX(NP),VV(NMAX)
	D=1.D0
	DO 12 I=1,N
	AAMAX=0.D0
	DO 11 J=1,N
	IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11	CONTINUE
	IF (AAMAX.EQ.0.D0) PAUSE 'Singular matrix.'
	VV(I)=1.D0/AAMAX
12	CONTINUE
	DO 19 J=1,N
	IF (J.GT.1) THEN
	DO 14 I=1,J-1
	SUM=A(I,J)
	IF (I.GT.1)THEN
	DO 13 K=1,I-1
	SUM=SUM-A(I,K)*A(K,J)
13	CONTINUE
	A(I,J)=SUM
	ENDIF
14	CONTINUE
	ENDIF
	AAMAX=0.D0
	DO 16 I=J,N
	SUM=A(I,J)
	IF (J.GT.1)THEN
	DO 15 K=1,J-1
	SUM=SUM-A(I,K)*A(K,J)
15	CONTINUE
	A(I,J)=SUM
	ENDIF
	DUM=VV(I)*ABS(SUM)
	IF (DUM.GE.AAMAX) THEN
	IMAX=I
	AAMAX=DUM
	ENDIF
16	CONTINUE
	IF (J.NE.IMAX)THEN
	DO 17 K=1,N
	DUM=A(IMAX,K)
	A(IMAX,K)=A(J,K)
	A(J,K)=DUM
17	CONTINUE
	D=-D
	VV(IMAX)=VV(J)
	ENDIF
	INDX(J)=IMAX
	IF(J.NE.N)THEN
	IF(A(J,J).EQ.0.D0) A(J,J)=TINY
	DUM=1.D0/A(J,J)
	DO 18 I=J+1,N
	A(I,J)=A(I,J)*DUM
18	CONTINUE
	ENDIF
19	CONTINUE
	IF(A(N,N).EQ.0.D0)A(N,N)=TINY
	RETURN
	END
C
C	==============================================================
C	Routine to solve the linear system
C	A. X = B
C	where A is the matrix already in LU decomposition form
C	(see DLUDCMP routine)
C	and B is the rhs vector.
C
C	On output, B contains the solution
C
C	The way of solving a system of linear equations is
C	Call DLUDCMP to change A into its LU decomposition
C	call afterwards DLUBCKS (backsubstitution)
C
C	A	dimensioned at (NP,NP), logical NxN
C	INDX	order of rows (see DLUDCMP), dimensioned to N
C	B	Column vector, dimensioned (at least) to N
C
C	A is not destroyed, so this routine may be called
C	many times, per one call of DLUDCMP
C	==============================================================
C
	SUBROUTINE DLUBKSB(A,N,NP,INDX,B)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(NP,NP),INDX(NP),B(NP)
	II=0
	DO 12 I=1,N
	LL=INDX(I)
	SUM=B(LL)
	B(LL)=B(I)
	IF (II.NE.0)THEN
	DO 11 J=II,I-1
	SUM=SUM-A(I,J)*B(J)
11	CONTINUE
	ELSE IF (SUM.NE.0.D0) THEN
	II=I
	ENDIF
	B(I)=SUM
12	CONTINUE
	DO 14 I=N,1,-1
	SUM=B(I)
	IF(I.LT.N)THEN
	DO 13 J=I+1,N
	SUM=SUM-A(I,J)*B(J)
13	CONTINUE
	ENDIF
	B(I)=SUM/A(I,I)
14	CONTINUE
	RETURN
	END
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
