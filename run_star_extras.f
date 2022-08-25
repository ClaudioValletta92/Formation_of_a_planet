! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
!
!   this file is  part of mesa.

!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! **********************************************************************



      module run_star_extras



      use star_lib
      use star_def
      use const_def
      use crlibm_lib

      implicit none

      integer :: time0, time1, clock_rate
      real(dp), parameter :: expected_runtime = 1 ! minutes
      integer :: numTacc,j                          !track location of accretion times read





      contains





      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
	 integer :: i
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_start_step => extras_start_step
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% other_opacity_factor => opacity_factor_routine
         s% other_energy => energy_routine
         s% other_kap_get_Type1 => kap_freedman_grains

         !s% other_eosDT_get => my_eosDT_ideal_gas_get


           !s% min_timestep = 10**(s% x_ctrl(7))*secyer
           !s% dt_next = 10**(s% x_ctrl(7))*secyer   !set as thermal timescale  s% cgrav(1)*s% mstar/s% r(1)*s% L(1)]


	       s% other_energy => energy_routine
         s% other_wind => planetary_accretion
      end subroutine extras_controls

            subroutine energy_routine(id, ierr)
               type (star_info), pointer :: s
               !use const_def, only: Rsun
               integer, intent(in) :: id
               integer, intent(out) :: ierr


               real(dp) initialposition, initialvelocity
               real(dp) size,sunr,deltam,beforesize
                real(dp) temperature, rhoair,onemass
                real(dp) remainingmass,beforemetallicity,beforemetallicity2
                real(dp) accretedmass, factor,dummy1,diff_fraction
                integer i, Layer_number, j,m
      	 !use star_lib, only: star_ptr
      	 double precision :: extraheat,junk,diag_heat
               integer :: k,n,counter,z,p,numOfAcc,zeta,jafter,jbefore,indexI
      	 double precision :: core_epsrad,Rpl,pressureDep,pressureHere,random_dp
      	 real(dp) :: tauHere,Vesc, KE,massTot,Cd,phi,Mpl,Rhopl,H,g,mH,Tacc,cv,rockfr
      	 !real(dp), DIMENSION(30000000) :: readR,readM,readT
               !REAL(dp), DIMENSION(:), ALLOCATABLE :: readR,readM,readT
               real(dp), DIMENSION(700) :: arrayKE
               real(dp), DIMENSION(5000) :: arrayI
               real(8) e0,ef,cp,tfus,fraction
               cp = 4.19*10**7
               e0 = 2.8D10
               ef = 3.30D9
               tfus = 373.



               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               ierr = 0


               sunr=6.59700*10**10
               ierr = 0

               m= 1

               size = s%xtra1
               onemass = 4./3.*pi*s%xtra17*size**3


               accretedmass = s% dt/secyer*5.97d27/(10**6)
               factor = accretedmass / onemass

               initialvelocity = 100000
               initialposition = 2*sunr*10**s% log_surface_radius
               fraction = 0




               do while(m .LE. s% nz)
                 deltam=0
                beforesize = size

                if(m .eq. s% nz)THEN
                  size = size +Integrand(s% r(m),s%mstar,initialvelocity,initialposition,exp(s%lnT(m)),exp(s%lnd(m)))&
                  &*(s%R_center-s%r(m))
                else
                  size = size +Integrand(s% r(m),s%mstar,initialvelocity,initialposition,exp(s%lnT(m)),exp(s%lnd(m)))&
                  &*(s%r(m+1)-s%r(m))
                ENDIF

                deltam=4./3.*pi*beforesize**3*s%xtra17-4./3.*pi*size**3*s%xtra17

                       fraction = fraction+deltam
                       diff_fraction = deltam/onemass*s% dt/secyer*s%xtra5


                s%xtra1_array(m) = deltam

                if(m .eq. s% nz)THEN
                  s%xtra2_array(m) = -Drag(velocity(s% r(m),s%mstar,initialvelocity,initialposition)&
                  &,size,exp(s%lnd(m)),exp(s%lnT(m)),exp(s%lnP(m)))*(s%R_center-s%r(m))+&
                  &+0.5*deltam*velocity(s% r(m),s%mstar,initialvelocity,initialposition)**2-&
                 &-(cp*tfus+e0)*deltam
                else
                  s%xtra2_array(m) =-Drag(velocity(s% r(m),s%mstar,initialvelocity,initialposition)&
                  &,size,exp(s%lnd(m)),exp(s%lnT(m)),exp(s%lnP(m)))*(s%r(m+1)-s%r(m))+&
                  &+0.5*deltam*velocity(s% r(m),s%mstar,initialvelocity,initialposition)**2-&
                  &-(cp*tfus+e0)*deltam
                ENDIF
                if(s%xtra2_array(m) .ge. 0 )then
                        s% extra_heat(m)=s%xtra2_array(m)/(s%dm(m)*s%dt)*diff_fraction
        else
                s%extra_heat(m)=0
        endif
             !   write(*,*) 'extra heat',s% extra_heat(m)
                m = m+1
        enddo

                 s%xtra3 = fraction/onemass

               s% extra_heat = 0.d0

               cv = 1d10
               k = s% nz
               extraheat = -cv * s% M_center * s% T(s% nz) * s% dlnT_dt(s% nz) / s% dm(s% nz) ! erg/g/sec
                 !assuming dlnT_dt is in s^-1

               ! E XTRA HEATING CONSTANT IN SPACE AND TIME
               ! Heat produced by radioactive decay due to 40K, 235U, 238U, 232Th, respectively
               rockfr = 1.0
               k = s% nz
               core_epsrad = 36.7d-8 * exp(-5.543d-10 * s% star_age) ! erg/g/sec   40K
               core_epsrad = core_epsrad + 4.6d-8 * exp(-9.8485d-10 * s% star_age) ! erg/g/sec  235U
               core_epsrad = core_epsrad + 2.3d-8 * exp( -1.5513d-10 * s% star_age)! erg/g/sec  238U
               core_epsrad = core_epsrad + 1.3d-8 * exp( -0.4948d-10 * s% star_age)! erg/g/sec  232Th

               s% extra_heat(k) =s% extra_heat(k)+ (extraheat+ rockfr * s% M_center * core_epsrad / s% dm(k)) ! erg/g/sec, core heat flux density


               write(*,*) 'extra heat',s% extra_heat(k)

            end subroutine energy_routine


                            function Drag(velocity,size,rho,temperature,P) result(ValueofDrag)
                              real(dp)  ValueofDrag,fpi,dkn,psi
                              real(dp) velocity,size,rho,temperature,P




                                !if(Reynold_Number(velocity,size,rho,temperature) .LE. 1) then
                                  fpi=1.019*10**9*rho
                                	dkn=1./(size*fpi);
                                	psi=1.+dkn*(1.25+.42*exp(-.87/dkn));
                                	ValueofDrag=18.85*Viscosity(temperature)*velocity*size/psi
                                !else if (Reynold_Number(velocity,size,rho,temperature) .GT. 1) then
                                  ValueofDrag=Dk(velocity,size,rho,temperature,P)*size*size*PI*rho*velocity*velocity
                                  !return



                             !ENDIF


                           end function Drag


                                              function Dk(velocity,size,rho,temperature,P) result(ValueofDk)
                                                real(dp)  valueofdk

                                                real(dp)  velocity,size,rho,temperature,P


                                                      if(Reynold_Number(velocity,size,rho,temperature) .LE. 1) then
                                                   valueofdk = 1
                                                   		return

                                                    else if(Reynold_Number(velocity,size,rho,temperature) .LE. 10**3) then
                                                      if(Mach_Number(P,rho,velocity) .LE. 1) then
                                                   		valueofdk = 6/sqrt(Reynold_Number(velocity,size,rho,temperature))
                                                      return
                                                    else if(Mach_Number(P,rho,velocity) .GT. 1) then
                                                   			valueofdk = 1.1-log(Reynold_Number(velocity,size,rho,temperature))/6.;
                                                        return
                                                      endif

                                                    else if(Reynold_Number(velocity,size,rho,temperature) .LE. 10**5) then

                                                   		valueofdk=0.2
                                                   		if(Mach_Number(P,rho,velocity)>1) then

                                                   			ValueofDk=0.5
                                                   		endif

                                                    else	if(Reynold_Number(velocity,size,rho,temperature) .GT. 10**5) then

                                                   		valueofdk=0.15
                                                  		if(Mach_Number(P,rho,velocity)>1) then

                                                   			valueofdk=0.5;
                                                   		  return
                                                      endif
                                                    ENDIF



                                              end function dk

                                              function Mach_Number(P,rho,velocity) result(Mach)
                                                real(dp)  Mach,P,rho,velocity

                                                Mach = velocity/SQRT(1.4*P/rho)

                                              end function Mach_Number

                                              function Reynold_Number(velocity,size,rho,temperature) result(Reynold)
                                                real(dp)  Reynold,velocity,rho,size,temperature

                                                 Reynold = velocity*rho*size/Viscosity(temperature)

                                              end function Reynold_Number
                                              function Viscosity(temperature) result(ValueOfViscosity)
                                                real(dp)  ValueOfViscosity,temperature


                                                ValueOfViscosity = 4.3*10**(-6)*SQRT(temperature)

                                              end function Viscosity

      !#########module to test extra heating sensitivity######
         SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)

            DEALLOCATE(seed)
          END SUBROUTINE
     !!!!!##############################################




      subroutine planetary_accretion(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
         !use star_defe
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf

          real core_avg_rho! surface values (cgs)
         !double precision              :: DM !The change of the mass of the core per timestep,Msun units
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         real(8) oldradius,newradius,maxabsdot,middleradius,sunr,soundspeed,Bondi,righthandside

         real(8) AccretedLissauer,EarthRadius,orbitaldistance,PeriodinYS
         real(8) ds, Hillradius, Fg,mdotz,sigmaz,sigmazinit,deltasigmaz
         real(8) captradius,initial_M_center,accretedmass
         real(8) Planetesimalradius,Planetesimaldensity,Japanese
         real(8) initial_accretion_radius1,accretion_radius1,initial_available_mass
         real(8) initial_accretion_radius2,accretion_radius2
         real(8) available_mass,Derivative_RoverCore,Derivative_RoverXY
         real(8) normalizingvalue,orbitaldistanceinau
         real(8) alphavisc,sigma0,tdep,hp,di,kappa,sigmaacc,sigmaun
         real(8) smallm,smallh,smallt,t1,t2,gapfactor
         integer Layer_number, j, initialmodelnumber,m


         ierr = 0
         call get_star_ptr(id, s, ierr)


         s% num_accretion_species = 3



        !---------------------------------------------------!
         !xtra10 to check if i am in the runaway phase or not!
         if (ierr /= 0) return
         s% xtra10 = 0


         if(s% xtra10 ==1)then
           s% xtra10 =1
         else
         s% xtra10 = 0
       ENDIF
       orbitaldistance = s%xtra4
       orbitaldistanceinau = orbitaldistance / (1.496*10**13)
       periodinYS = s%xtra9

       !---------------------------------------------------!
       !xtra2 and xtra7 for the Hilld and Bondi radius
       Hillradius = orbitaldistance*(s% mstar / (3*Msun))**(1./3.)
       soundspeed = 0.4*(orbitaldistanceinau / 10.)**(3./14.)*10**5
       Bondi = 6.67/(10**8)*s%mstar/(2*soundspeed**2+4*6.67/(10**8)*s%mstar/Hillradius)
       s% xtra2 = Bondi
       s%xtra7 = Hillradius
       !---------------------------------------------------!


        if(s%xtra16 .ge. 1)then
                s%xtra16 = 1
        endif

        sunr=6.59700*10**10
        newradius=sunr*10**s% log_surface_radius



        s%xtra16 = 0
        s%xtra3 = 0

        w = 0



        call SolidAccretionRate(id)


        s% M_center =s% M_center + s% dt/secyer * s%xtra5 * (1-s%xtra16)

        s%xtra25 = s%xtra25 + s% dt/secyer*s%xtra5

        call PlanetaryRadius(id)
        s% L_center=3./5.*6.67d-8*s%M_center/s%R_center*s%xtra5*3.17d-8*(1-s%xtra3)
      if(s%xtra3 .ge. 0.99)then
               s% L_center = 1d18
        endif


        s% mstar = s% M_center + s% xmstar

        !-------------------------------------------------------!

         Hillradius = orbitaldistance*(s% mstar / (3*Msun))**(1./3.)
         Bondi = 6.67/(10**8)*s%mstar/(2*soundspeed**2+4*6.67/(10**8)*s%mstar/Hillradius)
         s% xtra2 = Bondi
         s%xtra7 = Hillradius


          s% mass_change =s%xtra21*4./3.*pi*(-newradius**3+s% xtra2**3)/(Msun*s%dt/secyer)

         if(s% mass_change .ge.1d-5)then
                 s%mass_change = 0
         endif

       if(s% mass_change .le. 0)then
           s% mass_change = 0
        ENDIF



     if(s% mstar .ge. 100*5.97d+27)then

        call GasAccretionRate(id)
        s% xtra10 = 1
         s% mass_change =5.97d27*(1-s%mstar/(318*5.97d27))*s%xtra11/Msun
        !s% mass_change = AccretedLissauer/msun
       ENDIF


       write(*,*) 'Planetary accretion sub---Model number is t1', s%model_number
       write(*,*) 'Planetary accretion sub---Star age is', s% star_age/1000000
       write(*,*) 'Planetary accretion sub---Mass of the envelope  in plantary accretion', s% xmstar/(5.97*10**27)
       write(*,*) 'Planetary accretion sub---Mass of the star  in plantary accretion', s% mstar/(5.97*10**27)
       write(*,*) 'Planetary accretion sub---Mass of the core  in plantary accretion', s% M_center/(5.97*10**27)



       s%mass_change = s%mass_change+s%xtra5*(s%xtra16)/Msun



end subroutine planetary_accretion


      subroutine SolidAccretionRate(id)
        integer, intent(in) :: id
        integer :: ierr






        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return



       !---------------------------------------------------!
       !Masakazu Shiraishi and Shigeru Ida
       !call ShiraishiIda(id)
       !---------------------------------------------------!



       !Pollack 1996
       ! ------------------------------------------------------!
       call Pollack96Solidaccretionrate(id)
       ! ------------------------------------------------------!

       !constant solid accretion rate
       ! ------------------------------------------------------!
       !call Constantsolidaccretionrate(id)
       ! ------------------------------------------------------!

       !Rafikov
       ! ------------------------------------------------------!
       !call Rafikovsolidaccretionrate(id)
       !-------------------------------------------------------!


       !Oligarchic planetesimal accretion and giant planet formation!
       !---------------------------------------------------!
       !call Fortieretal2007(id)
       !---------------------------------------------------!

       !---------------------------------------------------!




      end subroutine SolidAccretionRate


      subroutine solidsurfacedensity(id)
        integer, intent(in) :: id
        integer :: ierr

        real(8) mdotz,orbitaldistance,sigmazinit,initial_M_center,captradius
        real(8) initial_accretion_radius1,initial_accretion_radius2,fg,alpha
        real(8) sigmaz,PeriodinYS,addedmass1,addedmass2,totaladded,totalmass
        real(8) sigma0,sigmagas,sigmasolid,sigmainit,orbitaldistanceinau
        real(8) accretion_radius1,accretion_radius2,initial_available_mass,available_mass





        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        au = 1.496*10**13

        orbitaldistance = s%xtra4

        alpha = 1.5
        sigma0 = 9.5

        orbitaldistanceinau = orbitaldistance / au

        initial_M_center = 5.9700000889599331D+26
        captradius = s%xtra6
        periodinYS = s%xtra9


        sigmainit = sigma0*(orbitaldistanceinau / 5.2)**(-alpha)



        initial_accretion_radius1 = orbitaldistance+SQRT(12.)*orbitaldistance*(initial_M_center / (3*Msun))**(1./3.)
        initial_accretion_radius2 = orbitaldistance-SQRT(12.)*orbitaldistance*(initial_M_center / (3*Msun))**(1./3.)

        s%xtra12 = orbitaldistance+SQRT(12.)*s%xtra7
        s%xtra14 = orbitaldistance-SQRT(12.)*s%xtra7
        addedmass1 = sigmainit*pi*(s%xtra12**2-s%xtra12_old**2)
        addedmass2 = abs(sigmainit*pi*(s%xtra14**2-s%xtra14_old**2))

        totaladded = addedmass1+addedmass2



        if(s% model_number .le. 290)then
          s%xtra15 =sigmainit*pi*(initial_accretion_radius1**2-initial_accretion_radius2**2)
        ENDIF

        totalmass = s%xtra15 + totaladded - s%xtra5*s% dt/secyer

        sigmasolid = totalmass/(pi*(s%xtra12**2-s%xtra14**2))

        s%xtra15 = totalmass



        !---------------------------------------------------!


        

        if(sigmasolid .le. 0) then
          sigmasolid =  s%xtra13
        endif

        write(*,*) 'solidsurfacedensity sub---rho is', sigmasolid

        !---------------------------------------------------!
        s%xtra13 = sigmasolid


      end subroutine solidsurfacedensity



      subroutine Rafikovsolidaccretionrate(id)
              integer, intent(in) :: id
              integer :: ierr

              real(8) mdotz,sigmasolid


              !Rafikov solid accretion rate, the highest
              !Formula taken from Helled & Venturini, minineptunes

              type (star_info), pointer :: s
              ierr = 0
              call star_ptr(id, s, ierr)
              if (ierr /= 0) return

              call solidsurfacedensity(id)
              sigmasolid = s%xtra13

              mdotz = 6.47*2*pi/s%xtra9*SQRT(s%xtra6/s%xtra7)*sigmasolid*s%xtra7*s%xtra7
              s%xtra5 = mdotz
              write(*,*) 'Rafikov sub---mdotz is', mdotz

      end subroutine Rafikovsolidaccretionrate


      subroutine Constantsolidaccretionrate(id)
              integer, intent(in) :: id
              integer :: ierr

              real(8) mdotz

              type (star_info), pointer :: s
              ierr = 0
              call star_ptr(id, s, ierr)
              if (ierr /= 0) return

              s%xtra5 = 5.97d27/(10**6)
              write(*,*) 'constant sub---mdotz is', s%xtra5

      end subroutine Constantsolidaccretionrate



      subroutine Pollack96Solidaccretionrate(id)
              integer, intent(in) :: id
              integer :: ierr

              real(8) mdotz,sigmasolid,fg,captradius

              !Standard solid accretion rate, as taken from Pollack 1996
              type (star_info), pointer :: s
              ierr = 0
              call star_ptr(id, s, ierr)
              if (ierr /= 0) return

              call solidsurfacedensity(id)
              sigmasolid = s%xtra13
              call fg2(id)
              fg = s%xtra8
              captradius = s%xtra6

              mdotz = 2*pi*pi*fg*sigmasolid/(s%xtra9)*captradius**2
              !if(s%M_center .le. 9*5.97d27)THEN
               ! mdotz = 6.47*2*pi/s%xtra9*SQRT(s%xtra6/s%xtra7)*sigmasolid*s%xtra7*s%xtra7
              !ENDIF
              s%xtra5 = mdotz
              write(*,*) 'P96 sub---mdotz is', mdotz

      end subroutine Pollack96Solidaccretionrate



      subroutine Fortieretal2007(id)
              integer, intent(in) :: id
              integer :: ierr

              real(8) mdotz,sigmasolid,fg,orbitaldistance,periodinYS
              real(8) Planetesimalradius, Planetesimalmass,beta,cd,diskdensity
              real(8) eccentricity,inclination,planscaleheight,vrel,planetesimaldensity

              type (star_info), pointer :: s
              ierr = 0
              call star_ptr(id, s, ierr)
              if (ierr /= 0) return



              !Oligcachic growth of the planet as implemeneted by
              !Forteir et al 2007,
              !Fortier, A., Benvenuto, O. G. & Brunini, A. 2007, A&A, 473, 311.
              orbitaldistance = s%xtra4
              call solidsurfacedensity(id)
              sigmasolid = s%xtra13
              periodinYS = s%xtra9


              Planetesimalradius = s%xtra1
              planetesimaldensity = s%xtra17
              planetesimalmass = 4./3.*pi*planetesimaldensity*Planetesimalradius**3
              beta = 10
              cd = 1
              diskdensity = 5.2/(10**11)

              eccentricity = 1.7*planetesimalmass**(1./15.)*s%mstar**(1./3.)*planetesimaldensity**(2./15.)*&
              &1/(beta**(1./5.)*diskdensity**(1./5.)*Msun**(1./3.)*orbitaldistance**(1./5.)*cd**(1./5.))

              inclination = eccentricity / 2.
              planscaleheight = inclination*orbitaldistance
              vrel =orbitaldistance*2*pi/periodinYS*sqrt(eccentricity*eccentricity + inclination*inclination)

              mdotz = 3*sigmasolid*pi/(2*planscaleheight)*s%xtra6**2*vrel

              s%xtra5 = mdotz
              write(*,*) 'Fortier sub---mdotz is', mdotz

      end subroutine Fortieretal2007


      subroutine ShiraishiIda(id)
          integer, intent(in) :: id
          integer :: ierr
          real(8) mdotz,planetradius,tdamp,vscat
          real(8) beta,beta2,alpha,alpha2,fd,mdotey,tgacc,epsilon,orbitaldistance
          real(8) eta,h,vh,scat,mdotlargeeta,mdotsmalleta,planetaradius,kilometer,periodinYS
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

                    !Oligcachic growth of the planet as implemeneted by
                    !Masakazu Shiraishi and Shigeru Ida
                    !https://arxiv.org/pdf/0805.2200.pdf

          planetradius = (10**(s% log_surface_radius)*Rsun)
          kilometer = 100000
          beta  = -6.
          beta2  = -6.
          tdamp = 1d5
          alpha2 = 1.4
          alpha  = 0.8
          fd = 0.7
          orbitaldistance = s%xtra4
          periodinYS = s%xtra9

          mdotey = (s% mstar-s% mstar_older)/(s%dt/secyer)

          tgacc = s%mstar / mdotey
          epsilon = 4.1*(orbitaldistance/(5*1.496*10**13))**(3./2.)/((s%mstar/(5.97d72))**(-1./3.)*(tgacc/10000))
          eta = 0.8*(orbitaldistance/(5*1.496*10**13))**(3./4.)*(tdamp/10000)**(1./2.)/((tgacc/10000)*(s%mstar/(5.97d27))**(1./6.))
          h = (s%mstar/(3*MSun))**(1./3.)
          vh = 8/tgacc
          vscat = 0.22*h/periodinYS
          mdotlargeeta = 10**beta*(planetradius/(6378*kilometer))**2*fd*epsilon**alpha
          mdotsmalleta = 10**beta2*(planetradius/(6378*kilometer))**2*fd*eta**alpha2

         if(eta .gt. 1)then
           mdotz = mdotlargeeta*5.97d27
         endif
         if(eta .le. 1)THEN
           mdotz = mdotsmalleta*5.97d27
         endif

         if(s%model_number .le. 400)then
           mdotz = 5.97d27/(1d6)
         ENDIF
         if(mdotz .eq. 0 )then
           mdotz = 5.97d19
        ENDIF
         s%xtra5 = mdotz
         write(*,*) 'Shiraishi sub---mdotz is', mdotz
         write(*,*) 'Shiraishi sub---ey is',(s% mstar-s% mstar_older)/(s%dt/secyer)
         write(*,*) 'Shiraishi sub---eta is',eta
         write(*,*) 'Shiraishi sub---tdamp is',tdamp
         write(*,*) 'Shiraishi sub---epsilon is',epsilon

     end subroutine ShiraishiIda


            subroutine GasAccretionRate(id)
              integer, intent(in) :: id
              integer :: ierr


              type (star_info), pointer :: s
              ierr = 0
              call star_ptr(id, s, ierr)
              if (ierr /= 0) return

              !gas accretion rate


              !----------------------------------------------------------------------------!
              !Lissauer


              !----------------------------------------------------------------------------!
               call Lissauer(id)
              !----------------------------------------------------------------------------!

              !Tanagawa
              !----------------------------------------------------------------------------!
                !call tanagawa(id)
              !----------------------------------------------------------------------------!


              !ginzbug chiange
              !----------------------------------------------------------------------------!
                !call GinzburgChiang(id)
              !----------------------------------------------------------------------------!
            end subroutine GasAccretionRate


           subroutine GinzburgChiang(id)
                          integer, intent(in) :: id
                          integer :: ierr

                          real(8) smallm,smallh,smallt,t1,t2,t0,m1,m2
                          real(8) gapfactor,tdim,accretedlissauer,alphavisc,hp
                          real(8) newtonconstant,periodinYS,orbitaldistance
                          type (star_info), pointer :: s
                          ierr = 0
                          call star_ptr(id, s, ierr)
                          if (ierr /= 0) return



                          newtonconstant = 6.67259/(10**8)
                          PeriodinYS = s%xtra9
                          orbitaldistance = s%xtra4

                          !Gas accretion rate as implemented by GinzburgChiang
                              alphavisc = 1./(10**4)
                                smallm = s%mstar/(2d33)
                                smallh = hp/orbitaldistance
                                smallt = 2*pi*s%star_age/PeriodinYS
                                t1 = (smallh**2/alphavisc**7)**(1./5.)
                                t2 = (smallh**4/alphavisc**9)**(1./5.)
                                t0 = pi/(PeriodinYS**2)*1/(secyer*secyer*newtonconstant*5.2/10**11)
                                m1 = smallh**(58./21.)/(alphavisc**(2./21.)*t0**(1./3.))
                                m2 = smallh**(14./5.)/(alphavisc**(2./15.)*t0**(1./3.))
                                write(*,*) 't0',t0
                                write(*,*) 't1',t1
                                write(*,*) 't2',t2
                                write(*,*) 'm',smallm
                                write(*,*) 'm1',m1
                                write(*,*) 'm2',m2


                                if(t0 .ge. t2)then
                                  gapfactor = smallm**4*smallh**(-53./5.)*alphavisc**(-7./5.)
                                endif
                                if(t0.le. t2)then
                                  gapfactor = smallm**4*smallh**(-391./35.)*smallt**(5./7.)*alphavisc**(-4./35.)
                                endif
                                if(t0 .le. t1)then
                                  gapfactor = smallm**4*smallh**(-549./49.)*smallt**(39./49.)
                                endif


                                if(smallm .ge. m2)then
                                   tdim = periodinYS/(2*pi)*t0*smallm**3*smallh**(-38./5.)*alphavisc**(-7./5.)
                                endif
                                if(smallm .le. m2)then
                                  tdim = periodinYS/(2*pi)*t0**(7./2.)*smallm**(21./2.)*smallh**(-143./5.)*alphavisc**(-2./5.)
                                endif
                                if(smallm .le. m1)then
                                  tdim = periodinYS/(2*pi)*t0**(49./10.)*smallm**(147./10.)*smallh**(-201./5.)
                                endif

                                AccretedLissauer = s%mstar/tdim
                                  write(*,*) 'GinzburgChiang sub - gas accretion rate is ', AccretedLissauer

                              s%xtra11 = AccretedLissauer


          end subroutine GinzburgChiang

          subroutine tanagawa(id)
                          integer, intent(in) :: id
                          integer :: ierr

                          real(8) alphavisc,sigma0,tdep,hp,di,kappa,kappaprime,sigmaun
                          real(8) sigmaacc,gapdensity,deltar1,deltar2,accretedlissauer,orbitaldistance,periodinYS



                          type (star_info), pointer :: s
                          ierr = 0
                          call star_ptr(id, s, ierr)
                          if (ierr /= 0) return

                          orbitaldistance = s%xtra4
                          periodinYS = s%xtra9
                          !tanagawa
                          alphavisc = 1./(10**4)
                          sigma0 = 700
                          tdep = 3.
                          hp = (orbitaldistance*5.2**0.25)/(10**1.5)
                          di = (s%mstar/(2d33))**(4./3.)*0.29*(orbitaldistance/hp)**2*orbitaldistance**2*2*pi/PeriodinYS
                          kappa = (orbitaldistance/hp)**5*(s%mstar/(2d33))**2/alphavisc
                          kappaprime = (orbitaldistance/hp)**3*(s%mstar/(2d33))**2/alphavisc
                          sigmaun = sigma0*exp(-s%star_age/(1000000*tdep))
                          sigmaacc = sigmaun/(1+0.034*kappa)
                          gapdensity = 1/(1+0.034*kappa)
                          deltar1 = (0.25/(1+0.034*kappa)+0.08)*kappaprime**(0.25)*s%xtra4
                          deltar2 = 0.33*kappaprime**(0.25)*s%xtra4

                          AccretedLissauer=di*sigmaacc
                          s%xtra11 = AccretedLissauer
                          write(*,*) 'Tanagawa sub - gas accretion rate is ', AccretedLissauer
                          !----------------------------------------------------------------------------!

          end subroutine tanagawa


subroutine Lissauer(id)
    integer, intent(in) :: id
    integer :: ierr
    real(8) c0,c1,c2,righthandside,diskdensity,AccretedLissauer
    real(8) sigmaacc,gapdensity,deltar1,deltar2
    real(8) orbitaldistance,periodinYS


    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return


    c0 = -18.67
    c1 = -8.97
    c2 = -1.23
    orbitaldistance = s%xtra4
    periodinYS = s%xtra9

    righthandside = c0 + c1 * LOG10(s%mstar/Msun) +c2 * (LOG10(s%mstar/Msun))**2

    Diskdensity = 700*((3.-s%star_age/1000000))/3.
    Diskdensity = 700


   AccretedLissauer=Diskdensity*orbitaldistance**2/(PeriodinYS)*10**(righthandside)/5.97d+2
   write(*,*) 'Lissauer sub - gas accretion rate is ', AccretedLissauer
   s%xtra11 = AccretedLissauer


end subroutine Lissauer


subroutine PlanetaryRadius(id)
  integer, intent(in) :: id
  integer :: ierr
  real :: con1
  real(8) core_avg_rho,rho0,c,n,r1,GravConstant,m1
  real k1,k2,k3
  real(8) EarthRadius,radius4,centerm
  real(8) c1,c2,c3,c4,c5,c6,c7,c8

  type (star_info), pointer :: s
  ierr = 0
  call star_ptr(id, s, ierr)
  if (ierr /= 0) return
  EarthRadius = 6.3781*10**8
  GravConstant = 6.67259/(10**(11))
  n = 0.528
  c = 0.00349
  rho0 = 8300
  k1 = -0.20945
  k2 = 0.0804
  k3 = 0.394
  r1 = GravConstant**(-0.5)*rho0**(1./(2.*n)-1)*c**(-1./(2*n))
  m1 = rho0*r1**3

  c1 = 9.68787696e-01
  c2 =3.01194386e-01
  c3 =-2.50831307e-02
  c4 =1.16031438e-03
  c5 = -2.88618946e-05
  c6 =3.88103876e-07
  c7 =-2.65821050e-09
  c8 = 7.26617658e-12
  centerm = s%M_center/(5.97d27)
  core_avg_rho = 3.2
  radius4 = c1+c2*centerm+c3*centerm**2+c4*centerm**3+&
          & c5*centerm**4+c6*centerm**5+c7*centerm**6+c8*centerm**7
  s% R_center = pow_cr(s% M_center/(core_avg_rho*4.*pi/3.),1.d0/3.d0)
  !PROBABILISTIC MASS-RADIUS RELATIONSHIP FOR SUB-NEPTUNE-SIZED PLANETS
if(s% M_center .ge. 1*5.97d27)then
 ! s% R_center = EarthRadius/(2.7**(1./1.3))*(s% M_center/5.97d27)**(1./1.3)
!     s% R_center = 100*r1*10**(k1 + 1./3.*LOG10(s% M_center/1000/m1)-k2*(s%M_center/1000/m1)**k3) done
  !s% R_center = 10**(k1 + 1./3.*LOG10(s% M_center/1000/m1) -k2*(s%M_center/1000/m1)**k3)
s%R_center = radius4*EarthRadius

endif



end subroutine PlanetaryRadius


      subroutine fg2(id)
        integer, intent(in) :: id
        integer :: ierr
        real :: con1,eh,ds
        real(8) GravConstant,ih,Hillradius,orbitaldistance
        real PeriodinYS
        real(8) EarthRadius,vesc,Planetesimalradius,Planetesimaldensity,Planetesimalmass
        real(8) s1,s2,s3,eps1,eps2,gamma,t1,t2,t3,c,beta,b
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        EarthRadius = 6.3781*10**8
        GravConstant = 6.67259/(10**(8))
        orbitaldistance = s%xtra4

        Hillradius = s%xtra7
        PeriodinYS = s%xtra9*3.154e+7
        Planetesimalradius = s%xtra1

        Planetesimaldensity = 2.5
        Planetesimalmass = 4./3.*pi*Planetesimaldensity*Planetesimalradius**3
        ds = (s% R_center + Planetesimalradius)/Hillradius
        vesc = sqrt(2*GravConstant*Planetesimalmass/Planetesimalradius)
        ih = vesc*PeriodinYS/(sqrt(3.)*Hillradius*2*pi)
        if(ih .ge. 2)then
          eh = 2*ih
        else
          eh = 2
        endif
        eps1 = 0.2*ih
        if(ds .le. 0.0284)then
          eps2 = eps1+sqrt(2.2*ds)/ih
          s2 = 6.69/(ds*ih)*(erf(0.25/ih)-erf(sqrt(2.2*ds)/ih))
        else if(ds .gt. 0.0284)then
          eps2 = eps1+0.25/ih
          s2 = 0
        ENDIF
        s1 = 4.*exp(eps1**2)/ds**(3./2)*(exp(-eps1**2)-exp(-eps2**2)+eps1*sqrt(pi)*(erf(eps1)-erf(eps2)))
        s3 = 1+3.31/(ds*ih**2)

         beta =0.5*pi/(1+exp((LOG(ih)-0.131)/0.343+1./6.*((LOG(ih)-0.131)/0.343)**3))
        gamma =0.5*pi/(1+exp((LOG(ih)+0.218)/0.244+1./6.*((LOG(ih)-3.8)/2)**3))
        t1 = ds**(0.04*(ih**0.25-1))*s1
        t2 = s2
        t3 = s3*ds**(0.08*(1-ih))
        c = 0.666 + 0.78*exp(-0.5*((LOG(ih)-0.172)/0.437)**2)
        B = 1. + 0.422*exp(-1./24.*((LOG(ih)+0.109)/0.231)**4)


        if(ih .ge. 2)then
          s% xtra8 = B*((s1+s2)*SIN(beta) +s3*COS(beta))
        else
          s% xtra8 = c*(t3*COS(gamma)+SIN(gamma)*(t1+t2))
        ENDIF



      end subroutine fg2


      function Velocity(position,mass,initialvelocity,initialposition) result(ValueofVelocity)
        real(dp)  ValueofVelocity
        real(dp)  position,mass,initialvelocity,initialposition

        ValueofVelocity = sqrt(initialvelocity**2+2*6.67/(10**8)*mass/position-2*mass*6.67/(10**8)/initialposition)

        return

     end function Velocity

     function Integrand(position,mass,initialvelocity,initialposition,temperature,rhoair) result(ValueofIntegrand)

       real(dp)  ValueofIntegrand
       real(dp)  position,mass,initialvelocity,initialposition,temperature,rhoair
       real(dp)  density,cp,tfus,e0

       cp = 4.19*10**7
       e0 = 2.8D10
       tfus = 373.
       density = 1.

       ValueofIntegrand = 1/(density*(cp*tfus+e0))*&
       &(5.67/10**5*temperature**4/velocity(position,mass,initialvelocity,initialposition)+&
       0.25*0.25*rhoair*velocity(position,mass,initialvelocity,initialposition)**2)

       return

    end function Integrand


    subroutine MetalsDeposition(id, ierr)

       integer, intent(in) :: id
       integer, intent(out) :: ierr
       type (star_info), pointer :: s
       real(dp) initialposition, initialvelocity,materialstrength
       real(dp) size,sunr,deltam,beforesize,fraction
        real(dp) temperature, rhoair,onemass
        real(dp) remainingmass,beforemetallicity,beforemetallicity2
        real(dp) accretedmass, factor,dummy1,diff_fraction,totalcomp

       integer i, Layer_number, j,m

       sunr=6.59700*10**10
       ierr = 0
       call star_ptr(id, s, ierr)
       if (ierr /= 0) return
       m= 1

       size = s%xtra1
       onemass = 4./3.*pi*s%xtra17*size**3


       accretedmass = s% dt/secyer*5.97d27/(10**6)
       factor = accretedmass / onemass

       initialvelocity = 100000
       initialposition = 2*sunr*10**s% log_surface_radius
       fraction = 0

       materialstrength = 1d6



       do while(m .LE. s% nz)


         deltam=0
        beforesize = size

        if(m .eq. s% nz)THEN
          size = size +Integrand(s% r(m),s%mstar,initialvelocity,initialposition,exp(s%lnT(m)),exp(s%lnd(m)))&
          &*(s%R_center-s%r(m))
        else
          size = size +Integrand(s% r(m),s%mstar,initialvelocity,initialposition,exp(s%lnT(m)),exp(s%lnd(m)))&
          &*(s%r(m+1)-s%r(m))
        ENDIF

        if(0.5*exp(s%lnd(m))*velocity(s%r(m),s%mstar,initialvelocity,initialposition)**2 .ge. materialstrength)then
          !size = size/2
        size = 0
          !materialstrength = materialstrength * (4./3.*pi*beforesize**3*Halley%density)/(4./3.*pi*size**3*Halley%density)
        endif

        deltam=4./3.*pi*beforesize**3*s%xtra17-4./3.*pi*size**3*s%xtra17



        if(beforesize .le. 0) then
          deltam = 0
        ENDIF

        s%xtra1_array(m) = deltam

        fraction = fraction + deltam

        diff_fraction = deltam/onemass*s% dt/secyer*s%xtra5

        beforemetallicity = s%xa(6,m)
        s%xa(6,m)=(s%xa(6,m)*s%dm(m)+diff_fraction)/s%dm(m)
        if(s%xa(6,m) .ge. 1)then
          remainingmass = (s%xa(6,m) -1)*s%dm(m)
          s%xa(6,m) = 1

          j = m+1
          do while(j .LE. s%nz)
            !write(*,*) "-------start-------"
            !write(*,*) "remaining mass ", remainingmass
            !write(*,*) "dm ", s%dm(j)
            !write(*,*) "m ", m
            !write(*,*) "j ", j
            beforemetallicity2 = s%xa(6,j)
            s%xa(6,j)=(s%xa(6,j)*s%dm(j)+remainingmass)/s%dm(j)
            !write(*,*) "before metallicity",beforemetallicity2
            !write(*,*) "after metallicity",s%xa(6,j)
            if(s%xa(6,j) .ge. 1)then
              s%xa(6,j) = 1
            ENDIF
            remainingmass = remainingmass - (s%xa(6,j)-beforemetallicity2)*s%dm(j)
            if(remainingmass .le. 0)then
              !write(*,*) "remaining mass le0", remainingmass
              exit
            ENDIF

            if(j .eq. s%nz)then
              !write(*,*) "-----------------diocane---------",remainingmass
              if(remainingmass .ge. 0)then
                !write(*,*) "-----------------core added---------"
                !s% xmstar = s% xmstar + remainingmass
                fraction = fraction - onemass*remainingmass/s%xtra5*secyer/s% dt
              ENDIF
            ENDIF
            j = j+1
          enddo
        ENDIF
            s%xa(1,m)=0.720000000000*(1-s%xa(6,m))
        s%xa(2,m)=0.0
        s%xa(3,m)=1-s%xa(6,m)-s%xa(1,m)
        s%xa(4,m)=0.0
        s%xa(5,m)=0.0
        s%xa(7,m)=0.0
        s%xa(8,m)=0.0

        totalcomp = s%xa(1,m)+s%xa(2,m)+s%xa(3,m)+s%xa(4,m)+s%xa(5,m)&
          &+s%xa(6,m)+s%xa(7,m)+s%xa(8,m)
!        write(*,*) "Metallicity",totalcomp


  m = m+1
enddo




        s%xtra16 = fraction/onemass



            s% accrete_same_as_surface = .false.
            s% accrete_given_mass_fractions = .true.



            dummy1 = 0.03


             !write(*,*) "metals deposition------mass change metals deposition",dummy1

            s% num_accretion_species = 3

            s% accretion_species_id(1) = 'h1   '
            !s% accretion_species_xa(1) = 0.1
            s% accretion_species_xa(1) =0.72*(1-dummy1)

            s% accretion_species_id(2) = 'he4  '
            !s% accretion_species_xa(2) = 0.1
            s% accretion_species_xa(2) = 0.28*(1-dummy1)

            s% accretion_species_id(3) = 'o16  '
            !s% accretion_species_xa(3) = 0.8
            s% accretion_species_xa(3) =  dummy1

!          write(*,*) "metals deposition------mass change metals deposition",&
!          &s%xtra16*s% dt/secyer*5.97d27/(10**6)/(s%xmstar - s%xmstar_older)





        !write(*,*) "Metallicity",fraction/onemass*accretedmass*1/(s%mass_change_older)

  end subroutine MetalsDeposition





         integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr,Layers_number,m
         real(8) Hillradius,Japanese,orbitaldistance,Planetesimaldensity,Planetesimalradius
         real(8) w,Compare,bondi,soundspeed,orbitaldistanceinau
         real(8) periodinYS,t0,p0,tsurf,psurf


         type (star_info), pointer :: s


         ierr = 0
         call star_ptr(id, s, ierr)
                 !s% xtra9 = 0
            if (ierr /= 0) return

        orbitaldistanceinau = 5.2

        orbitaldistance = orbitaldistanceinau*1.496*10**13

        PeriodinYS = 11.86*(orbitaldistanceinau/5.2)**(1.5)

        Planetesimalradius = 10**7
        Planetesimaldensity = 2.5

        s%xtra17 = Planetesimaldensity
        s%xtra4 = orbitaldistance
        s%xtra9 = PeriodinYS
        s%xtra1 = Planetesimalradius

        t0 = 150
        p0 = 0.27

        tsurf =t0*(orbitaldistanceinau / 5.2)**(-3./7.)
        psurf = p0*(orbitaldistanceinau / 5.2)**(-45./14.)

        call captureradius(id)

        s%which_atm_option = 'fixed_Psurf_and_Tsurf'
        s%atm_fixed_Psurf = psurf
           s%atm_fixed_Tsurf = tsurf
            s%use_atm_PT_at_center_of_surface_cell = .true.


        write(*,*) 'T surf',tsurf
        write(*,*) 'P surf',psurf
        write(*,*) 'start'


      end function extras_start_step



            subroutine captureradius(id)
                    integer, intent(in) :: id
                    integer :: ierr
                    integer :: Layers_number,m,j
                    real(8) Hillradius,Japanese,Planetesimalradius,planetesimaldensity
                    real(8) alpha,nablaad,opacity,denominator,rrcb
                    real(8) claudiodensity,gamma,rhorrcb,exponent,rprime 

                    type (star_info), pointer :: s
                    ierr = 0
                    call star_ptr(id, s, ierr)
                    if (ierr /= 0) return

                    gamma = 1.66
                    nablaad = (gamma - 1.)/gamma
                    opacity = 0.01
                    denominator = 64*PI*6.6743e-8*s%mstar*s%xtra26**4*5.6704e-5
                    exponent = 1./(gamma-1.)
                    rprime = s% xtra2*nablaad/1.53


                    planetesimalradius = 1d5
                    planetesimaldensity =1.
                    Hillradius = s%xtra7
                    Layers_number = s% nz

                    alpha = 6.6743e-8*s%mstar*s%xtra21/(s%xtra22*s% xtra2)
                    
                    rrcb = 1
                   Japanese = 3./2.*Hillradius/Planetesimalradius*1/Planetesimaldensity
                   rhorrcb = s%xtra21*exp(alpha*(s%xtra2/rrcb-1))

                   


                   j = 1
                   m = j
                   do while(Japanese * exp(s%lnd(j)) - 1 .LE. 0)
          

                         j = j+1
                         if (j .EQ. Layers_number-1) then
                          exit
                         ENDIF
                       enddo
                       s% xtra20 = s%R(j)



                   j = 1
                   m = j
                   claudiodensity = s%xtra21
                   
                   do while(Japanese * claudiodensity - 1 .LE. 0)

                         j = j+1
                         claudiodensity=s%xtra21*exp(alpha*(s%xtra2/s%R(j)-1))
                         if(s%R(j) .le. rrcb) then
                          write(*,*) 'OOOOOOOO', rrcb/(1d10)
                          claudiodensity = rhorrcb*(1+rprime/s%R(j)-rprime/rrcb)**exponent
                         endif

                         if (j .EQ. Layers_number-1) then
                          exit
                         ENDIF
                       enddo
                       s% xtra18 = s%R(j)




                       j = Layers_number
                       m = j
                       planetesimalradius = 1d7
     
                       Japanese = 3./2.*Hillradius/Planetesimalradius*1/Planetesimaldensity
     
                       j = 1
                       m = j
                       do while(Japanese * exp(s%lnd(j)) - 1 .LE. 0)
              
     
                             j = j+1
                             if (j .EQ. Layers_number-1) then
                              exit
                             ENDIF
                           enddo
                           s%xtra28 = s%R(j)
     
     






                  j = Layers_number
                  m = j
                  planetesimalradius = 1d7

                  Japanese = 3./2.*Hillradius/Planetesimalradius*1/2.5

                  j = 1
                  m = j
                  do while(Japanese * exp(s%lnd(j)) - 1 .LE. 0)
         

                        j = j+1
                        if (j .EQ. Layers_number-1) then
                         exit
                        ENDIF
                      enddo
                      s%xtra6 = s%R(j)



                  j = 1
                  m = j
                
                  claudiodensity = s%xtra21
                  
                  do while(Japanese * claudiodensity - 1 .LE. 0)

                        j = j+1
                        claudiodensity=s%xtra21*exp(alpha*(s%xtra2/s%R(j)-1))
                        if(s%R(j) .le. rrcb) then
                          claudiodensity = rhorrcb*(1+rprime/s%R(j)-rprime/rrcb)**exponent
                         endif

                        if (j .EQ. Layers_number-1) then
                         exit
                        ENDIF
                      enddo
                      s% xtra23 = s%R(j)





                      j = Layers_number
                      m = j
                      planetesimalradius = 1d3
    
                      Japanese = 3./2.*Hillradius/Planetesimalradius*1/Planetesimaldensity
    
                      j = 1
                      m = j
                      do while(Japanese * exp(s%lnd(j)) - 1 .LE. 0)
             
    
                            j = j+1
                            if (j .EQ. Layers_number-1) then
                             exit
                            ENDIF
                          enddo
                          s%xtra19 = s%R(j)
    
    
    
                      j = 1
                      m = j
                      claudiodensity = s%xtra21
                      
                      do while(Japanese * claudiodensity - 1 .LE. 0)
    
                            j = j+1
                            claudiodensity=s%xtra21*exp(alpha*(s%xtra2/s%R(j)-1))
                            if(s%R(j) .le. rrcb) then
                              claudiodensity = rhorrcb*(1+rprime/s%R(j)-rprime/rrcb)**exponent
                             endif
    
                            if (j .EQ. Layers_number-1) then
                             exit
                            ENDIF
                          enddo
                          s% xtra24 = s%R(j)
    
                          write(*,*) 'Rcrb ', rrcb/(1d10)
                          write(*,*) 'Lcenter ', s%L_center





            end subroutine captureradius



      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         s%xtra5 =0.
         s%xtra25 =0.



         call system_clock(time0,clock_rate)
         if (.not. restart) then
          call alloc_extra_info(s)
        else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup


      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         character (len=strlen) :: test
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call system_clock(time1,clock_rate)
         dt = dble(time1 - time0) / clock_rate / 60
         call GET_ENVIRONMENT_VARIABLE( &
            "MESA_TEST_SUITE_CHECK_RUNTIME", test, status=ierr, trim_name=.true.)
         if (ierr == 0 .and. trim(test) == 'true' .and. dt > 1.5*expected_runtime) then
            write(*,'(/,a70,2f12.1,99i10/)') &
               'failed: EXCESSIVE runtime, prev time, retries, backups, steps', &
               dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
         else
            write(*,'(/,a50,2f12.1,99i10/)') 'runtime, prev time, retries, backups, steps', &
               dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
         end if
         ierr = 0
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr,m
         real(8) sunr,Hillradius,orbitaldistance,tollerance
         real(8) soundspeed,Bondi,newradius,kilometer,orbitaldistanceinau
         type (star_info), pointer :: s
         ierr = 0
           call star_ptr(id, s, ierr)



          kilometer = 10**10


	       sunr=6.59700*10**10
         newradius=sunr*10**s% log_surface_radius

        orbitaldistance = s%xtra4
        orbitaldistanceinau = orbitaldistance / (1.496*10**13)
        Hillradius = orbitaldistance*(s% mstar / (3*Msun))**(1./3.)
        sunr=6.59700*10**10


        soundspeed = 0.4*(orbitaldistanceinau / 10.)**(3./14.)*10**5
        Bondi = 6.67/(10**8)*s%mstar/(2*soundspeed**2+4*6.67/(10**8)*s%mstar/Hillradius)

        s% xtra2 = Bondi

         !write(*,*) 'Check before if--mass change is', s% mass_change
         write(*,*) 'Check before if--ratio is', s% xtra2/newradius
        if(s% model_number .le. 310)then
                tollerance = 10000000.
        else
                tollerance = 1.05
         endif
         if(S% XTRA10 == 0)THEN
          if(s% xtra2/newradius .ge.tollerance)then

             write(*,*) 'Check--ratio greater than 1.10'
             write(*,*) 'Check--ratio is',s% xtra2/newradius
             write(*,*) '------------Redo the step--------------'
             s% mass_change = s% mass_change *4


             s% dt = s% dt/1.5


             extras_check_model = retry



          else
            if(s% model_number .le. 1000)then
               s% dt = 1.5*s% dt
            else
              s%dt = s%dt
            endif
          extras_check_model = keep_going
          ENDIF
        else
          extras_check_model = keep_going

        ENDIF







        if(s% mstar.GE. 100*5.97d27) then
          extras_check_model = terminate
          endif

          !if(s% star_age/1000000 .GT. 3.)then
          !              !extras_check_model = terminate
          !endif

 if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


            integer function how_many_extra_history_columns(id, id_extra)
               integer, intent(in) :: id, id_extra
               integer :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               how_many_extra_history_columns = 16
            end function how_many_extra_history_columns


            subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
               integer, intent(in) :: id, id_extra, n
               character (len=maxlen_history_column_name) :: names(n)
               real(dp) :: vals(n)
               integer :: m
               integer, intent(out) :: ierr
               real(8) :: c0,c1,c2,EarthRadius
               real(8) :: PeriodinYS,orbitaldistance,Hillradius,ds,Fg,mdotz,sunr,newradius,sigmaz
               real(8) :: sigmazinit,accretedmass,initialMcenter,kilometer,soundspeed,bondi
               real(8) :: zcorei,zcoreii,xcorei,xcoreii,mhydro,mhydro2
               real(8) :: m_0,m_1,m_2,m_3,m_4,m_5,mixing_type,gapdensity,deltar2
               real(8) :: regions,region_0,region_1,region_2,region_3,region_4,region_5
               real(8) :: outrho,outp,outt,rhoderivative,deltar1
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return

                      kilometer = 10**5
                      orbitaldistance =s%xtra4
                      Hillradius = s%xtra7
                      mhydro = 0
                      mhydro2 = 0


                      outrho = exp(s%lnd(1))
                      outp = exp(s%lnp(1))
                      outt = exp(s%lnt(1))

                      rhoderivative = (exp(s%lnd(1)) - exp(s%lnd(2)))/(s%R(1)-s%R(2))


                      soundspeed = 0.4*0.52**(3./14.)*10**5
                      Bondi = 6.67/(10**8)*s%mstar/(2*soundspeed**2+4*6.67/(10**8)*s%mstar/Hillradius)

                      s% xtra2 = Bondi


               	       names(1)= "Radius_Jupiterunits"
               	       vals(1)=  (10**(s% log_surface_radius)*Rsun)/(69911*kilometer)

                       names(2)= "Rcap_Jupiterunits"
                       vals(2)= s%xtra6/(69911*kilometer)

               	       names(3) = "Planet_Mass_earthunits"
               	       vals(3) =  s% mstar/5.97e27

                       names(3) = "Planet_Mass_earthunits"
                       vals(3) =  s% mstar/5.97e27

                       names(4) = "Mdot(Eart/yr)"
                       vals(4) =  (s% mstar-s% mstar_old)/(5.97e27*s%dt/secyer)

                       names(5) = "captradius1d7"
                       vals(5) =  s%xtra28

                       names(6) = "captradius1d7claudio"
                       vals(6) =  s%xtra23


                       names(7) = "captradius1d3"
                       vals(7) =  s%xtra19

                       names(8) = "captradius1d3claudio"
                       vals(8) =  s%xtra24


                        names(9) = "captradius1d5"
                       vals(9) =  s%xtra20

                       names(10) = "captradius1d5claudio"
                       vals(10) =  s%xtra18



                       names(11) = "deltar2(au)"
                       vals(11) =  deltar2/(1.496*10**13)

                       names(12) = "solidaccretedmass"
                       vals(12) =  s% xtra25/5.97e27

                       names(13) = "Mcore"
                       vals(13) =  s% M_center/5.97e27


                       names(14) = "Menvelope"
                       vals(14) =  s% xmstar/5.97e27

                       names(15) = "mdotsolid"
                       vals(15) =  s% xtra5

                       names(16) = "Menvelope"
                       vals(16) =  s% xtra5 * s%xtra16



            end subroutine data_for_extra_history_columns



      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return



      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr,m,j
        real(8) new_core_mass,core_avg_rho,core_avg_eps,sunr,oldradius,newradius,maxabsdot,middleradius,oldradius2,newradius2
        real(8) normalizingvalue,initialposition,initialvelocity
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)

        s% use_other_wind = .true.
        s% use_other_energy = .false.


         sunr=6.59700*10**10



         m =1
                do while(m .LE. s% nz)

                 if(0.5*exp(s%lnd(m))*velocity(s%r(m),s%mstar,initialvelocity,initialposition)**2 .ge. 1d6)then
                   write(*,*) 'fragmentation'
                   !s%xtra24 = s%R(m)
                   exit
                 endif




           m = m+1
         enddo

                  j =1
                do while(j .LE. s% nz)

                 if(0.5*exp(s%lnd(j))*velocity(s%r(j),s%mstar,initialvelocity,initialposition)**2 .ge. 1d8)then
                   write(*,*) 'fragmentation'
                   !s%xtra21 = s%R(j)
                   exit
                 endif




           j = j+1
         enddo



           !call MetalsDeposition(id, ierr)

         s%xtra21 = exp(s%lnd(1))
         s%xtra22 = exp(s%lnp(1))
         s%xtra26 = exp(s%lnt(1))



         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)


      end function extras_finish_step




      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
           end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      subroutine enxa ( n, x, en )

!*****************************************************************************80
!
!! ENXA computes the exponential integral En(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) EN(0:N), the function values.
!
        implicit none

        integer ( kind = 4 ) n

        real ( kind = 8 ) e1
        real ( kind = 8 ) ek
        real ( kind = 8 ) en(0:n)
        integer ( kind = 4 ) k
        real ( kind = 8 ) x

        en(0) = exp ( - x ) / x
        call e1xb ( x, e1 )

        en(1) = e1
        do k = 2, n
           ek = ( exp ( - x ) - x * e1 ) / ( k - 1.0D+00 )
           en(k) = ek
           e1 = ek
        end do

        return
      end subroutine enxa

      subroutine e1xb ( x, e1 )

!*****************************************************************************80
!
!! E1XB computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
!
        implicit none

        real ( kind = 8 ) e1
        real ( kind = 8 ) ga
        integer ( kind = 4 ) k
        integer ( kind = 4 ) m
        real ( kind = 8 ) r
        real ( kind = 8 ) t
        real ( kind = 8 ) t0
        real ( kind = 8 ) x

        if ( x == 0.0D+00 ) then

           e1 = 1.0D+300

        else if ( x <= 1.0D+00 ) then

           e1 = 1.0D+00
           r = 1.0D+00

           do k = 1, 25
              r = -r * k * x / ( k + 1.0D+00 )**2
              e1 = e1 + r
              if ( abs ( r ) <= abs ( e1 ) * 1.0D-15 ) then
                 exit
              end if
           end do

           ga = 0.5772156649015328D+00
           e1 = - ga - log ( x ) + x * e1

        else

           m = 20 + int ( 80.0D+00 / x )
           t0 = 0.0D+00
           do k = m, 1, -1
              t0 = k / ( 1.0D+00 + k / ( x + t0 ) )
           end do
           t = 1.0D+00 / ( x + t0 )
           e1 = exp ( -x ) * t

        end if

        return
      end subroutine e1xb


            subroutine opacity_factor_routine(id, ierr)
               use const_def, only: Rsun
               integer, intent(in) :: id
              integer, intent(out) :: ierr
               type (star_info), pointer :: s
               integer :: k
               ierr = 0
               call star_ptr(id, s, ierr)
              if (ierr /= 0) return
                s% extra_opacity_factor(1:s% nz) = s% opacity_factor

            end subroutine opacity_factor_routine


            ! adds contribution of grains to the opacity using the Valencia+ (2013) description
               subroutine kap_freedman_grains(id, k, handle, zbar, X, Z, logRho, logT, &
                                     species, chem_id, net_iso, xa, lnfree_e, &
                                     d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                                     kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

                   use eos_lib
                   use eos_def
                   use kap_lib, only: kap_get_Type1
                   use chem_def

                   integer, intent(in) :: id, k, handle, species
                   real(dp), intent(in) :: zbar, x, Z, logRho, logT, &
                                           lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, xa(:)
                   real(dp), intent(out) :: kap, dln_kap_dlnRho, dln_kap_dlnT  ! total opacity
                   integer, intent(out) :: ierr
                   real(dp) :: kap_gas, dln_kap_gas_dlnRho, dln_kap_gas_dlnT  ! gas opacity
                   real(dp) :: logT_eps, logRho_eps, eps  ! used for derivatives
                   real(dp) :: kapRho, kapT  ! used for derivatives
                   real(dp) :: T_6, logR_bar, logT1_star, logT2_star  ! used for the grain opacity
                   real(dp) :: dumvar1, dumvar2  ! dummy variables
                   real(dp) :: abar  ! average atomic weight
                   integer :: i
                   integer, pointer :: chem_id(:), net_iso(:)
                   logical, parameter :: debug = .false.
                   type (star_info), pointer :: s
                   call star_ptr(id, s, ierr)

                   ! calculate average atomic weight
           !$omp parallel default(shared)
           !$omp do reduction(+:abar)
                   do i = 1, species
                       abar = abar + xa(i) / (chem_isos % N (chem_id(i)) + chem_isos % Z (chem_id(i)))
                   end do
           !$omp end do
           !$omp end parallel
                   abar = 1d0 / abar

                   eps = 1d-3
                   logRho_eps = logRho + eps
                   logT_eps = logT + eps

                   T_6 = exp10_cr(logT) / (10**6)
                   logR_bar = log10_cr(exp10_cr(logRho) / (T_6**3))
                   logT1_star = 0.0245 * logR_bar + 1.971
                   logT2_star = 0.0245 * logR_bar + 3.221
                   ! typical values are logT1_star = 2.0 and logT2_star = 3.2

                   ! get the opacity without grains
                   call kap_get_Type1(s% kap_handle, zbar, X, Z, logRho, logT, &
                                      lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                                      kap_gas, dln_kap_gas_dlnRho, dln_kap_gas_dlnT, ierr)
                   if(ierr /= 0) write(*, *) 'kap_freedman_grains failed to get kappa_ng'

                   ! get the results
                   ! full grain contribution in this range
                   if(logT .lt. logT1_star) then
                       call kap_get_Type1(s% kap_handle, zbar, X, Z, logRho_eps, logT, &
                                          lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                                          kapRho, dumvar1, dumvar2, ierr)
                       if(ierr /= 0) write(*, *) 'kap_freedman_grains failed to get kapRho'
                       call kap_get_Type1(s% kap_handle, zbar, X, Z, logRho, logT_eps, &
                                          lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                                          kapT, dumvar1, dumvar2, ierr)
                       if(ierr /= 0) write(*, *) 'kap_freedman_grains failed to get kapT'

                       kap = kap_gas + get_kap_grains(logT)
                       kapRho = kapRho + get_kap_grains(logT)
                       kapT = kapT + get_kap_grains(logT_eps)
                       dln_kap_dlnRho = log_cr(kapRho / kap) / (ln10 * eps)
                       dln_kap_dlnT = log_cr(kapT / kap) / (ln10 * eps)
                   ! no grains if the temperature is high
                   else if(logT .gt. logT2_star) then
                       kap = kap_gas
                       dln_kap_dlnRho = dln_kap_gas_dlnRho
                       dln_kap_dlnT = dln_kap_gas_dlnT
                   ! linearly interpolate in the intermediate region
                   else
                       kap = get_kap_linear(logT, logRho)
                       kapRho = get_kap_linear(logT, logRho_eps)
                       kapT = get_kap_linear(logT_eps, logRho)
                       dln_kap_dlnRho = log_cr(kapRho / kap) / (ln10 * eps)
                       dln_kap_dlnT = log_cr(kapT / kap) / (ln10 * eps)
                   end if

                   contains


                   real(dp) function get_kap_grains(logT_in) result(kap_out)
                       real(dp), intent(in) :: logT_in
                       kap_out = s% x_ctrl(10) * exp10_cr(0.430 + 1.3143 * (logT_in - 2.85))
                   end function get_kap_grains


                   real(dp) function get_kap_linear(logT_in, logRho_in) result(kap_out)
                       real(dp), intent(in) :: logT_in, logRho_in
                       real(dp) :: T_6, logR_bar, logT1_star, logT2_star
                       real(dp) :: kapT1_gas, dln_kapT1_ng_dlnRho, dln_kapT1_ng_dlnT
                       real(dp) :: kapT2_gas, dln_kapT2_ng_dlnRho, dln_kapT2_ng_dlnT
                       real(dp) :: kapT1, kapT2
                       real(dp) :: dummy1, dummy2

                       T_6 = exp10_cr(logT_in) / (10**6)
                       logR_bar = log10_cr(exp10_cr(logRho_in) / (T_6**3))
                       logT1_star = 0.0245 * logR_bar + 1.971
                       logT2_star = 0.0245 * logR_bar + 3.221

                       call kap_get_Type1(s% kap_handle, zbar, X, Z, logRho_in, logT1_star, &
                                          lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                                          kapT1_gas, dummy1, dummy2, ierr)
                       if(ierr /= 0) write(*, *) 'kap_freedman_grains failed to get kapT1_gas'

                       call kap_get_Type1(s% kap_handle, zbar, X, Z, logRho_in, logT2_star, &
                                          lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                                          kapT2_gas, dummy1, dummy2, ierr)
                       if(ierr /= 0) write(*, *) 'kap_freedman_grains failed to get kapT2_gas'

                       kapT1 = kapT1_gas + s% x_ctrl(10) * exp10_cr(0.430 + 1.3143 * (logT1_star - 2.85))
                       kapT2 = kapT2_gas
                       kap_out = kapT1 + ((logT_in - logT1_star) * (kapT2 - kapT1) / (logT2_star - logT1_star))
                   end function get_kap_linear


               end subroutine kap_freedman_grains



      end module run_star_extras
