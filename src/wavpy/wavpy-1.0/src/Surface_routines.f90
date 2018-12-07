
 !   Copyright (C) 2017  Fran Fabra

 !   This program is free software: you can redistribute it and/or modify
 !   it under the terms of the GNU General Public License as published by
 !   the Free Software Foundation, either version 3 of the License, or
 !   (at your option) any later version.

 !   This program is distributed in the hope that it will be useful,
 !   but WITHOUT ANY WARRANTY; without even the implied warranty of
 !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 !   GNU General Public License for more details at:

 !   <http://www.gnu.org/licenses/>

subroutine compute_epsilon_sea_water(sal, temp, freq_GHz, epsilon_real, epsilon_imag)

  implicit none
  real*8 :: sal, temp, freq_GHz, epsilon_real, epsilon_imag
  real*8 :: freq,E0,Eswinf,a,Esw0,Esw00,Tausw0,b,Tausw,si1,inc,phi,si,Epimn,Epimd,C,pi

  freq = freq_GHz*1d9
  pi = acos(-1.d0)
  E0 = 8.854d-12
  Eswinf = 4.9d0
  Esw00 = 87.174-(1.949d-1)*temp-(1.279d-2)*(temp**2d0)+(2.491d-4)*(temp**3d0)
  a = 1.0+(1.613d-5)*temp*sal-(3.656d-3)*sal+(3.21e-5)*(sal**2d0)-(4.232d-7)*(sal**3d0)
  Esw0 = Esw00*a*1d0
  Tausw0 = (1.1109d-10)-(3.824d-12)*temp+6.238d-14*(temp**2d0)-(5.096d-16)*(temp**3d0)
  b = 1.0+(2.282d-5)*temp*sal-(7.638d-4)*sal-(7.760d-6)*(sal**2d0)+(1.105d-8)*(sal**3d0)
  Tausw = (Tausw0/(2d0*pi))*b
  si1 = sal*(0.18252-1.4619d-3*sal+2.093d-5*(sal**2d0)-1.282d-7*(sal**3d0))
  inc = 25-temp
  phi = inc*(2.033d-2+1.266d-4*inc+2.464d-6*(inc**2d0)-sal*(1.849d-5-2.551d-7*inc+2.551d-8*(inc**2d0)))
  si = si1*exp(-phi)
  Epimn = (2*pi*freq*Tausw*(Esw0-Eswinf))
  Epimd = 1d0+((2*pi*freq*Tausw)**2d0)
  C = (si/(2*pi*E0*freq))

  epsilon_real = Eswinf+((Esw0-Eswinf)/(1d0+((2d0*pi*freq*Tausw)**2)))
  epsilon_imag = ((Epimn/(Epimd*1d0))+C)

end subroutine compute_epsilon_sea_water


subroutine compute_epsilon_sea_ice(v_b, epsilon_real, epsilon_imag)

  implicit none 
  real*8 :: v_b, epsilon_real, epsilon_imag

  !HALLIKAINEN AND WINEBRENNER
  epsilon_real = 3.12 + 0.009*v_b  !Vb < 70, 1 GHz 
  epsilon_imag = 0.04 + 0.005*v_b  !Vb < 70, 1 GHz 

end subroutine compute_epsilon_sea_ice


subroutine compute_epsilon_dry_snow(rho_s, epsilon_real, epsilon_imag)

  implicit none 
  real*8 :: rho_s, epsilon_real, epsilon_imag
  real*8 :: epsr_ice, epsi_ice, v_i

  !rho_s = 0.3       !snow density in g cm**(-3)   (it rarely exceeds 0.5)        300 Kg/m**(-3)  in DOMEX 2004 
  epsr_ice = 2.95   !TABLE E.1 Ulaby (2028)  L-band  ->  2.95 
  v_i = rho_s/0.916 !volume fraction of ice in the dry snow, 0.916 is the density of pure ice. Ulaby (2061) 
  epsi_ice = 0.001  !from fig E.3, 0.000405  from E.31  Ulaby (2026)

  epsilon_real = (1.0 + 0.51*rho_s)**3d0 
  epsilon_imag = 3.0*v_i*epsi_ice*(epsilon_real**2d0)*(2.0*epsilon_real + 1.0)/ & 
((epsr_ice + 2.0*epsilon_real)*(epsr_ice + 2.0*(epsilon_real**2d0))) !Fung (420) 

end subroutine compute_epsilon_dry_snow
 

subroutine compute_epsilon_wet_snow(rho_s, m_v, freq_GHz, epsilon_real, epsilon_imag)

  implicit none 
  real*8 :: rho_s, m_v, freq_GHz, epsilon_real, epsilon_imag
  real*8 :: freq_0_GHz

  freq_0_GHz = 9.07 
  !m_v = volumen fraction of liquid water in %. Ranges from 1 to 12 %
  epsilon_real = 1.0 + 1.83*rho_s + 0.02*(m_v**1.015) + (0.073*(m_v**1.31))/(1.0 + (freq_GHz/freq_0_GHz)**2d0)  !Fung (420)
  epsilon_imag = (0.073*(freq_GHz/freq_0_GHz)*(m_v**1.31))/(1.0 + (freq_GHz/freq_0_GHz)**2d0)  !Fung (420) 

end subroutine compute_epsilon_wet_snow


subroutine compute_fresnel_linear_refl(theta, epsilon1_real, epsilon1_imag, epsilon2_real, &
epsilon2_imag, rvv_real, rvv_imag, rhh_real, rhh_imag)

  implicit none
  real*8 :: theta, epsilon1_real, epsilon1_imag, epsilon2_real, epsilon2_imag
  real*8 :: rvv_real, rvv_imag, rhh_real, rhh_imag
  complex*16 epsilon1, epsilon2, rvv, rhh

  epsilon1 = cmplx(epsilon1_real,epsilon1_imag)
  epsilon2 = cmplx(epsilon2_real,epsilon2_imag)

  rhh = (epsilon1*cos(theta) - sqrt(epsilon1*epsilon2 - (epsilon1*sin(theta))**2))/&
(epsilon1*cos(theta) + sqrt(epsilon1*epsilon2 - (epsilon1*sin(theta))**2))
  rvv = (epsilon2*cos(theta) - sqrt(epsilon1*epsilon2 - (epsilon1*sin(theta))**2))/&
(epsilon2*cos(theta) + sqrt(epsilon1*epsilon2 - (epsilon1*sin(theta))**2))

  rvv_real = real(rvv)
  rvv_imag = aimag(rvv)
  rhh_real = real(rhh)
  rhh_imag = aimag(rhh)

end subroutine compute_fresnel_linear_refl


subroutine compute_fresnel_linear_trans(theta, epsilon1_real, epsilon1_imag, epsilon2_real, &
epsilon2_imag, tvv_real, tvv_imag, thh_real, thh_imag)

  implicit none
  real*8 :: theta, epsilon1_real, epsilon1_imag, epsilon2_real, epsilon2_imag
  real*8 :: tvv_real, tvv_imag, thh_real, thh_imag
  complex*16 epsilon1, epsilon2, tvv, thh

  epsilon1 = cmplx(epsilon1_real,epsilon1_imag)
  epsilon2 = cmplx(epsilon2_real,epsilon2_imag)

  thh = 2.0*cos(theta)/(cos(theta) + sqrt((epsilon2/epsilon1) - (sin(theta)**2)))
  tvv = 2.0*cos(theta)/(cos(theta) + sqrt((epsilon1/epsilon2) - (((epsilon1/epsilon2)*sin(theta))**2)))

  tvv_real = real(tvv)
  tvv_imag = aimag(tvv)
  thh_real = real(thh)
  thh_imag = aimag(thh)

end subroutine compute_fresnel_linear_trans


subroutine sea_spectrum_Elfouhaily(psi, kx, ky, deltak, n, U10, thetaD, Ome)
!--------------------------------------------
! For a given set of sea state conditions,
! this program generates the 
! OMNIDIRECTIONAL sea spectra, S(k),
! as suggested in [Elfouhaily et al., 1997] 
! Also the DIRECTIONAL spectrum using 'unified'
! 'angular spreading functions'.
! 
! ASSUMPTIONS:
!   -Deep waters (c=sqrt(g/k))
!--------------------------------------------
!
! psi(n,n)=spectra (IN empty/OUT)
! kx(n),ky(n)=arrays of kx and ky (IN empty/OUT)
! deltak=step in kx and ky (SAME)
! n= half-length of spectra (i.e.length of pos or/and negative)
! U10: wind speed at 10m above the surface [m/s]
! theta: angle between wind and waves [deg]
! omega: wave age > 0.83
!---------------------------------------------
! Estel Cardellach, 12/17/2002
!                 , 02/03/2007, now subroutine
!---------------------------------------------
  implicit none
  real*8, INTENT(IN) :: U10, thetaD, Ome, deltak
  integer, INTENT(IN) :: n ! n is half-length of spectra (pos/neg)
  real*8, dimension(2*n,2*n), INTENT(INOUT) :: psi
  real*8, dimension(2*n), INTENT(INOUT) :: kx, ky
  real*8, parameter :: g = 9.81     ! gravitational acceleration
  integer i, ii, j, ix, iy, length_u, length_t, length_o
  real*8 S, Bl, Bh, Del, Phi, theta
  real*8 k, c, cp, alphap, kp, jonswap, gamma, Oc
  real*8 cm, alpham, km, ufric, sigma2
  real*8 xp, xm, yp, ym, z, y, w
  real*8 a0, ap, am, ph

  !---------- Get arguments:
  theta = thetaD*3.14159/180.0
  !---------- Wind dependent parameters, k-independent  
  cp = U10/(Ome*1.0)
  kp = g/(cp**2)     !**ATTENTION: assuming Deep Waters**  
  km = 370.0
  !ufric = sqrt(1e-3*(1.10 + 0.04*U10))*U10  !** from web page reference AOMIP
  ufric = sqrt(1e-3*(0.81 + 0.065*U10))*U10  !** from Elfouhaily's code
  alphap = 6e-3*sqrt(Ome)
  Oc = Ome*cos(theta)
  if ((Oc >= 0.84) .and. (Oc < 1.0)) gamma = 1.7
  if ((Oc >= 1.0) .and. (Oc < 5.0)) gamma = 1.7 + 6.0*log10(Oc)
  sigma2 = (0.08*(1.0 + 4.0*(Oc)**(-3)))**2
  cm = 0.23
  if (ufric < cm) alpham = 1e-2*(1.0+log(ufric/cm))
  if (ufric > cm) alpham = 1e-2*(1.0+3.0*log(ufric/cm))
  !for spreading functions...
  a0 = log(2.0)/4.0
  ap = 4.0
  am = 0.13*ufric/cm
  !------------------ k=0
  psi(:,:) = 0.0d0
  psi(1,1) = 0.0d0
  !------------------LOOP IN kx
  do i=1,2*n
     if (i.le.(n + 1)) then
       ix = i
       kx(ix) = (i-1)*deltak
     else
       j = i-n-1
       ix = 2*n-(j-1)
       kx(ix) = -j*deltak
     endif
     !---------------LOOP IN ky
     do ii=1,2*n
        if ((ii.eq.1).and.(i.eq.1)) goto 5
        if (ii.le.(n+1)) then
          iy = ii
          ky(iy) = (ii-1)*deltak
        else
          j = ii-n-1
          iy = 2*n-(j-1)
          ky(iy) = -j*deltak
        endif
        k = sqrt(kx(ix)**2 + ky(iy)**2)
        ph = atan2(ky(iy), kx(ix))
        !--- now omnidirectional
        xp = sqrt(k/kp)
        xm = k/km
        c = sqrt(g/k*(1.0 + xm*xm))
        yp = Ome/sqrt(10.0)*(xp - 1.0)
        Bl = 0.50*alphap*cp/c*exp(-yp)
        ym = (xm - 1.0)**2/4.0
        Bh = 0.50*alpham*cm/c*exp(-ym)
        y = exp(-5.0/4.0 * (kp/k) * (kp/k))
        w = sqrt(k/kp) - 1.0
        z = exp(-0.5*w*w/sigma2)
        jonswap = y*gamma**z
        S = (Bl + Bh)*jonswap/k**3.0
        ! S is omnidirectional spectra
        !--- now directional factor
        Del = tanh(a0 + ap*(c/cp)**2.5 + am*(cm/c)**2.5)
        
        !--- TOTAL
        Phi = 1.0/(2.0*3.14159)*(1.0 + Del*cos(2.0*ph))
        psi(ix,iy) = 1.0/k *S*Phi
5       continue 
     end do
6 end do
end subroutine sea_spectrum_Elfouhaily

