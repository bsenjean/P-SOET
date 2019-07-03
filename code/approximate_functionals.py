#!/usr/bin/python 
# -*- coding: utf-8 -*-
# print "Hello, Python!";

# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


######################################
#
# This program provide all the derivations
# for the approximation in SOET:
# iBALDA, SIAMBALDA and 2LBALDA.
# See Ref. Senjean et al. Theoretical Chemistry Accounts (2018) 137:169
# for the derivations.
#
######################################
import numpy as np


##################
# Approximations #
##################

def correlation_BALDA(U,t,occ,beta,dbeta_dU):
   if occ <= 1.0:
       ec = - 2*t*beta*np.sin(np.pi*occ/beta)/np.pi + 4.0*t*np.sin(np.pi*occ/2.0)/np.pi - U*0.25*occ*occ
       dec_dU = dbeta_dU*(-2*t*np.sin(np.pi*occ/beta)/np.pi + 2*t*occ*np.cos(np.pi*occ/beta)/beta) - occ*occ*0.25
       dec_dn = - 2*t*np.cos(np.pi*occ/beta) + 2*t*np.cos(np.pi*occ/2.0) - U*0.5*occ
   else:
       ec = - 2*t*beta*np.sin(np.pi*(2-occ)/beta)/np.pi + U*(occ-1) + 4.0*t*np.sin(np.pi*(2-occ)/2.0)/np.pi - U*occ*occ*0.25
       dec_dU = dbeta_dU*(-2*t*np.sin(np.pi*(2-occ)/beta)/np.pi + 2*t*(2-occ)*np.cos(np.pi*(2-occ)/beta)/beta) + occ - 1 - occ*occ*0.25
       dec_dn = 2*t*np.cos(np.pi*(2.0-occ)/beta) - 2*t*np.cos(np.pi*(2.0-occ)/2.0) + U - U*0.5*occ
   dec_dt = ec/t - U*dec_dU/t
   return ec, dec_dU, dec_dt, dec_dn

def correlation_iBALDA(U,t,occ,beta,dbeta_dU):
   (Ecimp,dEcimp_dU,dEcimp_dt, dEcimp_dn)=correlation_BALDA(U,t,occ,beta,dbeta_dU)
   return Ecimp, dEcimp_dU, dEcimp_dt, dEcimp_dn

def correlation_SIAMBALDA(U,t,occ):
   if occ <=1 :
       Gamma = t*(1 + np.cos(np.pi*occ/2.0))/np.sin(np.pi*occ/2.0)
       dGamma_dn = - np.pi*Gamma/(2*np.sin(np.pi*occ/2.0))
   else:
       Gamma = t*(1 + np.cos(np.pi*(2-occ)/2.0))/np.sin(np.pi*(2-occ)/2.0)
       dGamma_dn = + np.pi*Gamma/(2*np.sin(np.pi*(2-occ)/2.0))
   dEcimp_dGamma = 0.0369*(U/Gamma)**2/np.pi - 3*0.0008*(U/Gamma)**4/np.pi**3
   dEcimp_dn = dGamma_dn*dEcimp_dGamma
   Ecimp = U**2*(-0.0369 + 0.0008*(U/(np.pi*Gamma))**2)/(np.pi*Gamma)
   dEcimp_dU = -2*0.0369*(U/Gamma)/np.pi + 4*0.0008*(U/Gamma)**3/(np.pi**3)
   dEcimp_dt = Ecimp/t - U*dEcimp_dU/t
   return Ecimp, dEcimp_dU, dEcimp_dt, dEcimp_dn

def correlation_2LBALDA(U,t,occ):
   """
   Second parametrization for the impurity Hubbard dimer,
   see Senjean et al. Theoretical Chemistry Accounts (2018) 137:169
   """
   occ_tmp = 0
   if U == 0:
     U =0.000001
   if occ == 1:
     occ = 0.99999
     occ_tmp = 1
   sign = lambda x: x and (1, -1)[x<0]
   rho = abs(occ - 1)
   # For the impurity case, u = U/4t and not U/2t as in the fully-interacting dimer.
   u = U/(4*t)
   # Impurity correlation energy
   a_12 = 0.5*(1 - rho)
   a_21 = 0.5*np.sqrt(rho*(1- rho)/2.0)
   a_11 = a_21*(1 + 1.0/rho)
   a_22 = a_12*0.5
   a_1 = a_11 + u*a_12
   a_2 = a_21 + u*a_22
   da12_drho = - 0.5
   da21_drho = (1 - 2*rho)/(8*np.sqrt((1 - rho)*rho/2))
   da11_drho = da21_drho*(1 + 1.0/rho) - a_21/rho**2
   da22_drho = - 0.25
   da1_drho = da11_drho + u*da12_drho
   da2_drho = da21_drho + u*da22_drho
   N = (1 - rho)*(1 + rho*(1 + (1 + rho)**3*u*a_1))
   D = 1 + (1 + rho)**3*u*a_2
   dN_drho = - 1 + (1 - 2*rho)*(1 + (1 + rho)**3*u*a_1) + rho*u*(1 - rho)*(1 + rho)**2*(3*a_1 + (1 + rho)*da1_drho)
   dD_drho = u*(1 + rho)**2*(3*a_2 + (1 + rho)*da2_drho)
   g0 = np.sqrt((1 - rho)*(1 + rho*(1 + (1 + rho)**3*u*a_1))/(1 + (1 + rho)**3*u*a_2))
   dg0_drho = (dN_drho - g0**2*dD_drho)/(2*g0*D)
   Y0 = np.sqrt(1 - g0**2 - rho**2)
   dh_dg0 = g0*(g0**4 + 3*g0**2*rho**2 + 2*rho**2*(rho**2 - 1 - Y0))/(2*(g0**2 + rho**2)**2*Y0)
   j = (1 - rho)*(1 + rho)**3*u**2
   k = (3*rho/2 - 1 + rho*(1 + rho)**3*u*a_2)*a_12 - rho*(1 + (1 + rho)**3*u*a_1)*a_22
   l = 2*g0*(1 + (1 + rho)**3*u*a_2)**2
   q = j*k/l
   g1 = g0 + (u*dh_dg0 - 1)*q
   Y1 = np.sqrt(1 - g1**2 - rho**2)
   h = (g1**2*(1 - Y1) + 2*rho**2)/(2*(g1**2 + rho**2))
   f = -2*t*g1 + 0.5*U*h
   Ts = - 2*t*np.sqrt(occ*(2 - occ))
   EHximp = U*0.5*(1 - occ*(1 - occ*0.5))
   Ecimp = f - Ts - EHximp

   # Derivative with respect to U
   da1_du = a_12
   da2_du = a_22
   v = 2*Y0*(g0**2 + rho**2)**2
   w = g0*(g0**4 + 3*g0**2*rho**2 + 2*rho**2*(rho**2 - 1 - Y0))
   dj_du = 2*(1 - rho)*(1 + rho)**3*u
   G = N/D
   dN_du = (1 - rho)*rho*(1 + rho)**3*(a_1 + u*da1_du)
   dD_du = (1 + rho)**3*(a_2 + u*da2_du)
   dG_du = (dN_du*D - N*dD_du)/D**2
   dg0_du = dG_du/(2*np.sqrt(G))
   dk_du = a_12*(rho*(1 + rho)**3*(a_2 + u*da2_du)) - a_22*(rho*(1 + rho)**3*(a_1 + u*da1_du))
   dl_du = 4*g0*(1 + (1 + rho)**3*u*a_2)*(1 + rho)**3*(a_2 + u*da2_du) + 2*dg0_du*(1 + (1 + rho)**3*u*a_2)**2
   dw_du = dg0_du*(g0**4 + 3*g0**2*rho**2 + 2*rho**2*(rho**2 - 1 - Y0) + g0*(4*g0**3 + 6*g0*rho**2 + 2*rho**2*g0/Y0))
   dv_du = g0*(g0**2 + rho**2)*dg0_du*(-2*(g0**2 + rho**2)/Y0 + 8*Y0)
   dh_dg1 = g1*(g1**4 + 3*g1**2*rho**2 + 2*rho**2*(rho**2 - 1 - Y1))/(2*(g1**2 + rho**2)**2*Y1)
   dq_du = ( (dj_du*k + j*dk_du)*l - j*k*dl_du )/l**2
   ddh_ddug0 = (dw_du*v - w*dv_du)/v**2
   dg1_du = dg0_du + (dh_dg0 + u*ddh_ddug0)*q + (u*dh_dg0 - 1)*dq_du
   dh_du = dg1_du*dh_dg1
   df_du = -2*t*(dg1_du - h - u*dh_du)
   dEcimp_dU = (1/(4*t))*df_du - 0.5*(1 + occ*(0.5*occ - 1))

   # Derivative with respect to t
   dEcimp_dt = Ecimp/t - U*dEcimp_dU/t   

   # Derivative with respect to n
   if occ_tmp == 1:
     dEcimp_dn = 0
   else:
     drho_dn = sign(occ - 1)
     dTs_dn = - 2*t*(1 - occ)/np.sqrt(occ*(2 - occ))
         # Note that P = j*k and Q = l in Eq. 155 of Ref. Senjean et al. Theoretical Chemistry Accounts (2018) 137:169
     dP_drho = (3*(1 - rho)*(1 + rho)**2 - (1 + rho)**3)*u**2*((3*rho/2 - 1 + rho*(1 + rho)**3*u*a_2)*a_12 \
               - rho*(1 + (1 + rho)**3*u*a_1)*a_22) \
               + (1 - rho)*(1 + rho)**3*u**2*((3/2.0 + 3*u*(1 + rho)**2*rho*a_2 \
               + u*(1 + rho)**3*(a_2 + rho*da2_drho))*a_12 \
               + (3*rho/2 - 1 + rho*(1 + rho)**3*u*a_2)*da12_drho - (rho*da22_drho + a_22)*(1 + (1 + rho)**3*u*a_1) \
               - rho*a_22*(3*(1 + rho)**2*u*a_1 + (1 + rho)**3*u*da1_drho))
     dQ_drho = 2*dg0_drho*(1 + (1 + rho)**3*u*a_2)**2 + 4*g0*(1 + (1 + rho)**3*u*a_2)*u*(3.0*(1 + rho)**2*a_2 + (1 + rho)**3*da2_drho)
     dq_drho = (dP_drho*l - j*k*dQ_drho)/l**2
     ddh_ddrhog0 = (- dg0_drho*(g0**2 + rho**2) + 4*g0*(g0*dg0_drho + rho))*(2*rho**2 + g0**2*(1 - Y0))/(g0**2 + rho**2)**3 \
                   - g0*(4*rho + 2*g0*dg0_drho*(1 - Y0) + g0**2*(g0*dg0_drho + rho)/Y0)/(g0**2 + rho**2)**2 \
                   - (g0*dg0_drho + rho)*(2*g0*(1 - Y0) + g0**3/Y0)/(g0**2 + rho**2)**2 \
                   + (2*dg0_drho*(1 - Y0) + 2*g0*rho/Y0 + 5*g0**2*dg0_drho/Y0 + g0**3*(g0*dg0_drho + rho)/Y0**3)/(2*(g0**2 + rho**2))
     dg1_drho = dg0_drho + (u*dh_dg0 - 1)*dq_drho + u*ddh_ddrhog0*q
     dh1_drho = (1/(2*(g1**2 + rho**2)))*(4*rho + 2*g1*dg1_drho*(1 - Y1) + g1**2*(g1*dg1_drho + rho)/Y1) \
               - (g1*dg1_drho + rho)*(2*rho**2 + g1**2*(1 - Y1))/(g1**2 + rho**2)**2
     df_drho = -2*t*dg1_drho + 0.5*U*dh1_drho
     dEHxcimp_dn = drho_dn*df_drho - dTs_dn + U/2.0
     dEcimp_dn = drho_dn*df_drho - dTs_dn + U/2.0 - U*occ*0.5
   return Ecimp, dEcimp_dU, dEcimp_dt, dEcimp_dn
