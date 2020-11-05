import thermo.unifac as th
import numpy as np
import scipy.optimize as sp
import math 
import matplotlib.pyplot as plt

chem_groups=[]
analyte={}
while True:
    try:
        subg=int(input('Subgroup: '))
        subg_num=int(input('N: '))
        analyte.update({subg:subg_num})
        del subg, subg_num
    except ValueError:
        break
chem_groups.append(analyte)
del analyte
chem_groups.append({16:1})
print(chem_groups)

temp=np.linspace(273.15,375.15,100)
gamma_list=np.linspace(0,1,100)
i=0
while i<len(temp):
    gamma=th.UNIFAC(T=temp[i], xs=[1e-15, 1], chemgroups=chem_groups)
    gamma_list[i]=gamma[0]
    del gamma
    i+=1
del i

popt, pcov = sp.curve_fit(lambda fx,a,b,c: a+(b/fx)+c*np.log(fx),  temp,  np.log(gamma_list), maxfev=1000)
print("ln(gamma)="+str(popt[0])+"+"+str(popt[1])+"/T + " + str(popt[2])+"ln(T)")
plt.style.use('ggplot')
plt.grid(False)
plt.scatter(temp, gamma_list)
plt.plot(temp, np.exp(popt[0]+popt[1]/temp+popt[2]*np.log(temp)))
plt.show()

print("Vm=18.02/(0.14395/(0.0112^(1+(1-temp/649.727)^0.05107)))")

i = input()

# print("gamma", gamma_list[0],"\n")
# density=0.14395/(0.0112**(1+((1-temp/649.727)**0.05107)))
# molar_volume_water=18.02/density
#
# print("log(p/kPa) = A1 - A2/(A3 + T/K)")
# A1=float(input("A1: "))
# A2=float(input("A2: "))
# A3=float(input("A3: "))
# sat_pressure=10**(A1-(A2/(A3+temp)))
# HLC_UNIFAC=sat_pressure*gamma_list*molar_volume_water/8.31446/temp
#
# print("\nHLC", HLC_UNIFAC[0], "\n")
#
# T_critical=float(input("Critical temperature (K): "))
# P_critical=float(input("Critical pressure (bar): "))
# T_boiling=float(input("Boiling temperature at np (K): "))
#
# Tbr=T_boiling/T_critical
# tau=1-Tbr
# f0=(-5.97616*tau+1.29874*(tau**1.5)-0.60394*(tau**2.5)-1.06841*(tau**5))/Tbr
# f1=(-5.03365*tau+1.11505*(tau**1.5)-5.41217*(tau**2.5)-7.46628*(tau**5))/Tbr
# W=-(f0+math.log(P_critical/1.01325, math.e))/f1
# del Tbr, tau, f0, f1
# d1=7.8149+(11.409*W)+(2.1674*(W**2))-(0.65342*(W**3))
# d2=0.81892-0.67637*W+1.2798*(W**2)-0.47594*(W**3)
# d3=-0.84408+1.8297*W-3.2435*(W**2)+1.449*(W**3)
# d4=0.41923-1.0892*W+1.9138*(W**2)-0.65758*(W**3)
#
# dHvap=[]
# i=0
# while i<len(temp):
#     term1=d4*((temp[i]/T_critical)**2)
#     term2=d3*(temp[i]/T_critical)
#     term3=1-(temp[i]/T_critical)
#     term4=term3**(d2+term2+term1)
#     dHvapi=term4*d1*8.31446*T_critical
#     dHvap.append(dHvapi)
#     del dHvapi, term1, term2, term3, term4
#     i+=1
# del d1, d2, d3, d4, i
#
# dHvap=np.array(dHvap)
#
# H_excess=-8.31*(temp**2)*(popt[1]*popt[0]*(temp**(popt[1]-1)))
# dH=dHvap-H_excess
#
# HLC_ref=float(input("HLC at 298.15 K (dimensionless): "))
#
# HLC_Morgan=HLC_ref*(998.946/density)*(298.15/temp)*np.exp((dH/8.31446)*(1/298.15-1/temp))

# print('available functional groups: \n   -CH3\n   >CH2\n   >CH-\n   >C<\n   -CH2-(ring)\n   >CH-(ring)\n   =CH-(ring)\n   =C<(ring)\n   -Cl\n   -OH\n   >CO')
# fgb=[]
# fgc=[]
# nfg=[]
# while 1==1:
#     fgname=input('functional group: ')
#     if fgname=='-CH3':
#         fgib=7775
#         fgic=-21.76
#     elif fgname=='>CH2':
#         fgib=1498.0
#         fgic=-3.79
#     elif fgname=='>CH-':
#         fgib=-5623
#         fgic=17.03
#     elif fgname=='>C<':
#         fgib=-16670
#         fgic=48.34
#     elif fgname=='-CH2- ring':
#         fgib=3436.99
#         fgic=-8.95
#     elif fgname=='>CH- ring':
#         fgib=-5298.99
#         fgic=14.32
#     elif fgname=='=CH- ring':
#         fgib=2179.99
#         fgic=-5.24
#     elif fgname=='=C< ring':
#         fgib=-3902.99
#         fgic=12.29
#     elif fgname=='-Cl':
#         fgib=4792
#         fgic=-10.65
#     elif fgname=='-OH':
#         fgib=11681
#         fgic=-21.5
#     elif fgname=='>CO':
#         fgib=-3310
#         fgic=19.5
#     elif fgname=='end':
#         break
#     else:
#         continue
#     fgb.append(fgib)
#     fgc.append(fgic)
#     nfgi=int(input('number: '))
#     nfg.append(nfgi)
# del fgib, fgic, fgname
# nfga=np.array(nfg)
# fgba=np.array(fgb)
# fgca=np.array(fgc)
# del nfg, fgb, fgc
# B=np.sum(fgba*nfga)
# C=np.sum(fgca*nfga)
# print('B:', B)
# print('C:', C)
# HLC_Lau=HLC_ref*(998.946/density)*((temp/298.15)**(C-1))*np.exp(B*(1/298.15-1/temp))
# plt.style.use("ggplot")
# plt.grid(False)
# plt.plot(temp, HLC_UNIFAC, label="UNIFAC")
# plt.plot(temp, HLC_Morgan, label="Морган&UNIFAC")
# plt.plot(temp, HLC_Lau, label="Лау")

# plt.xlabel("Температура, К")
# plt.ylabel("Константа Генри")
# plt.legend()
# an = input("analyte: ")
# plt.savefig(an+".png", dpi = 400)