import scipy.constants as scc
import numpy as np
import matplotlib.pyplot as plt

lande_factor_ground_state=2.0032
lande_factor_excited_state=1.335
target_center_z=0.5
bohr_magnetron = scc.physical_constants['Bohr magneton'][0]
hbar=scc.hbar
ERyd=scc.physical_constants['Rydberg constant times hc in J'][0]
A_factor_gs = 152e6*hbar#221.864e6*hbar*2*np.pi#152.136e6*hbar#(1 / 137)**2 * 1**4 * ERyd / (1 * (1 + 0.5) * (1 + 1) * 2**3) #1.1e6*hbar
A_factor_es = -1.15e6*hbar
mu_N=scc.physical_constants["nuclear magneton"][0]
g_I=-0.0004476540

g_J_gs=lande_factor_ground_state
#mu_prime_gs=mu_N*g_I/g_J_gs
mu_min_gs=1.0007*bohr_magnetron#bohr_magnetron-mu_prime_gs
mu_plus_gs=1.0016*bohr_magnetron#bohr_magnetron+mu_prime_gs

g_J_es=lande_factor_excited_state
mu_prime_es=mu_N*g_I/g_J_es
mu_min_es=bohr_magnetron-mu_prime_es
mu_plus_es=bohr_magnetron+mu_prime_es

E_a=[] #E_a: F=3/2; mf=3/2
E_b=[] #E_b: F=3/2; mf=1/2
E_c=[] #E_c: F=3/2; mf=-1/2
E_d=[] #E_d: F=3/2; mf=-3/2
E_e=[] #E_e: F=1/2; mf=1/2
E_f=[] #E_f: F=1/2; mf=-1/2
E_ugs=[]
E_lgs_a=[]
E_lgs_b=[]
E_ugs_a=[]
E_ugs_b=[]
E_ugs_c=[]
E_ugs_d=[]

B_field_array=np.arange(0,0.017,0.001) # 0 bis 0.1T entspricht 0 bis 1000G
for B_field in B_field_array:
    #E_a.append(A_factor / 2 + 2 * bohr_magnetron * B_field)
    #E_b.append(A_factor / 2 * (
     #       bohr_magnetron * B_field / A_factor - 0.5 + np.sqrt(
     #   (bohr_magnetron * B_field / A_factor + 0.5) ** 2 + 2)))
    #E_c.append(A_factor / 2 * (
    #        -bohr_magnetron * B_field / A_factor - 0.5 + np.sqrt(
    #    (bohr_magnetron * B_field / A_factor - 0.5) ** 2 + 2)))
    #E_d.append(A_factor / 2 - 2 * bohr_magnetron * B_field)
    #E_e.append(A_factor / 2 * (
    #        bohr_magnetron * B_field / A_factor - 0.5 - np.sqrt(
    #    (bohr_magnetron * B_field / A_factor + 0.5) ** 2 + 2)))
    #E_f.append(A_factor / 2 * (
    #        -bohr_magnetron * B_field / A_factor - 0.5 - np.sqrt(
    #    (bohr_magnetron * B_field / A_factor - 0.5) ** 2 + 2)))

    E_f.append(1/4*(2*B_field*mu_prime_es-A_factor_es-np.sqrt(9*A_factor_es**2-4*A_factor_es*B_field*(mu_plus_es+mu_prime_es)+4*B_field**2*(mu_plus_es+mu_prime_es)**2))) #E_a: F=3/2; mf=3/2
    E_e.append(1/4*(-2*B_field*mu_prime_es-A_factor_es-np.sqrt(9*A_factor_es**2+4*A_factor_es*B_field*(mu_plus_es+mu_prime_es)+4*B_field**2*(mu_plus_es+mu_prime_es)**2))) #E_b: F=3/2; mf=1/2
    E_d.append(1/2*A_factor_es-B_field*mu_min_es)#E_c: F=3/2; mf=-1/2
    E_c.append(1/4*(2*B_field*mu_prime_es-A_factor_es+np.sqrt(9*A_factor_es**2-4*A_factor_es*B_field*(mu_plus_es+mu_prime_es)+4*B_field**2*(mu_plus_es+mu_prime_es)**2)))#E_d: F=3/2; mf=-3/2
    E_b.append(1/4*(-2*B_field*mu_prime_es-A_factor_es+np.sqrt(9*A_factor_es**2+4*A_factor_es*B_field*(mu_plus_es+mu_prime_es)+4*B_field**2*(mu_plus_es+mu_prime_es)**2)))#E_e: F=1/2; mf=1/2
    E_a.append(1/2*A_factor_es+B_field*mu_min_es)#E_f: F=1/2; mf=-1/2



    E_lgs_a.append(1/4*(-2*(mu_plus_gs-bohr_magnetron)*B_field-A_factor_gs-np.sqrt(9*A_factor_gs**2+4*A_factor_gs*B_field*\
                                                (mu_plus_gs+bohr_magnetron)+4*B_field**2*(mu_plus_gs+bohr_magnetron)**2) ) ) #mf -1/2
    E_lgs_b.append(1/4*(-2*(mu_plus_gs-bohr_magnetron)*B_field-A_factor_gs+np.sqrt(9*A_factor_gs**2+4*A_factor_gs*B_field*\
                                                (mu_plus_gs+bohr_magnetron)+4*B_field**2*(mu_plus_gs+bohr_magnetron)**2) ) ) #mf 1/2

    E_ugs_a.append(1/2*A_factor_gs+B_field*mu_min_gs) #mf -3/2
    E_ugs_b.append(1/4*(2*B_field*(mu_plus_gs-bohr_magnetron)-A_factor_gs-np.sqrt(9*A_factor_gs**2-4*A_factor_gs*B_field*\
                                                        (mu_plus_gs+bohr_magnetron)+4*B_field**2*(mu_plus_gs+bohr_magnetron)**2) ) ) #mf -1/2
    E_ugs_c.append(1/4*(2*B_field*(mu_plus_gs-bohr_magnetron)-A_factor_gs+np.sqrt(9*A_factor_gs**2-4*A_factor_gs*B_field*\
                                                        (mu_plus_gs+bohr_magnetron)+4*B_field**2*(mu_plus_gs+bohr_magnetron)**2) ) ) #mf 1/2
    E_ugs_d.append(1/2*A_factor_gs-B_field*mu_min_gs) #mf 3/2
    #E_lgs_a.append( (bohr_magnetron * lande_factor_ground_state * -1/2 * B_field) + (-152e6*hbar) )  # -1/2
    #E_lgs_b.append( (bohr_magnetron * lande_factor_ground_state * 1/2 * B_field) + (-152e6*hbar) )  # 1/2
    #E_ugs_a.append( (bohr_magnetron * lande_factor_ground_state * -3/2 * B_field) +((152e6)/2)*hbar )  #-3/2
    #E_ugs_b.append( (bohr_magnetron * lande_factor_ground_state * -1/2 * B_field) +((152e6)/2)*hbar )  #-1/2
    #E_ugs_c.append( (bohr_magnetron * lande_factor_ground_state * 1/2 * B_field) +((152e6)/2)*hbar )   #1/2
    #E_ugs_d.append( (bohr_magnetron * lande_factor_ground_state * 3/2 * B_field) +((152e6)/2)*hbar )   #3/2

#E_state= 0#-ERyd*3**2/2**2
#E_FS_exc=0#E_state + 3**2*(1/137)**2/2 * (1/(3/2+1/2)-3/(4*2))
#E_FS_gs=0#E_state + 3**2*(1/137)**2/2 * (1/(1/2+1/2)-3/(4*2))
GHz=1e-9
MHz=1e-6
scale=MHz
omega_a=[scale*energy/hbar for energy in E_a]
omega_b=[scale*energy/hbar for energy in E_b]
omega_c=[scale*energy/hbar for energy in E_c]
omega_d=[scale*energy/hbar for energy in E_d]
omega_e=[scale*energy/hbar for energy in E_e]
omega_f=[scale*energy/hbar for energy in E_f]

omega_lgs_a=[scale*energy/hbar for energy in E_lgs_a]
omega_lgs_b=[scale*energy/hbar for energy in E_lgs_b]
omega_ugs_a=[scale*energy/hbar for energy in E_ugs_a]
omega_ugs_b=[scale*energy/hbar for energy in E_ugs_b]
omega_ugs_c=[scale*energy/hbar for energy in E_ugs_c]
omega_ugs_d=[scale*energy/hbar for energy in E_ugs_d]

#Plotting
#fig, ((ax1), (ax2)) = plt.subplots(2, 1, sharex=True, sharey=False)

#ax1.plot(B_field_array*1e4,omega_a,label="e.s.: F=3/2; mf=3/2", color="red")
#ax1.plot(B_field_array*1e4,omega_b,label="e.s.: F=3/2; mf=1/2", color="red")
#ax1.plot(B_field_array*1e4,omega_c,label="e.s.: F=3/2; mf=-1/2", color="red")
#ax1.plot(B_field_array*1e4,omega_d,label="e.s.: F=3/2; mf=-3/2", color="red")
#ax1.plot(B_field_array*1e4,omega_e,label="e.s.: F=1/2; mf=1/2", color="black")
#ax1.plot(B_field_array*1e4,omega_f,label="e.s.: F=1/2; mf=-1/2", color="black")
#ax1.legend()
#ax1.grid()

ax2=plt
ax2.plot(B_field_array*1e4,omega_ugs_a,label="ugs.: F=3/2; mf=3/2", color="red")
ax2.plot(B_field_array*1e4,omega_ugs_c,label="ugs.: F=3/2; mf=1/2", color="blue")
ax2.plot(B_field_array*1e4,omega_lgs_b,label="ugs.: F=3/2; mf=-1/2", color="green")
ax2.plot(B_field_array*1e4,omega_ugs_d,label="ugs.: F=3/2; mf=-3/2", color="orange")

ax2.plot(B_field_array*1e4,omega_ugs_b,label="lgs.: F=1/2; mf=-1/2", color="black")
ax2.plot(B_field_array*1e4,omega_lgs_a,label="lgs.: F=1/2; mf=1/2", color="magenta")
ax2.legend()
ax2.grid()

plt.xlabel("B/ G")
if scale==1e-9: Hz="GHz"
else: Hz="MHz"
plt.ylabel("energy shift/ {}".format(Hz))
#fig.text(0.06, 0.5, "energy shift/ {}".format(Hz), ha='center', va='center', rotation='vertical')

plt.show()