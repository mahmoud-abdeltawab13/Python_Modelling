import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *
import sympy as sp
import scipy
import Optimization1
#Paramater Defination


# def Plot(Gain):
##Inputs From User
Fref=160E6
PM_deg=60
Kvco=90E6
Fout=156.25E6
P_divider=24
Gain=105E-6

#################################################
#Automating Computation of X
PM_rad=np.pi/180 *float(PM_deg) #take as input from user 
val=sp.tan(PM_rad) #represent the value that will be used in eqn to compute X
y=sp.Symbol('y')
soln=sp.solve(sp.Eq(0.5*y-0.5*1/y,val),y)              
for i in range (1,len(soln)):
    if(soln[i]>0):
        z=soln[i]
# ####################################  

# #Calculating N_divider
N=(Fout/Fref)*P_divider


# ##MODELED PHASE NOISE

########################
#Modeled PN of Crystal Refrence
beta1=6.22036E-10
beta2=0.0248/(2*np.pi*1E3)
beta3=0.0248
NoiseFloor_REF=1E-16


num4=np.array([NoiseFloor_REF,beta1,beta2,beta3])
den4=np.array([1,0,0,0])
PN_REF=tf(num4,den4)
print(PN_REF)
mag_Ref,phase,w=bode(PN_REF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])

plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_Ref),label="Reference")
plt.grid(True,which="both",axis="both")
plt.title("Refrence Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")


# ###########
#Modeled PN of VCO
num5=np.array([124.685,39305502.88])
den5=np.array([1,0,0,0])
PN_VCO=tf(num5,den5)
print(PN_VCO)
mag_VCO,phase,w=bode(PN_VCO,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])

for i in range(len(w)):
    if 10*np.log10(mag_VCO[i])<= -145:
        mag_VCO[i:]=mag_VCO[i]

plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_VCO),label="vco")
plt.grid(True,which="both",axis="both")
plt.title("VCO Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
#####################

#Modeled PN of CP
alpha=(7.943E-24)*((300E-6)**2)

num6=np.array([alpha,alpha*(2*np.pi*2E6)])
den6=np.array([1,0])
PN_CP=tf(num6,den6)
print(PN_CP)
mag_CP,phase,w=bode(PN_CP,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])

#####################

#Modeled PN of Dividers and O/P Buffer
alpha1=10**(-160/10)#1.2566E-9
num7=np.array([alpha1,alpha1*(2*np.pi*2E6)])
den7=np.array([1,0])
PN_Dividers=tf(num7,den7)
print(PN_Dividers)
mag_Dividers,phase,w=bode(PN_Dividers,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_Dividers),label="Dividers")
plt.grid(True,which="both",axis="both")
plt.title("Dividers & O/P Buffer Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()
#####################

#DSM PHASE_NOISE
Sq=1/(12*Fref)
f=w/(2*np.pi)#np.logspace(3,7.60206,1000)
NTF=np.abs((2*np.sin(np.pi*f/Fref))**4)
Sg=Sq*NTF
S_phi_n=Sg*((Fref/f)**2)/(N**2)
plt.figure(figsize=(10, 8))
plt.semilogx(f,10*np.log10(S_phi_n),label="Reference")
plt.grid(True,which="both",axis="both")
plt.title("DSM PHASE NOISE")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")

#######################
#Modeled PN of VCXO I/P Referred Noise
alpha2=10**(-135/10)
num9=np.array([alpha2,alpha2*(2*np.pi*1E4)])
den9=np.array([1,0])
PN_Input_Referred=tf(num9,den9)
print(PN_Input_Referred)
mag_Input_Referred_Noise,phase,w=bode(PN_Input_Referred,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])

plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_Input_Referred_Noise),label="I/P Referred Noise")
plt.grid(True,which="both",axis="both")
plt.title("VCXO I/P Referred Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()


#Plot of the LPF TF

#Defining LPF TF H_LPF(s)
wp_lpf=2*np.pi*1E4
num10=np.array([1])
den10=np.array([1/wp_lpf,1])
H_LPF=tf(num10,den10)
print("H_LPF(S)=",H_LPF)
########################
mag_TF_LPF,phase,w=bode(H_LPF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_LPF**2))
plt.grid(True,which="both",axis="both")
plt.title("LPF Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
num=[2*np.pi*Fref/N]
den=[1,0]
integr=tf(num,den)
print("Integrator= ",integr)
mag_TF_integr,phase,w=bode(integr,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])

mag_VCXO=((Gain*mag_TF_LPF)**2)*(mag_Input_Referred_Noise)*(mag_TF_integr**2)
plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_VCXO))
plt.grid(True,which="both",axis="both")
plt.title("VCXO Magnitude Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")


####################



#####Getting Optimum Omega_u
PN_CLOSEIN=(mag_VCXO+mag_Ref+mag_Dividers+S_phi_n+(mag_CP*((2*np.pi)**2)))*(N**2)

DIFF=np.abs(PN_CLOSEIN-mag_VCO)
fu_inital=int(w[np.where(DIFF==np.min(DIFF))])/(2*np.pi) 
print("fu_initial= ",fu_inital)
plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(PN_CLOSEIN),label="CLOSEIN")
plt.semilogx(w/(2*np.pi),10*np.log10(mag_VCO),label="VCO")
plt.grid(True,which="both",axis="both")
plt.title("CLOSEIN & VCO ")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()



plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(DIFF),label="CloseIN-VCO")

plt.grid(True,which="both",axis="both")
plt.title("CLOSEIN-VCO ")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()


fu,Ip=Optimization1.Omega_Optimum(fu_inital,Fref,PM_deg,Kvco,Fout,P_divider,N,Gain) #Optimum fu
print("optimum fu=",fu*1E-6,"MHz")
print("optimum Ip=",Ip*1E6)
#######################################

#Calculating Caps and Resistors
x=float(z)
omega_u=2*np.pi*fu
omega_z=omega_u/x
omega_p=omega_u*x
K=Ip*Kvco
Beta=1/N



r1 = omega_u*N*(1+1/((x**2)-1))/K
print("R1=",r1*1E-3,"Kohm")
c1=1/(omega_z*r1)
print("C1=",c1*1E12,"pF")
c2 = c1 /((x**2) - 1)
print("C2=",c2*1E12,"pF")
c_sum=c1+c2
print("Csum=",c_sum*1E12,"pF")




#Modeled PN of CP
alpha=(7.943E-24)*((300E-6)**2)/(Ip**2)

num6=np.array([alpha,alpha*(2*np.pi*2E6)])
den6=np.array([1,0])
PN_CP=tf(num6,den6)
print(PN_CP)
mag_CP,phase,w=bode(PN_CP,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_CP),label="Charge Pump")
plt.grid(True,which="both",axis="both")
plt.title("CP Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
##########################3
#Modeled PN of Loop Filter
Vn=4*(1.38E-23)*300*r1
num8=np.array([Vn])
den8=np.array([1])
PN_LF=tf(num8,den8)
print(PN_LF)
mag_LF,phase,w=bode(PN_LF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_LF),label="loop filter")
plt.grid(True,which="both",axis="both")
plt.title("Loop filter Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
#####################


#DEFINING TF OF THE BLOCKS
#########################################
#Defining Transfer Functions of Z(s) & A(s)
num1=np.array([1])
den1=np.array([c_sum,0])
H1=tf(num1,den1)

num2=np.array([1/omega_z,1])
den2=np.array([1/omega_p,1])
H2=tf(num2,den2)

#Defining Z(s)
Z=series(H1,H2)
#######
print("Z=",Z)
##################



###
#Defining Open loop TF A(s)
num3=np.array([K])
den3=np.array([1,0])
H3=tf(num3,den3)

A=series(Z,H3)
print("A(S)=",A)
#################

###
#Defining loop TF T(s)
T=A*Beta
print("T(S)=",T)
#################



####
#Defining Close_loop TF H(s)
H=A/(1+T)
print("H(S)=",H)
################

###################################################################################



#Defining Reference TF H_REF(s)
H_REF=H/(P_divider)
print("H_REF(S)=",H_REF)
###########



#Defining VCO TF H_VCO(s)
H_VCO=1/(P_divider*(1+T))
print("H_VCP(S)=",H_VCO)
###########

#Defining CP TF H_CP(s)
H_CP=2*H*np.pi/(Ip*P_divider)
print("H_CP(S)=",H_CP)
###########


#Defining LF TF H_LF(s)

s=TransferFunction.s
H_LF=2*np.pi*Kvco/(P_divider*s*(1+T))
print("H_LF(S)=",H_LF)
###########


#Defining N-Divider TF H_N_Divider(s)
H_N_Divider=H/(P_divider)
print("H_N_Divider(S)=",H_N_Divider)
###########


#Defining DSM TF H_DSM(s)
H_DSM=H/(P_divider)
print("H_DSM(S)=",H_DSM)
###########


#Defining O/P Buffer & P-divider TF H_Buffer(s)
H_Buffer=1
print("H_Buffer(S)=",H_Buffer)
############################################


#####################
##Closed_loop Response
plt.figure(figsize=(10, 8))
bode(H,dB=True,deg=True,plot=True,Hz=True)
plt.title("CLOSED LOOP Response")
plt.grid(True,which="both",axis="both")


# #####################################
##LOOP Response
plt.figure(figsize=(10, 8))
bode(T,dB=True,deg=True,plot=True,margins=True,grid=True,Hz=True)
plt.title("LOOP RESPONSE")
plt.grid(True,which="both",axis="both")
gm,pm,w180,wc=margin(T)
GX=(wc/(2*np.pi))*10**-6
PX=w180/(2*np.pi)
gm_db=mag2db(gm)
print("GX=",f'{GX:0.3f}',"MHz")
print("PM=",f'{pm:0.3f}',"deg")
# #####################################



# ##PLOT OF TRANSFER FUNCTION OF EACH BLOCK
# ##########################################

#Plot of the TF seen by the Reference
mag_TF_REF,phase,w=bode(H_REF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_REF**2))
plt.grid(True,which="both",axis="both")
plt.title("Reference Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################


#Plot of the TF seen by the VCO
mag_TF_VCO,phase,w=bode(H_VCO,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_VCO**2))
plt.grid(True,which="both",axis="both")
plt.title("VCO Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################


#Plot of the TF seen by the Charge Pump (CP)
mag_TF_CP,phase,w=bode(H_CP,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_CP**2))
plt.grid(True,which="both",axis="both")
plt.title("Charge Pump (CP) Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################

#Plot of the TF seen by the Loop Filter (LF)
mag_TF_LF,phase,w=bode(H_LF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_LF**2))
plt.grid(True,which="both",axis="both")
plt.title("Loop Filter (LF) Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################


#####################
#Plot of the TF seen by the N Divider
mag_TF_R_Divider,phase,w=bode(H_N_Divider,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_R_Divider**2))
plt.grid(True,which="both",axis="both")
plt.title(" N-Divider Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################



#Plot of the TF seen by the O/P Buffer & P-Divider
mag_TF_Buffer_dB=10*np.log10(H_Buffer**2)*(w/w)


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),mag_TF_Buffer_dB)
plt.grid(True,which="both",axis="both")
plt.title("O/P Buffer & P-Divider Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################


#Plot of the TF seen by the DSM
mag_TF_DSM,phase,w=bode(H_DSM,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_DSM**2))
plt.grid(True,which="both",axis="both")
plt.title("DSM Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")

################################



########## PN OF EACH BLOCK @OUTPUT

#Plot of the PN of the Crystal Reference @Output
PN_REF_OUT_mag=(mag_TF_REF**2)*(mag_Ref)

PN_REF_OUT=10*np.log10(mag_TF_REF**2)+10*np.log10(mag_Ref)

plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_REF_OUT,label="REF")
plt.grid(True,which="both",axis="both")
plt.title("Reference Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()

#####################



#Plot of the PN of the VCO @Output
PN_VCO_OUT_mag=(mag_TF_VCO**2)*(mag_VCO)

PN_VCO_OUT=10*np.log10(mag_TF_VCO**2)+10*np.log10(mag_VCO)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_VCO_OUT,label="VCO")
plt.grid(True,which="both",axis="both")
plt.title("VCO Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()

#####################


#Plot of the PN of the CP @Output
PN_CP_OUT_mag=(mag_TF_CP**2)*(mag_CP)

PN_CP_OUT=10*np.log10(mag_TF_CP**2)+10*np.log10(mag_CP)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_CP_OUT,label="CP")
plt.grid(True,which="both",axis="both")
plt.title("Charge Pump(CP) Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()


#####################


#Plot of the PN of the LF @Output
PN_LF_OUT_mag=(mag_TF_LF**2)*(mag_LF)

PN_LF_OUT=10*np.log10(mag_TF_LF**2)+10*np.log10(mag_LF)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_LF_OUT,label="LF")
plt.grid(True,which="both",axis="both")
plt.title(" Loop Filter(LF) Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()

#####################


#Plot of the PN of the N-Divider @Output
PN_N_Divider_OUT_mag=(mag_TF_R_Divider**2)*(mag_Dividers)

PN_N_Divider_OUT=10*np.log10(mag_TF_R_Divider**2)+10*np.log10(mag_Dividers)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_N_Divider_OUT,label="N-DIV")
plt.grid(True,which="both",axis="both")
plt.title(" N-Divider Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()

#####################

#Plot of the PN of the P-Divider @Output
PN_P_Divider_OUT_mag=(10**(mag_TF_Buffer_dB/10))*(mag_Dividers)

PN_P_Divider_OUT=mag_TF_Buffer_dB+10*np.log10(mag_Dividers)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_P_Divider_OUT,label="P-DIV")
plt.grid(True,which="both",axis="both")
plt.title(" P-Divider Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()

#####################


#Plot of the PN of the Output Buffer @Output
PN_Buffer_OUT_mag=(10**(mag_TF_Buffer_dB/10))*(mag_Dividers)

PN_Buffer_OUT=mag_TF_Buffer_dB+10*np.log10(mag_Dividers)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_Buffer_OUT,label="O/P Buffer")
plt.grid(True,which="both",axis="both")
plt.title(" Buffer Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()

###########################

#Plot of the PN of the DSM @Output
PN_DSM_OUT_mag=(mag_TF_DSM**2)*(S_phi_n)

PN_DSM_OUT=10*np.log10(mag_TF_DSM**2)+10*np.log10(S_phi_n)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_DSM_OUT,label="DSM")
plt.grid(True,which="both",axis="both")
plt.title(" DSM Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()
# #####################


#Plot of the VCXO I/P Referred Noise @Output
PN_VCXO_Noise_OUT_mag=mag_VCXO*(mag_TF_DSM**2)

PN_VCXO_Noise_OUT=10*np.log10(mag_VCXO)+10*np.log10(mag_TF_DSM**2)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_VCXO_Noise_OUT,label="VCXO")
plt.grid(True,which="both",axis="both")
plt.title(" VCXO Output Phase Noise")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()

##########################
#TOTAL OUTPUT PHASE NOISE
Total_Out_PN_mag=PN_VCXO_Noise_OUT_mag+PN_DSM_OUT_mag+PN_Buffer_OUT_mag+PN_P_Divider_OUT_mag+PN_N_Divider_OUT_mag+PN_CP_OUT_mag+PN_LF_OUT_mag+PN_VCO_OUT_mag+PN_REF_OUT_mag
Total_Out_PN=10*np.log10(Total_Out_PN_mag)

plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),Total_Out_PN,label="Total")
plt.grid(True,which="both",axis="both")
plt.title("Phase Noise @ The Output")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")

plt.legend()

# ####################################################

# ####JITTER CALCULATION
# #TOTAL RMS JITTER
f1=np.logspace(np.log10(12E3), np.log10(20E6),1000)
Total_Integral=np.trapz(10**(Total_Out_PN/10),x=f1)
Total_Integral=2*Total_Integral
Total_RMS_Jitter=np.sqrt(Total_Integral)/(2*np.pi*Fout)
print("Total RMS Jitter= ",Total_RMS_Jitter*1E15,"fsec")


#REF RMS JITTER
REF_Integral=np.trapz(10**(PN_REF_OUT/10),x=f1)
REF_Integral=2*REF_Integral
REF_RMS_Jitter=np.sqrt(REF_Integral)/(2*np.pi*Fout)
print("Reference RMS Jitter= ",REF_RMS_Jitter*1E15,"fsec")


#VCO RMS JITTER
VCO_Integral=np.trapz(10**(PN_VCO_OUT/10),x=f1)
VCO_Integral=2*VCO_Integral
VCO_RMS_Jitter=np.sqrt(VCO_Integral)/(2*np.pi*Fout)
print("VCO RMS Jitter= ",VCO_RMS_Jitter*1E15,"fsec")



#CP RMS JITTER
CP_Integral=np.trapz(10**(PN_CP_OUT/10),x=f1)
CP_Integral=2*CP_Integral
CP_RMS_Jitter=np.sqrt(CP_Integral)/(2*np.pi*Fout)
print("Charge Pump(CP) RMS Jitter= ",CP_RMS_Jitter*1E15,"fsec")



#LF RMS JITTER
LF_Integral=np.trapz(10**(PN_LF_OUT/10),x=f1)
LF_Integral=2*LF_Integral
LF_RMS_Jitter=np.sqrt(LF_Integral)/(2*np.pi*Fout)
print("Loop Filter(LF) RMS Jitter= ",LF_RMS_Jitter*1E15,"fsec")


#N-Divider RMS JITTER
N_Divider_Integral=np.trapz(10**(PN_N_Divider_OUT/10),x=f1)
N_Divider_Integral=2*N_Divider_Integral
N_Divider_RMS_Jitter=np.sqrt(N_Divider_Integral)/(2*np.pi*Fout)
print("N-Divider RMS Jitter= ",N_Divider_RMS_Jitter*1E15,"fsec")


#P-Divider RMS JITTER
P_Divider_Integral=np.trapz(10**(PN_P_Divider_OUT/10),x=f1)
P_Divider_Integral=2*P_Divider_Integral
P_Divider_RMS_Jitter=np.sqrt(P_Divider_Integral)/(2*np.pi*Fout)
print("P-Divider RMS Jitter= ",P_Divider_RMS_Jitter*1E15,"fsec")


#Buffer RMS JITTER
Buffer_Integral=np.trapz(10**(PN_Buffer_OUT/10),x=f1)
Buffer_Integral=2*Buffer_Integral
Buffer_RMS_Jitter=np.sqrt(Buffer_Integral)/(2*np.pi*Fout)
print("Buffer RMS Jitter= ",Buffer_RMS_Jitter*1E15,"fsec")


#DSM RMS JITTER
DSM_Integral=np.trapz(10**(PN_DSM_OUT/10),x=f1)
DSM_Integral=2*DSM_Integral
DSM_RMS_Jitter=np.sqrt(DSM_Integral)/(2*np.pi*Fout)
print("DSM RMS Jitter= ",DSM_RMS_Jitter*1E15,"fsec")


#VCXO RMS JITTER
VCXO_Integral=np.trapz(10**(PN_VCXO_Noise_OUT/10),x=f1)
VCXO_Integral=2*VCXO_Integral
VCXO_RMS_Jitter=np.sqrt(VCXO_Integral)/(2*np.pi*Fout)
print("VCXO RMS Jitter= ",VCXO_RMS_Jitter*1E15,"fsec")
########################################################

#Contribution Percentage of each block in the total output integrated RMS jitter

#REF Contribution
REF_Contribution_Precentage=(REF_Integral/Total_Integral)*100
print("Reference Contribution Precentage= ",REF_Contribution_Precentage,"%")


#VCO Contribution
VCO_Contribution_Precentage=(VCO_Integral/Total_Integral)*100
print("VCO Contribution Precentage= ",VCO_Contribution_Precentage,"%")

#CP Contribution
CP_Contribution_Precentage=(CP_Integral/Total_Integral)*100
print("Charge Pump (CP) Contribution Precentage= ",CP_Contribution_Precentage,"%")


#LF Contribution
LF_Contribution_Precentage=(LF_Integral/Total_Integral)*100
print("Loop Filter (LF) Contribution Precentage= ",LF_Contribution_Precentage,"%")



#N-Divider Contribution
N_Divider_Contribution_Precentage=(N_Divider_Integral/Total_Integral)*100
print("N-Divider Contribution Precentage= ",N_Divider_Contribution_Precentage,"%")


#P-Divider Contribution
P_Divider_Contribution_Precentage=(P_Divider_Integral/Total_Integral)*100
print("P-Divider Contribution Precentage= ",P_Divider_Contribution_Precentage,"%")

#Buffer Contribution
Buffer_Contribution_Precentage=(Buffer_Integral/Total_Integral)*100
print("Buffer Contribution Precentage= ",Buffer_Contribution_Precentage,"%")

#DSM Contribution
DSM_Contribution_Precentage=(DSM_Integral/Total_Integral)*100
print("DSM Contribution Precentage= ",DSM_Contribution_Precentage,"%")

#VCXO Contribution
VCXO_Contribution_Precentage=(VCXO_Integral/Total_Integral)*100
print("VCXO Contribution Precentage= ",VCXO_Contribution_Precentage,"%")

plt.show()





