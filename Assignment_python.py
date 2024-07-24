import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *
import sympy as sp
import scipy
import Optimization
#Paramater Defination

##Inputs From User
Fref=1E6*float(input("Enter Reference Freq(in MHz):"))
PM_deg=input("Enter the Phase Margin required(in deg):")
Ip=1E-6*float(input("Enter the Charge pump current(Ip)(in microAmp):"))
Kvco=1E6*float(input("Enter the VCO Gain(Kvco)(in MHz/V):"))
Fout=1E6*float(input("Enter Required Output Freq(in MHz):"))

#P-divider value based on the Fout

if Fout==560E6:
    P_divider=6
elif Fout==140E6:
    P_divider=24
else:
    print("Please eneter Valid value for Fout(either 560 OR 140)")
#################################################

#Automating Computation of X
PM_rad=np.pi/180 *float(PM_deg) #take as input from user 
val=sp.tan(PM_rad) #represent the value that will be used in eqn to compute X
y=sp.Symbol('y')
soln=sp.solve(sp.Eq(0.5*y-0.5*1/y,val),y)              
for i in range (1,len(soln)):
    if(soln[i]>0):
        z=soln[i]
        
####################################  

#Calculating N_divider
R_divider=2    
N=(Fout/Fref)*R_divider*P_divider


##MODELED PHASE NOISE

########################
#Modeled PN of Refrence
Noise_floor_Ref=1E-15
num4=np.array([4.95458E-3,246.089])
den4=np.array([1,0,0,0])
PN_REF=tf(num4,den4)
print(PN_REF)
mag_Ref,phase,w=bode(PN_REF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])
mag_Ref=mag_Ref+Noise_floor_Ref

plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_Ref),label="Reference")
plt.grid(True,which="both",axis="both")
plt.title("Refrence Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")

###########

#Modeled PN of VCO
Noise_floor_VCO=pow(10,-14.5)
num5=np.array([124.685,39305502.88])
den5=np.array([1,0,0,0])
PN_VCO=tf(num5,den5)
print(PN_VCO)
mag_VCO,phase,w=bode(PN_VCO,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])
mag_VCO=mag_VCO+Noise_floor_VCO


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_VCO),label="vco")
plt.grid(True,which="both",axis="both")
plt.title("VCO Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
#####################

#Modeled PN of CP
alpha=(2.6476E-20)*(Ip)

num6=np.array([alpha,alpha*(2*np.pi*2E6)])
den6=np.array([1,0])
PN_CP=tf(num6,den6)
print(PN_CP)
mag_CP,phase,w=bode(PN_CP,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_CP),label="Charge Pump")
plt.grid(True,which="both",axis="both")
plt.title("CP Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
#####################

#Modeled PN of Dividers and O/P Buffer
alpha1=1.2566E-9
num7=np.array([alpha1/(2*np.pi*2E6),alpha1])
den7=np.array([1,0])
PN_Dividers=tf(num7,den7)
print(PN_Dividers)
mag_Dividers,phase,w=bode(PN_Dividers,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_Dividers),label="Dividers")
plt.grid(True,which="both",axis="both")
plt.title("Dividers & O/P Buffer Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
plt.legend()
#####################




#####Getting Optimum Omega_u
PN_CLOSEIN=((mag_Ref/(R_divider**2))+mag_Dividers+mag_Dividers+(mag_CP*((2*np.pi/Ip)**2)))*(N**2)
DIFF=np.abs(PN_CLOSEIN-mag_VCO)
fu_inital=int(w[np.where(DIFF==np.min(DIFF))])/(2*np.pi) 
print("fu_initial= ",fu_inital*1E-6,'MHz')
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


fu=Optimization.Omega_Optimum(fu_inital,Fref,PM_deg,Ip,Kvco,Fout,R_divider,P_divider,N)  #Optimum fu
print("optimum fu=",fu*1E-6,"MHz")
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
print("Icp=",Ip*1E6,"microAmp")

##########################3
#Modeled PN of Loop Filter
Vn=4*(1.38E-23)*300*r1
num8=np.array([Vn])
den8=np.array([1])
PN_LF=tf(num8,den8)
print(PN_LF)
mag_LF,phase,w=bode(PN_LF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_LF),label="loop filter")
plt.grid(True,which="both",axis="both")
plt.title("Loop filter Phase Noise Curve")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")
#####################


##DEFINING TF OF THE BLOCKS
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


#Defining Reference TF H_REF(s)
H_REF=H/(R_divider*P_divider)
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


#Defining REF_Divider & N-Divider TF H_R_Divider(s)
H_R_Divider=H/(P_divider)
print("H_REF(S)=",H_R_Divider)
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



##PLOT OF TRANSFER FUNCTION OF EACH BLOCK
##########################################

#Plot of the TF seen by the Reference
mag_TF_REF,phase,w=bode(H_REF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_REF**2))
plt.grid(True,which="both",axis="both")
plt.title("Reference Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################


#Plot of the TF seen by the VCO
mag_TF_VCO,phase,w=bode(H_VCO,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_VCO**2))
plt.grid(True,which="both",axis="both")
plt.title("VCO Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################


#Plot of the TF seen by the Charge Pump (CP)
mag_TF_CP,phase,w=bode(H_CP,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_CP**2))
plt.grid(True,which="both",axis="both")
plt.title("Charge Pump (CP) Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################

#Plot of the TF seen by the Loop Filter (LF)
mag_TF_LF,phase,w=bode(H_LF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_LF**2))
plt.grid(True,which="both",axis="both")
plt.title("Loop Filter (LF) Transfer Function")
plt.xlabel("Frequency(Hz)")
plt.ylabel("(dB)")
#####################


#####################
#Plot of the TF seen by the Reference Divider & N Divider
mag_TF_R_Divider,phase,w=bode(H_R_Divider,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*12E3,2*np.pi*20E6])


plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),10*np.log10(mag_TF_R_Divider**2))
plt.grid(True,which="both",axis="both")
plt.title("(Reference & N) Divider Transfer Function")
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




########## PN OF EACH BLOCK @OUTPUT

#Plot of the PN of the Reference @Output
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

#Plot of the PN of the R-Divider @Output
PN_R_Divider_OUT_mag=(mag_TF_R_Divider**2)*(mag_Dividers)

PN_R_Divider_OUT=10*np.log10(mag_TF_R_Divider**2)+10*np.log10(mag_Dividers)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),PN_R_Divider_OUT,label="R-DIV")
plt.grid(True,which="both",axis="both")
plt.title(" R-Divider Output Phase Noise")
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


#####################


#TOTAL OUTPUT PHASE NOISE
Total_Out_PN_mag=PN_Buffer_OUT_mag+PN_P_Divider_OUT_mag+PN_N_Divider_OUT_mag+PN_R_Divider_OUT_mag+PN_CP_OUT_mag+PN_LF_OUT_mag+PN_VCO_OUT_mag+PN_REF_OUT_mag
Total_Out_PN=10*np.log10(Total_Out_PN_mag)

#plt.figure(figsize=(10, 8))
plt.semilogx(w/(2*np.pi),Total_Out_PN,label="Total")
plt.grid(True,which="both",axis="both")
plt.title("Phase Noise @ The Output")
plt.xlabel("Frequency(Hz)")
plt.ylabel("Phase Noise(dBc/Hz)")

plt.legend()

####################################################

####JITTER CALCULATION
#TOTAL RMS JITTER
f=np.logspace(np.log10(12E3), np.log10(20E6),1000)
Total_Integral=np.trapz(10**(Total_Out_PN/10),x=f)#x=(w/(2*np.pi)))
Total_Integral=2*Total_Integral
Total_RMS_Jitter=np.sqrt(Total_Integral)/(2*np.pi*Fout)
print("Total RMS Jitter= ",Total_RMS_Jitter*1E15,"fsec")


#REF RMS JITTER
REF_Integral=np.trapz(10**(PN_REF_OUT/10),x=f)
REF_Integral=2*REF_Integral
REF_RMS_Jitter=np.sqrt(REF_Integral)/(2*np.pi*Fout)
print("Reference RMS Jitter= ",REF_RMS_Jitter*1E15,"fsec")


#VCO RMS JITTER
VCO_Integral=np.trapz(10**(PN_VCO_OUT/10),x=f)
VCO_Integral=2*VCO_Integral
VCO_RMS_Jitter=np.sqrt(VCO_Integral)/(2*np.pi*Fout)
print("VCO RMS Jitter= ",VCO_RMS_Jitter*1E15,"fsec")



#CP RMS JITTER
CP_Integral=np.trapz(10**(PN_CP_OUT/10),x=f)
CP_Integral=2*CP_Integral
CP_RMS_Jitter=np.sqrt(CP_Integral)/(2*np.pi*Fout)
print("Charge Pump(CP) RMS Jitter= ",CP_RMS_Jitter*1E15,"fsec")



#LF RMS JITTER
LF_Integral=np.trapz(10**(PN_LF_OUT/10),x=f)
LF_Integral=2*LF_Integral
LF_RMS_Jitter=np.sqrt(LF_Integral)/(2*np.pi*Fout)
print("Loop Filter(LF) RMS Jitter= ",LF_RMS_Jitter*1E15,"fsec")


#N-Divider RMS JITTER
N_Divider_Integral=np.trapz(10**(PN_N_Divider_OUT/10),x=f)
N_Divider_Integral=2*N_Divider_Integral
N_Divider_RMS_Jitter=np.sqrt(N_Divider_Integral)/(2*np.pi*Fout)
print("N-Divider RMS Jitter= ",N_Divider_RMS_Jitter*1E15,"fsec")


#R-Divider RMS JITTER
R_Divider_Integral=np.trapz(10**(PN_R_Divider_OUT/10),x=f)
R_Divider_Integral=2*R_Divider_Integral
R_Divider_RMS_Jitter=np.sqrt(R_Divider_Integral)/(2*np.pi*Fout)
print("R-Divider RMS Jitter= ",R_Divider_RMS_Jitter*1E15,"fsec")


#P-Divider RMS JITTER
P_Divider_Integral=np.trapz(10**(PN_P_Divider_OUT/10),x=f)
P_Divider_Integral=2*P_Divider_Integral
P_Divider_RMS_Jitter=np.sqrt(P_Divider_Integral)/(2*np.pi*Fout)
print("P-Divider RMS Jitter= ",P_Divider_RMS_Jitter*1E15,"fsec")


#Buffer RMS JITTER
Buffer_Integral=np.trapz(10**(PN_Buffer_OUT/10),x=f)
Buffer_Integral=2*Buffer_Integral
Buffer_RMS_Jitter=np.sqrt(Buffer_Integral)/(2*np.pi*Fout)
print("Buffer RMS Jitter= ",Buffer_RMS_Jitter*1E15,"fsec")
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

#R-Divider Contribution
R_Divider_Contribution_Precentage=(R_Divider_Integral/Total_Integral)*100
print("R-Divider Contribution Precentage= ",R_Divider_Contribution_Precentage,"%")


#P-Divider Contribution
P_Divider_Contribution_Precentage=(P_Divider_Integral/Total_Integral)*100
print("P-Divider Contribution Precentage= ",P_Divider_Contribution_Precentage,"%")

#Buffer Contribution
Buffer_Contribution_Precentage=(Buffer_Integral/Total_Integral)*100
print("Buffer Contribution Precentage= ",Buffer_Contribution_Precentage,"%")


plt.show()




