import numpy as np
import control
import matplotlib.pyplot as plt
pi=np.pi
def Bode_10dB(Noise_floor,H,block_Name,Hz,TF):
    mag, phase, w = control.bode_plot(H,plot=False,omega_limits=(100,10**9))
    mag=(mag+Noise_floor)*TF
    magnitude = 10 * np.log10(mag)
    ww=w
    if Hz==1:
        ww=w/(2*pi)
    # Plot the Bode plot
    ##plt.figure(figsize=(10, 8))

    plt.semilogx(ww,magnitude,label=block_Name)
    plt.xlabel('Frequency (rad/sec)')
    if Hz==1:
        plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dBc/Hz)')
    plt.legend(fontsize="13")
    plt.grid(True,which="both",axis="both")
    return mag, phase, w
def Bode2_10dB(mag,w,block_Name,Hz):
    # Calculate magnitude in dB using 10log(x)
    magnitude = 10 * np.log10(mag)
    ww=w
    if Hz==1:
        ww=w/(2*pi)
    # Plot the Bode plot
    ##plt.figure(figsize=(10, 8))

    plt.semilogx(ww,magnitude,label=block_Name)
    plt.xlabel('Frequency (rad/sec)')
    plt.legend(fontsize="13")
    plt.grid(True,which="both",axis="both")
    if Hz==1:
        plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dBc/Hz)')


s=control.tf('s')
Fref=(150*10**6)
Fref2=(100*10**6)
Tref=1/Fref
Tref2=1/Fref2
TDC_res=0.1*10**-12
TDC_res2=10*10**-12
Fc=1.5*10**6 #assume for now
PM=pi/2.6 #phase margin = 69 degrees
P=4
N=64.5
N1=18
N2=5
Kvco=55*10**6 #in Hz/V
Kdco=(1.2/2**N1)*Kvco*2*pi #in rad/s/LSB
print(Kdco)
Fout=Fref*N/P
Fvco=Fref*N
DTC_res=1/((2**10)*Fvco)

## CALCULATING ALPHA AND BETA
B_A=(2*pi*Fc*Tref)/np.tan(PM)
Beta=(2*pi*TDC_res*((2*pi*Fc)**2)*N)/(np.sqrt((np.tan(PM)**2)+1)*Kdco)
alpha=Beta*(1/B_A)
print('Beta = ',Beta)
print('Alpha = ',alpha)

## LOOP GAIN EQUATION
Wz=(Beta*Fref)/alpha
b=(Kdco*Beta)/(2*pi*TDC_res*N)
LG=b*(1+s/Wz)*(1/(s*s))
#plt.figure(figsize=(10, 8))

control.bode_plot(LG, dB=True,Hz=1,omega_limits=(100,10**9),plot=1)
plt.title('Open loop transfer function')
x=control.margin(LG)
print('Phase Margin = ',x[1])
print('Bandwidth = ',x[3]/(2*pi*10**6),' Mhz')


## CLOSED LOOP GAIN EGUATION
CL=LG*N/(1+LG)
#plt.figure(figsize=(10, 8))

control.bode_plot(CL, dB=True,Hz=1,omega_limits=(100,10**9),plot=1)
plt.title('Closed loop transfer function')


###########################.noise transfer function of each block.####################################

#Refrence
Ref=CL/P
##plt.figure(figsize=(10, 8))

mag_Ref, phase, w =control.bode_plot(Ref,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of Refrence')

#N-Divider
N_div=CL/P
##plt.figure(figsize=(10, 8))

mag_N_div, phase, w = control.bode_plot(N_div,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of N-Divider')

#VCO
VCO=(1/(1+LG))*(1/P)
##plt.figure(figsize=(10, 8))

mag_VCO, phase, w = control.bode_plot(VCO,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of VCO')

#TDC
TDC=(CL*(2*pi*TDC_res/Tref))/P
##plt.figure(figsize=(10, 8))

mag_TDC, phase, w = control.bode_plot(TDC,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of TDC')

#DAC
DAC=(1/(1+LG))*(Kvco/s)*(1/P)
##plt.figure(figsize=(10, 8))

mag_DAC, phase, w = control.bode_plot(DAC,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of DCO QN')

#DTC
DTC=CL/P
##plt.figure(figsize=(10, 8))

mag_DTC, phase, w = control.bode_plot(DTC,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of DTC')

######################################.loop equation of 2nd loop.#####################################

Fref2=(100*10**6)
Tref2=1/Fref2
Fc2=100 #assume for now
N_div2=96.75
Nsm=24
Kvco2=150*10**6 #in Hz/V
Kdco2=(1/2**Nsm)*Kvco2*2*pi #in rad/s/LSB
## CALCULATING ALPHA AND BETA
B_A2=(2*pi*Fc2*Tref2)/np.tan(PM)
Beta2=(2*pi*TDC_res2*((2*pi*Fc2)**2)*N_div2)/(np.sqrt((np.tan(PM)**2)+1)*Kdco2)
alpha2=Beta2*(1/B_A2)
print('Beta = ',Beta2)
print('Alpha = ',alpha2)

LG2=(Tref2/(2*pi*TDC_res2))*(alpha2+(Beta2/(s*Tref2)))*(Kdco2/s)*1/N_div2
#plt.figure(figsize=(10, 8))
x=control.margin(LG2)
print('Phase Margin = ',x[1])
print('Bandwidth = ',x[3]/(2*pi*10**3),' khz')
control.bode_plot(LG2,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Loop-gain 2')

CL2=(N_div2*LG2)/(1+LG2)

#plt.figure(figsize=(10, 8))

control.bode_plot(CL2,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Closed-loop 2')

###############################.Transfer function of each block.#######################################

plt.figure(figsize=(10, 8))

#IP
IP=CL2/P
mag_IP, phase, w =control.bode_plot(IP,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of IP')

#TDC2
TDC2=(CL2*(2*pi*TDC_res2/Tref2))/P
plt.figure(figsize=(10, 8))

mag_TDC2, phase, w = control.bode_plot(TDC2,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of TDC2')

#DCO
DCO=(1/(1+LG2))*(1/P)
plt.figure(figsize=(10, 8))

mag_DCO, phase, w = control.bode_plot(DCO,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of DCO')

#plt.show()
#################################.modeled phase noise of each block.#####################################
##plt.figure(figsize=(10, 8))


#DSM_Noise
Sg=1/(12*Fref)
num=(4*(pi**2))*(s**2)
den=(N**2)*(Fref**2)#*(2**10)
DSM=Sg*num/den
Bode_10dB(0,DSM,'DSM noise spectrum',1,1)

#DAC_QN
delta_f=(1.2*Kvco)/(2**N2)
DAC_QN=((4*pi**2)*(s**2)*delta_f**2)/((Fref**5)*12)
Bode_10dB(0,DAC_QN,'DCO QN noise spectrum',1,1)

#IP_Noise
plt.figure(figsize=(10, 8))

A=1.17588
w11=2*pi*124.53
w22=2*pi*257.08
w21=2*pi*1*10**5
num3=(1+s/w11)*(1+s/w11)*(1+s/w11)
den3=(s*s*s)*(1+s/w22)*(1+s/w21)*(1+s/w21)
IP_N=(num3/den3)*A
Noise_floor_IP=10**(-11.4)
Bode_10dB(Noise_floor_IP,IP_N,'Jittery_I/P noise spectrum',1,1)

#TDC1_QN
plt.figure(figsize=(10, 8))

TDC_QN=(((2*pi*TDC_res)/Tref)**2)*(s**0)/(12*Fref)
Bode_10dB(0,TDC_QN,'TDC QN noise spectrum',1,1)

#TDC2_QN
TDC2_QN=(((2*pi*TDC_res2)/Tref2)**2)*(s**0)/(12*Fref2)
Bode_10dB(0,TDC2_QN,'TDC QN noise spectrum',1,1)

#DTC_QN
DTC_QN=((2*pi*DTC_res*Fvco)**2)*(Tref/12)*(s**2)*(Tref**2)/(N**2)
Bode_10dB(0,DTC_QN,'DTC QN noise spectrum',1,1)

#Buffer & Dividers
num1 =(1+(s/(2*pi*2*(10**6))))
den1 =s
H =num1/den1
H=H*1.22*10**-9
#print ('H(s) =', H)
Bode_10dB(0,H,'Buffer & Dividers noise spectrum',1,1)

#VCO
w1=2*pi*10**6
w2=pi*2*3*10**7
num2=(1+s/w1)*(1+s/w2)*(1+s/w2)
den2=(s*s*s)
g=(num2/den2)*10**9.05
Bode_10dB(0,g,'VCO noise spectrum',1,1)

#Crystal_Refrence
w11=2*pi*10**3
w22=2*pi*30*10**3
num3=(1+s/w11)*(1+s/w22)*(1+s/w22)
den3=(s*s*s)
G=(num3/den3)*10**-1.85
Bode_10dB(0,G,'Crystal_Refrence noise spectrum',1,1)

#DTC and TDC noise
Conv_noise = 2.50732E-9 * ((s/(2*np.pi*2E6)) + 1)/s
Bode_10dB(0,Conv_noise,'Data Converters Noise spectrum',1,1)

plt.show()
#####################################.close noise VS VCO noise.#########################################
##plt.figure(figsize=(10, 8))


#P-divider phase noise curve
mag1_Pdiv,phase1,w=Bode_10dB(0,H,'P-divider phase noise curve',Hz=1,TF=N**2)

#N-divider phase noise curve
mag1_Ndiv,phase1,w=Bode_10dB(0,H,'N-divider phase noise curve',Hz=1,TF=N**2)

#buffer phase noise curve
mag1_buff,phase1,w=Bode_10dB(0,H,'Buffer phase noise curve',Hz=1,TF=N**2)

#Refrence phase noise curve
mag1_Ref,phase1,w=Bode_10dB(0,G,'Crystal-Refrence phase noise curve',Hz=1,TF=N**2)

#DCO QN curve
mag1_DAC_QN,phase1,w=Bode_10dB(0,DAC_QN,'DCO QN curve',Hz=1,TF=N**2)

#TDC QN curve
mag1_TDC_QN,phase1,w=Bode_10dB(0,TDC_QN,'TDC QN noise curve',Hz=1,TF=N**2)

#DTC QN curve
mag1_DTC_QN,phase1,w=Bode_10dB(0,DTC_QN,'DTC QN noise curve',Hz=1,TF=N**2)

#VCO phase noise curve
mag1_VCO,phase,W=Bode_10dB(0,g,'VCO phase noise curve',Hz=1,TF=1)

#TDC and DTC phase noise curve
mag1_Conv_N,phase1,w=Bode_10dB(0,Conv_noise,'TDC and DTC phase noise curve',Hz=1,TF=N**2)

#DSM phase noise curve
mag1_DSM,phase1,w=Bode_10dB(0,DSM,'DSM phase noise curve',Hz=1,TF=N**2)


Mag1_close_noise=2*mag1_Conv_N+mag1_Ref+mag1_Ndiv+mag1_buff+mag1_Pdiv+mag1_TDC_QN+mag1_DAC_QN+mag1_DTC_QN#+mag1_DSM

#plt.figure(figsize=(10, 8))

plt.title('Wu Optimized at intersection ')
Bode2_10dB(Mag1_close_noise,W,'Close noise',Hz=1)
Bode2_10dB(mag1_VCO,W,'VCO',Hz=1)

##################################.Total Noise at output.############################################
##plt.figure(figsize=(10, 8))

#plt.title('Phase Noise of each block at output')

#P-divider phase noise curve
#plt.figure(figsize=(10, 8))

mag2_Pdiv,phase1,w=Bode_10dB(0,H,'P-divider phase noise curve',Hz=1,TF=1)

#N-divider phase noise curve
#plt.figure(figsize=(10, 8))

mag2_Ndiv,phase1,w=Bode_10dB(0,H,'N-divider phase noise curve',Hz=1,TF=mag_N_div**2)

#buffer phase noise curve
#plt.figure(figsize=(10, 8))

mag2_buff,phase1,w=Bode_10dB(0,H,'Buffer phase noise curve',Hz=1,TF=1)

#Refrence phase noise curve
#plt.figure(figsize=(10, 8))

mag2_Ref,phase1,w=Bode_10dB(0,G,'Crystal-Refrence phase noise curve',Hz=1,TF=mag_Ref**2)

#DCO QN curve
#plt.figure(figsize=(10, 8))

mag2_DAC_QN,phase1,w=Bode_10dB(0,DAC_QN,'DCO QN curve',Hz=1,TF=mag_DAC**2)

#TDC QN curve
#plt.figure(figsize=(10, 8))

mag2_TDC_QN,phase1,w=Bode_10dB(0,TDC_QN,'TDC QN noise curve',Hz=1,TF=mag_DTC**2)

#DTC QN curve
#plt.figure(figsize=(10, 8))

mag2_DTC_QN,phase1,w=Bode_10dB(0,DTC_QN,'DTC QN noise curve',Hz=1,TF=mag_DTC**2)

#VCO phase noise curve
#plt.figure(figsize=(10, 8))

mag2_VCO,phase,W=Bode_10dB(0,g,'VCO phase noise curve',Hz=1,TF=mag_VCO**2)

#DTC phase noise curve
#plt.figure(figsize=(10, 8))

mag2_Conv_N,phase1,w=Bode_10dB(0,Conv_noise,'DTC phase noise curve',Hz=1,TF=mag_DTC**2)

#TDC phase noise curve
#plt.figure(figsize=(10, 8))

mag22_Conv_N,phase1,w=Bode_10dB(0,Conv_noise,'TDC phase noise curve',Hz=1,TF=mag_TDC**2)

#DSM1 phase noise curve
#mag2_DSM,phase1,w=Bode_10dB(0,DSM,'DSM phase noise curve',Hz=1,TF=mag_N_div**2)
#plt.figure(figsize=(10, 8))

#N-divider2 phase noise curve
mag2_Ndiv2,phase1,w=Bode_10dB(0,H,'N-divider phase noise curve',Hz=1,TF=mag_IP**2)

#IP phase noise curve
mag2_IP,phase1,w=Bode_10dB(Noise_floor_IP,IP_N,'IP phase noise curve',Hz=1,TF=mag_IP**2)

#TDC2 phase noise curve
mag22_Conv_N2,phase1,w=Bode_10dB(0,Conv_noise,'TDC phase noise curve',Hz=1,TF=mag_TDC2**2)

#TDC2 QN noise curve
mag2_TDC_QN2,phase1,w=Bode_10dB(0,TDC2_QN,'TDC phase noise curve',Hz=1,TF=mag_TDC2**2)

#DSM2 phase noise curve
#plt.figure(figsize=(10, 8))

mag2_DSM,phase1,w=Bode_10dB(0,DSM,'DSM phase noise curve',Hz=1,TF=mag_IP**2)


Mag2_close_noise=mag22_Conv_N+mag2_Conv_N+mag2_Ref+mag2_Ndiv+mag2_buff+mag2_TDC_QN+mag2_Pdiv+mag2_DAC_QN+mag2_DTC_QN#+mag2_DSM
Mag_total=Mag2_close_noise+mag2_VCO
#Mag_total2=Mag_total*mag_DCO**2

Mag_loop2=mag2_Ndiv2+mag2_IP+mag22_Conv_N2+mag2_TDC_QN2+mag2_DSM

#plt.figure(figsize=(10, 8))


#Bode2_10dB(Mag_loop2,W,'output reffered noise of loop 2',Hz=1)

plt.figure()
plt.title('noise at output')
#Bode2_10dB(Mag_total,W,'output reffered noise of internal loop',Hz=1)

Mag_total2=Mag_total+Mag_loop2
Bode2_10dB(Mag_total2,W,'total output referred noise',Hz=1)
Bode2_10dB(mag2_IP,W,'jittery I/P noise output referred ',Hz=1)
Bode_10dB(Noise_floor_IP,IP_N,'Jittery I/P noise spectrum',1,1)
plt.show()

filterr=np.zeros(1000)
filterr[411:872]=1
W_Int=W*filterr/(2*pi)
Mag_Int=Mag_total*filterr
RMS_sq=np.trapz(Mag_Int,W_Int)
RMS_Jitter=np.sqrt(2*(RMS_sq))/(2*pi*2.42*10**9)

Mag_Int2=Mag_loop2*filterr
RMS_sq2=np.trapz(Mag_Int2,W_Int)
RMS_Jitter2=np.sqrt(2*(RMS_sq2))/(2*pi*2.42*10**9)
print('\n')
print('Total RMS jitter from VCXO is',RMS_Jitter)
print('\n')
print('Total RMS jitter from Loop 2 is ',RMS_Jitter2)
RMS_Final=RMS_sq+RMS_sq2
RMS_finall=np.sqrt(2*(RMS_Final))/(2*pi*2.42*10**9)#edit##
print('\n')
print('Total RMS jitter is ',RMS_finall)
RMS_sq_tot=RMS_sq+RMS_sq2
RMS_total=RMS_Jitter+RMS_Jitter2
#############################################.percentage of each block.########################################
#p-divider
P_divider=mag2_Pdiv*filterr
pdiv=np.trapz(P_divider,W_Int)
RMS_pdiv=np.trapz(Mag_Int2,W_Int)
RMS_Jitter2=np.sqrt(2*(RMS_pdiv))/(2*pi*2.42*10**9)
print('\n')
print('RMS jitter from P-divider is',RMS_Jitter)


#N-divider
N_divider=mag2_Ndiv*filterr
ndiv=np.trapz(N_divider,W_Int)

#Crystal-Refrence
refrence=mag2_Ref*filterr
ref=np.trapz(refrence,W_Int)

#VCO
VCO_out=mag2_VCO*filterr
vco=np.trapz(VCO_out,W_Int)

#DAC
DACC=mag2_DAC_QN*filterr
dac=np.trapz(DACC,W_Int)

#TDC
TDCC=(mag2_TDC_QN+mag22_Conv_N)*filterr
tdc=np.trapz(TDCC,W_Int)

#DTC
DTCC=(mag2_DTC_QN+ mag2_Conv_N)*filterr
dtc=np.trapz(DTCC,W_Int)

#Buffer
Buffer=mag2_buff*filterr
buff=np.trapz(Buffer,W_Int)

#N-Divider2
N_divider2=mag2_Ndiv2*filterr
N__div2=np.trapz(N_divider2,W_Int)

#DSM2
DSM2=mag2_DSM*filterr
DSM22=np.trapz(DSM2,W_Int)

#TDC2
TDC22=(mag22_Conv_N2+mag2_TDC_QN2)*filterr
TDCC2=np.trapz(TDC22,W_Int)

#IP
IP_JIT=mag2_IP*filterr
IP_JIT2=np.trapz(IP_JIT,W_Int)


print('\n')
print('Total output integrated RMS jitter is ', RMS_total)
print('\n')
print('P-divider percentage in the total output integrated RMS jitter is ',(pdiv/RMS_sq_tot)*100,'%')
print('N-divider percentage in the total output integrated RMS jitter is ',(ndiv/RMS_sq_tot)*100,'%')
print('Buffer percentage in the total output integrated RMS jitter is ',(buff/RMS_sq_tot)*100,'%')
print('DAC percentage in the total output integrated RMS jitter is ',(dac/RMS_sq_tot)*100,'%')
print('TDC percentage in the total output integrated RMS jitter is ',(tdc/RMS_sq_tot)*100,'%')
print('DTC+DSM 1 percentage in the total output integrated RMS jitter is ',(dtc/RMS_sq_tot)*100,'%')
print('VCO percentage in the total output integrated RMS jitter is ',(vco/RMS_sq_tot)*100,'%')
print('Refrence percentage in the total output integrated RMS jitter is ',(ref/RMS_sq_tot)*100,'%')
print('N-divider 2 percentage in the total output integrated RMS jitter is ',(N__div2/RMS_sq_tot)*100,'%')
print('DSM 2 percentage in the total output integrated RMS jitter is ',(DSM22/RMS_sq_tot)*100,'%')
print('TDC 2 percentage in the total output integrated RMS jitter is ',(TDCC2/RMS_sq_tot)*100,'%')
print('Jittery I/P percentage in the total output integrated RMS jitter is ',(IP_JIT2/RMS_sq_tot)*100,'%')
