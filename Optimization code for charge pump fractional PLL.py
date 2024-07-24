import numpy as np
from control.matlab import *
import sympy as sp

def Omega_Optimum(fu_inital,Fref,PM_deg,Kvco,Fout,P_divider,N,Gain):

    Total_RMS_Jitter_final=list()
    j_final=list()
    l_final=list()
    Ip=np.linspace(600E-6, 700E-6,30)
    f_start=fu_inital/1E6
    fu=np.linspace(f_start-.55,f_start+0.5,40)
    omega_u=1E6*fu*2*np.pi
    for y,j in enumerate(omega_u):
        for z,l in enumerate(Ip):
            #Automating Computation of X
            PM_rad=np.pi/180 *float(PM_deg) #take as input from user 
            val=sp.tan(PM_rad) #represent the value that will be used in eqn to compute X
            y=sp.Symbol('y')
            soln=sp.solve(sp.Eq(0.5*y-0.5*1/y,val),y)              
            for i in range (1,len(soln)):
                if(soln[i]>0):
                    z=soln[i]
                    
            
            #Calculating Caps and Resistors
            x=float(z)
            omega_z=j/x
            omega_p=j*x
            K=l*Kvco
            Beta=1/N
        
            r1 = j*N*(1+1/((x**2)-1))/K
            c1=1/(omega_z*r1)
            c2 = c1 /((x**2) - 1)
            c_sum=c1+c2
            if c_sum>1E-9:
                continue
            
            
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
            ##################
            
            
            
            ###
            #Defining loop TF A(s)
            num3=np.array([K])
            den3=np.array([1,0])
            H3=tf(num3,den3)
            
            A=series(Z,H3)
            #################
            
            ###
            #Defining loop TF T(s)
            T=A*Beta
            #################    
            
            ####
            #Defining Close_loop TF H(s)
            H=A/(1+T)
            ################
            
            
            #Defining Reference TF H_REF(s)
            H_REF=H/(P_divider)
            #print("H_REF(S)=",H_REF)
            ###########
            
            #Defining VCO TF H_VCO(s)
            H_VCO=1/(P_divider*(1+T))
            ###########
            
            #Defining CP TF H_CP(s)
            H_CP=2*H*np.pi/(l*P_divider)
            ###########
            
            
            #Defining LF TF H_LF(s)
            
            s=TransferFunction.s
            H_LF=2*np.pi*Kvco/(P_divider*s*(1+T))
            ###########
            
            #Defining N-Divider TF H_N_Divider(s)
            H_N_Divider=H/(P_divider)
            ###########
            
            #Defining DSM TF H_DSM(s)
            H_DSM=H/(P_divider)
            ###########
            
            #Defining O/P Buffer & P-divider TF H_Buffer(s)
            H_Buffer=1
            
            
            # #Defining LPF TF H_LPF(s)
            wp_lpf=2*np.pi*10E3
            num10=np.array([1])
            den10=np.array([1/wp_lpf,1])
            H_LPF=tf(num10,den10)
            ############################################
        
            
            
            ##MODELED PHASE NOISE
            
            ########################
            #Modeled PN of Crystal Refrence
            beta1=6.22036E-10
            beta2=0.0248/(2*np.pi*1E3)
            beta3=0.0248
            NoiseFloor_REF=1E-16
            
            num4=np.array([NoiseFloor_REF,beta1,beta2,beta3])
            den4=np.array([1,0,0,0])
            PN_REF=tf(num4,den4)
            mag_Ref,phase,w=bode(PN_REF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            
            ###########
            
            #Modeled PN of VCO
            num5=np.array([124.685,39305502.88])
            den5=np.array([1,0,0,0])
            PN_VCO=tf(num5,den5)
            mag_VCO,phase,w=bode(PN_VCO,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            for i in range(len(w)):
                if 10*np.log10(mag_VCO[i])<= -145:
                    mag_VCO[i:]=mag_VCO[i]
        
            #####################
            
            #Modeled PN of CP
            alpha=(7.943E-24)*((300E-6)**2)/(l**2)
            num6=np.array([alpha,alpha*((2*np.pi*2E6))])
            den6=np.array([1,0])
            PN_CP=tf(num6,den6)
            mag_CP,phase,w=bode(PN_CP,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            # #####################
            
            
            
            #Modeled PN of Loop Filter
            Vn=4*1.38E-23*300*r1
            num8=np.array([Vn])
            den8=np.array([1])
            PN_LF=tf(num8,den8)
            mag_LF,phase,w=bode(PN_LF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            #####################
            
            
            #Modeled PN of Dividers and O/P Buffer
            alpha1=1.2566E-9
            num7=np.array([alpha1/(2*np.pi*2E6),alpha1])
            den7=np.array([1,0])
            PN_Dividers=tf(num7,den7)
            mag_Dividers,phase,w=bode(PN_Dividers,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            #####################
            
            #DSM PHASE_NOISE
            Sq=1/(12*Fref)
            f=w/(2*np.pi)#np.logspace(2,7.60206,1000)
            NTF=np.abs((2*np.sin(np.pi*f/Fref))**4)
            Sg=Sq*NTF
            S_phi_n=Sg*((Fref/f)**2)/(N**2)
            
            
            #Modeled PN of I/P Referred Noise
            alpha2=10**(-135/10)
            num9=np.array([alpha2,alpha2*(2*np.pi*1E4)])
            den9=np.array([1,0])
            PN_Input_Referred=tf(num9,den9)
            mag_Input_Referred_Noise,phase,w=bode(PN_Input_Referred,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
    
            
            ##PLOT OF TRANSFER FUNCTION OF EACH BLOCK
            ##########################################
            
            #Plot of the TF seen by the Reference
            mag_TF_REF,phase,w=bode(H_REF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            #####################
            
            
            #Plot of the TF seen by the VCO
            mag_TF_VCO,phase,w=bode(H_VCO,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            #####################
            
            
            #Plot of the TF seen by the Charge Pump (CP)
            mag_TF_CP,phase,w=bode(H_CP,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            #####################
            
            #Plot of the TF seen by the Loop Filter (LF)
            mag_TF_LF,phase,w=bode(H_LF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            #####################
            
            
            #####################
            #Plot of the TF seen by the N Divider
            mag_TF_N_Divider,phase,w=bode(H_N_Divider,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            #####################
            
            #Plot of the TF seen by DSM
            mag_TF_DSM,phase,w=bode(H_DSM,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            
            
            #Plot of the TF seen by the O/P Buffer & P-Divider
            mag_TF_Buffer_dB=10*np.log10(H_Buffer**2)*(w/w)
            
            #####################
            #Plot of the LPF TF
            mag_TF_LPF,phase,w=bode(H_LPF,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])
            num=[2*np.pi*Fref/N]
            den=[1,0]
            integr=tf(num,den)
            mag_TF_integr,phase,w=bode(integr,dB=True,deg=True,plot=False,Hz=False,omega_limits=[2*np.pi*1E3,2*np.pi*40E6])

            mag_VCXO=((Gain*mag_TF_LPF)**2)*(mag_Input_Referred_Noise)*(mag_TF_integr**2)    
            
            ########## PN OF EACH BLOCK @OUTPUT
            
            # PN of the Reference @Output
            PN_REF_OUT_mag=(mag_TF_REF**2)*(mag_Ref)
            
            PN_REF_OUT=10*np.log10(mag_TF_REF**2)+10*np.log10(mag_Ref)
            
        
            #####################
            
            
            
            #PN of the VCO @Output
            PN_VCO_OUT_mag=(mag_TF_VCO**2)*(mag_VCO)
            
            PN_VCO_OUT=10*np.log10(mag_TF_VCO**2)+10*np.log10(mag_VCO)
            
            #####################
            
            
            #PN of the CP @Output
            PN_CP_OUT_mag=(mag_TF_CP**2)*(mag_CP)
            
            PN_CP_OUT=10*np.log10(mag_TF_CP**2)+10*np.log10(mag_CP)
            
            #####################
            
            
            #PN of the LF @Output
            PN_LF_OUT_mag=(mag_TF_LF**2)*(mag_LF)
            
            PN_LF_OUT=10*np.log10(mag_TF_LF**2)+10*np.log10(mag_LF)    
            #####################
            
            
            #Plot of the PN of the N-Divider @Output
            PN_N_Divider_OUT_mag=(mag_TF_N_Divider**2)*(mag_Dividers)
            
            PN_N_Divider_OUT=10*np.log10(mag_TF_N_Divider**2)+10*np.log10(mag_Dividers)
            
            #####################
            
            #Plot of the PN of the DSM @Output
            PN_DSM_OUT_mag=(mag_TF_DSM**2)*(S_phi_n)
            
            PN_DSM_OUT=10*np.log10(mag_TF_DSM**2)+10*np.log10(S_phi_n)
                    
            
            #####################
            
            
            #Plot of the PN of the P-Divider @Output
            PN_P_Divider_OUT_mag=(10**(mag_TF_Buffer_dB/10))*(mag_Dividers)
            
            PN_P_Divider_OUT=mag_TF_Buffer_dB+10*np.log10(mag_Dividers)
            
            #####################
            
            
            #Plot of the PN of the Output Buffer @Output
            PN_Buffer_OUT_mag=(10**(mag_TF_Buffer_dB/10))*(mag_Dividers)
            
            PN_Buffer_OUT=mag_TF_Buffer_dB+10*np.log10(mag_Dividers)
            
            #####################
            
            #Plot of the VCXO I/P Referred Noise @Output
            PN_VCXO_Noise_OUT_mag=mag_VCXO*(mag_TF_DSM**2)
            
            PN_VCXO_Noise_OUT=10*np.log10(mag_VCXO)+10*np.log10(mag_TF_DSM**2)
            
            ################################3
            #TOTAL OUTPUT PHASE NOISE
            Total_Out_PN_mag=PN_VCXO_Noise_OUT_mag+PN_DSM_OUT_mag+PN_Buffer_OUT_mag+PN_P_Divider_OUT_mag+PN_N_Divider_OUT_mag+PN_CP_OUT_mag+PN_LF_OUT_mag+PN_VCO_OUT_mag+PN_REF_OUT_mag
            Total_Out_PN=10*np.log10(Total_Out_PN_mag)
            
            ####################################################
            
            ####JITTER CALCULATION
            #TOTAL RMS JITTER
            f1=np.logspace(np.log10(12E3), np.log10(20E6),1000)
            Total_Integral=np.trapz(10**(Total_Out_PN/10),x=f1)
            Total_Integral=2*Total_Integral
            Total_RMS_Jitter=np.sqrt(Total_Integral)/(2*np.pi*Fout)
            
            Total_RMS_Jitter_final.append(Total_RMS_Jitter)
            j_final.append(j/(2*np.pi*1E6))
            l_final.append(l)



    min_jitter=np.min(Total_RMS_Jitter_final)
    index_of_min_jitter=Total_RMS_Jitter_final.index(min_jitter)
    return [j_final[index_of_min_jitter]*1E6,l_final[index_of_min_jitter]]




