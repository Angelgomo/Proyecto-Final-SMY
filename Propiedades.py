import numpy as np


def krw_point(Sw, Srw, Sro, ommega):
    
    Se=(Sw-Srw)/(1-Srw-Sro)
    krw=Se**ommega
    
    return (krw)


def kro_point(Sw, Srw, Sro, ommega):
    
    Se=(Sw-Srw)/(1-Srw-Sro)
    kro=(1-Se)**ommega
    
    return (kro)


def krw_2D(Sw_e, Sw_w, Sw_n, Sw_s, Srw, Sro, ommega):
    
    Se_e=(Sw_e-Srw)/(1-Srw-Sro)
    krw_e=Se_e**ommega

    Se_w=(Sw_w-Srw)/(1-Srw-Sro)
    krw_w=Se_w**ommega

    Se_n=(Sw_n-Srw)/(1-Srw-Sro)
    krw_n=Se_n**ommega

    Se_s=(Sw_s-Srw)/(1-Srw-Sro)
    krw_s=Se_s**ommega

    return (krw_e, krw_w, krw_n, krw_s)



def kro_2D(Sw_e, Sw_w, Sw_n, Sw_s, Srw, Sro, ommega):
    
    Se_e=(Sw_e-Srw)/(1-Srw-Sro)
    kro_e=(1-Se_e)**ommega

    Se_w=(Sw_w-Srw)/(1-Srw-Sro)
    kro_w=(1-Se_w)**ommega


    Se_n=(Sw_n-Srw)/(1-Srw-Sro)
    kro_n=(1-Se_n)**ommega

    Se_s=(Sw_s-Srw)/(1-Srw-Sro)
    kro_s=(1-Se_s)**ommega

    return (kro_e, kro_w, kro_n, kro_s)


def upwind_2D(Nx,Ny,i, j, P, Sw):
    if((i>0) and (i<Nx-1) and (j>0) and (j<Ny-1)):
        
        if(P[i][j]<=P[i+1][j]):
            Sw_e=Sw[i+1][j]
        else:
            Sw_e=Sw[i][j]
    
        if(P[i][j]<=P[i-1][j]):
            Sw_w=Sw[i-1][j]
        else:
            Sw_w=Sw[i][j]
    
    
        if(P[i][j]<=P[i][j+1]):
            Sw_n=Sw[i][j+1]
        else:
            Sw_n=Sw[i][j]
    
        if(P[i][j]<=P[i][j-1]):
            Sw_s=Sw[i][j-1]
        else:
            Sw_s=Sw[i][j]
    
    
    if ((i==0)and(j>0)and(j<Ny-1)):
        if(P[i][j]<=P[i][j+1]):
            Sw_n=Sw[i][j+1]
        else:
            Sw_n=Sw[i][j]
    
        if(P[i][j]<=P[i][j-1]):
            Sw_s=Sw[i][j-1]
        else:
            Sw_s=Sw[i][j]
        if(P[i][j]<=P[i+1][j]):
            Sw_e=Sw[i+1][j]
        else:
            Sw_e=Sw[i][j]
            
        Sw_w = Sw [i][j]
    
    if((i==Nx-1)and(j>0)and(j<Ny-1)):
        if(P[i][j]<=P[i-1][j]):
            Sw_w=Sw[i-1][j]
        else:
            Sw_w=Sw[i][j]
    
    
        if(P[i][j]<=P[i][j+1]):
            Sw_n=Sw[i][j+1]
        else:
            Sw_n=Sw[i][j]
    
        if(P[i][j]<=P[i][j-1]):
            Sw_s=Sw[i][j-1]
        else:
            Sw_s=Sw[i][j]
        
        Sw_e = Sw [i][j]
    
    if((i>0)and(i<Nx-1)and(j==0)):
        if(P[i][j]<=P[i+1][j]):
            Sw_e=Sw[i+1][j]
        else:
            Sw_e=Sw[i][j]
    
        if(P[i][j]<=P[i-1][j]):
            Sw_w=Sw[i-1][j]
        else:
            Sw_w=Sw[i][j]
    
    
        if(P[i][j]<=P[i][j+1]):
            Sw_n=Sw[i][j+1]
        else:
            Sw_n=Sw[i][j]
        
        Sw_s = Sw[i][j]
        
    if((i>0)and(i<Nx-1)and(j==Ny-1)):
        if(P[i][j]<=P[i+1][j]):
            Sw_e=Sw[i+1][j]
        else:
            Sw_e=Sw[i][j]
    
        if(P[i][j]<=P[i-1][j]):
            Sw_w=Sw[i-1][j]
        else:
            Sw_w=Sw[i][j]
    
        if(P[i][j]<=P[i][j-1]):
            Sw_s=Sw[i][j-1]
        else:
            Sw_s=Sw[i][j]
    
        Sw_n = Sw[i][j]
    
#Esquinas
#(0,0)
    if((i==0)and(j==0)):
        if(P[i][j]<=P[i+1][j]):
            Sw_e=Sw[i+1][j]
        else:
            Sw_e=Sw[i][j]
        if(P[i][j]<=P[i][j+1]):
            Sw_n=Sw[i][j+1]
        else:
            Sw_n=Sw[i][j]
            
        Sw_w = Sw[i][j]
        Sw_s = Sw[i][j]
#(0,Ny-1)
    if((i==0)and(j==Ny-1)):
        if(P[i][j]<=P[i+1][j]):
            Sw_e=Sw[i+1][j]
        else:
            Sw_e=Sw[i][j]
        if(P[i][j]<=P[i][j-1]):
            Sw_s=Sw[i][j-1]
        else:
            Sw_s=Sw[i][j]
        
        Sw_n = Sw[i][j]
        Sw_w = Sw[i][j]
        
#(Nx-1,Ny-1)
    if((i==Nx-1)and(j==Ny-1)):
        if(P[i][j]<=P[i][j-1]):
            Sw_s=Sw[i][j-1]
        else:
            Sw_s=Sw[i][j]
        if(P[i][j]<=P[i-1][j]):
            Sw_w=Sw[i-1][j]
        else:
            Sw_w=Sw[i][j]
            
        Sw_n = Sw[i][j]
        Sw_e = Sw[i][j]

#(Nx-1,0)
    if((i==Nx-1)and(j==0)):
        if(P[i][j]<=P[i-1][j]):
            Sw_w=Sw[i-1][j]
        else:
            Sw_w=Sw[i][j]
        if(P[i][j]<=P[i][j+1]):
            Sw_n=Sw[i][j+1]
        else:
            Sw_n=Sw[i][j]
        
        Sw_e = Sw[i][j]
        Sw_s =Sw[i][j]
        
        
    return (Sw_e, Sw_w, Sw_n, Sw_s)




















krw=0.22
kro=20
Srw=0.20
Sro=0.15
Sw=np.linspace(Srw,1-Sro, 20)
ommega=1.7


# plt.figure("Permeabilidades relativas")
# plt.title("Krw-Kro- Vs Sw")
# plt.plot(Sw, krw_point(Sw, Srw, Sro, ommega), 'b--')
# plt.plot(Sw, kro_point(Sw, Srw, Sro, ommega), 'r--')
# plt.grid()
# plt.xlabel("$S_{w}$")
# plt.ylabel("$k_{rw}$")

# fig, ax1 = plt.subplots()
# ax1.plot(Sw, krw(Sw, Srw, Sro, ommega), 'b--')
# ax1.set_ylabel("$k_{rw}$")
# ax1.set_xlabel("$S_{w}$")
# ax2 = ax1.twinx();
# ax2.plot(Sw, kro(Sw, Srw, Sro, ommega), 'r--')
# ax2.set_ylabel("$k_{ro}$")
# plt.grid()

