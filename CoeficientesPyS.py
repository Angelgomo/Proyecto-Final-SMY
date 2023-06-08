import Propiedades as pp
import numpy as np
import Propiedades as fP


def coefPres_2D(Nx, Ny, x_centers, y_centers,  x_grid, y_grid, delta_z, rwell_inj, coordInj_i, coordInj_j, Pbh_inj, Sw_inj,
                rwell_prod, coordProd_i, coordProd_j, Pbh_prod, Sw_prod, Srw, Sro, ommega, Muw, Muo, kx, ky, P, Sw, AP, AE, AW, AN, AS, B):
    
    for i in range (0,Nx):    
        for j in range (0,Ny):
            
            if((i>0)and(i<Nx-1)and(j>0)and(j<Ny-1)):       #Distancias de volumenes y centros nodos internos          
                deltax_e=x_centers[i+1]-x_centers[i]      #Distancia del nodo P al nodo E
                deltax_w=x_centers[i]-x_centers[i-1]      #Distancia del nodo W al nodo P
                deltay_n=y_centers[j+1]-y_centers[j]      #Distancia del nodo P al nodo N
                deltay_s=y_centers[j]-y_centers[j-1]      #Distancia del nodo S al nodo P
                
                delta_x=x_grid[i]-x_grid[i-1]
                delta_y=y_grid[j]-y_grid[j-1]
                
                kx_e=kx; kx_w=kx   #Permeabilidad en x evaluada en las caras E y W
                ky_n=ky; ky_s=ky   #Permeabilidad en y evaluada en las caras N y S
            
          
## condiciones en las frontera cerradas (no flujo)
    
            if((i==0)and(j>0)and(j<Ny-1)):  #Distancias de volumenes y centros en frontera oeste (w)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=deltax_e   #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=y_centers[j]-y_centers[j-1]
                 
                 delta_x=x_grid[1]-x_grid[0]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=0.0  #Permeabilidad en x evaluada en las caras E y W (Frontera cerrada oeste)
                 ky_n=ky; ky_s=ky   #Permeabilidad en y evaluada en las caras N y S
            
            if((i==Nx-1)and(j>0)and(j<Ny-1)):  #Distancias de volumenes y centros en frontera este (e)
                 deltax_w=x_centers[i]-x_centers[i-1]
                 deltax_e=deltax_w   #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=y_centers[j]-y_centers[j-1]
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=0.0; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W (Frontera cerrada este)
                 ky_n=ky; ky_s=ky   #Permeabilidad en y evaluada en las caras N y S
                 
            if((i>0)and(i<Nx-1)and(j==0)):  #Distancias de volumenes y centros en frontera sur (s)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=x_centers[i]-x_centers[i-1] 
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=deltay_n #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[1]-y_grid[0]
                 
                 kx_e=kx; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W 
                 ky_n=ky; ky_s=0.0  #Permeabilidad en y evaluada en las caras N y S (Frontera cerrada sur)
                 
            if((i>0)and(i<Nx-1)and(j==Ny-1)):  #Distancias de volumenes y centros en frontera norte (n)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=x_centers[i]-x_centers[i-1] 
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W 
                 ky_n=0.0; ky_s=ky #Permeabilidad en y evaluada en las caras N y S (Frontera cerrada norte)
                 
            if((i>0)and(i<Nx-1)and(j==Ny-1)):  #Distancias de volumenes y centros en frontera norte (n)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=x_centers[i]-x_centers[i-1] 
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W 
                 ky_n=0.0; ky_s=ky #Permeabilidad en y evaluada en las caras N y S (Frontera cerrada norte)
           
            if((i==0)and(j==0)):  #Nodo de esquina (0,0)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=deltax_e
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=deltay_n #Se clona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[1]-x_grid[0]
                 delta_y=y_grid[1]-y_grid[0]
                 
                 kx_e=kx; kx_w=0.0  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=ky; ky_s=0.0  #Permeabilidad en y evaluada en las caras N y S 
                 
            if((i==Nx-1)and(j==0)):  #Nodo de esquina (Nx-1,0)
                 deltax_w=x_centers[i]-x_centers[i-1]
                 deltax_e=deltax_w
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=deltay_n 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[1]-y_grid[0]
                 
                 kx_e=0.0; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=ky; ky_s=0.0  #Permeabilidad en y evaluada en las caras N y S 
                 
            if((i==0)and(j==Ny-1)):  #Nodo de esquina (0,Ny-1)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=deltax_e
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s
                 
                 delta_x=x_grid[1]-x_grid[0]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=0.0  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=0.0; ky_s=ky  #Permeabilidad en y evaluada en las caras N y S 
                 
            if((i==Nx-1)and(j==Ny-1)):  #Nodo de esquina (Nx-1,Ny-1)
                 deltax_w=x_centers[i]-x_centers[i-1]
                 deltax_e=deltax_w
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=0.0; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=0.0; ky_s=ky  #Permeabilidad en y evaluada en las caras N y S 
                 
            Sw_e, Sw_w, Sw_n, Sw_s=pp.upwind_2D(Nx, Ny, i, j, P, Sw)   #Funciones que llaman a otro archivo
            krw_e, krw_w, krw_n, krw_s=pp.krw_2D(Sw_e, Sw_w, Sw_n, Sw_s, Srw, Sro, ommega)    
            kro_e, kro_w, kro_n, kro_s=pp.kro_2D(Sw_e, Sw_w, Sw_n, Sw_s, Srw, Sro, ommega)
            
            lamdaw_w=krw_w/Muw ;  lamdaw_e=krw_e/Muw; lamdaw_n=krw_n/Muw; lamdaw_s=krw_s/Muw
            
            lamdao_w=kro_w/Muo; lamdao_e=kro_e/Muo; lamdao_n=kro_n/Muo;  lamdao_s=kro_s/Muo
            
            lamdaT_w=lamdaw_w+lamdao_w; lamdaT_e=lamdaw_e+lamdao_e  
            lamdaT_n=lamdaw_n+lamdao_n; lamdaT_s=lamdaw_s+lamdao_s
            
            lamdaw=pp.krw_point(Sw[i][j], Srw, Sro, ommega)/Muw
            lamdao=pp.kro_point(Sw[i][j], Srw, Sro, ommega)/Muo
            
            re=0.14*((ky/kx)**0.5*delta_x**2.0+(kx/ky)**0.5*delta_y**2.0)**0.5/(0.5*((ky/kx)**0.25+(kx/ky)**0.25))
            
            if ((coordInj_i==i ) and (coordInj_j==j)):
                WI_inj=(2*np.pi*delta_z*(kx*ky)**0.5)/(np.log(re/rwell_inj))
            else:
                WI_inj=0.0
            if ((coordProd_i==i) and (coordProd_j==j)):
                WI_prod=(2*np.pi*delta_z*(kx*ky)**0.5)/(np.log(re/rwell_prod))
            else:
                WI_prod=0.0
                    
            AP[i][j]=((kx_e*lamdaT_e/deltax_e+(kx_w*lamdaT_w)/deltax_w)*delta_y*delta_z)+((ky_n*lamdaT_n/deltay_n+ky_s*lamdaT_s/deltay_s)*delta_x*delta_z) 
            AE[i][j]=(kx_e*lamdaT_e/deltax_e)*delta_y*delta_z 
            AW[i][j]=(kx_w*lamdaT_w/deltax_w)*delta_y*delta_z
            AN[i][j]=(ky_n*lamdaT_n/deltay_n)*delta_x*delta_z 
            AS[i][j]=(ky_s*lamdaT_s/deltay_s)*delta_x*delta_z 
            B[i][j]=lamdaw*WI_inj*Pbh_inj+lamdaw*WI_prod*Pbh_prod+lamdao*WI_prod*Pbh_prod
               
            
                 
             
                    
           
                



     
def coefSat_2D(Nx, Ny, delta_t, x_centers, y_centers, x_grid, y_grid, delta_z, rwell_inj, coordInj_i, coordInj_j, pbh_inj, Sw_inj, rwell_prod, coordProd_i, coordProd_j, pbh_prod, Sw_prod, phi, Srw, Sro, ommega, Muw, Muo, kx, ky, P, Sw, Sw_old, AP, AE, AW, AN, AS, B):

    for i in range (0,Nx):    
        for j in range (0,Ny):
            
            if((i>0)and(i<Nx-1)and(j>0)and(j<Ny-1)):       #Distancias de volumenes y centros nodos internos          
                deltax_e=x_centers[i+1]-x_centers[i]      #Distancia del nodo P al nodo E
                deltax_w=x_centers[i]-x_centers[i-1]      #Distancia del nodo W al nodo P
                deltay_n=y_centers[j+1]-y_centers[j]      #Distancia del nodo P al nodo N
                deltay_s=y_centers[j]-y_centers[j-1]      #Distancia del nodo S al nodo P
                
                delta_x=x_grid[i]-x_grid[i-1]
                delta_y=y_grid[j]-y_grid[j-1]
                
                kx_e=kx; kx_w=kx   #Permeabilidad en x evaluada en las caras E y W
                ky_n=ky; ky_s=ky   #Permeabilidad en y evaluada en las caras N y S
            
          
## condiciones en las frontera cerradas (no flujo)
    
            if((i==0)and(j>0)and(j<Ny-1)):  #Distancias de volumenes y centros en frontera oeste (w)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=deltax_e   #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=y_centers[j]-y_centers[j-1]
                 
                 delta_x=x_grid[1]-x_grid[0]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=0.0  #Permeabilidad en x evaluada en las caras E y W (Frontera cerrada oeste)
                 ky_n=ky; ky_s=ky   #Permeabilidad en y evaluada en las caras N y S
            
            if((i==Nx-1)and(j>0)and(j<Ny-1)):  #Distancias de volumenes y centros en frontera este (e)
                 deltax_w=x_centers[i]-x_centers[i-1]
                 deltax_e=deltax_w   #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=y_centers[j]-y_centers[j-1]
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=0.0; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W (Frontera cerrada este)
                 ky_n=ky; ky_s=ky   #Permeabilidad en y evaluada en las caras N y S
                 
            if((i>0)and(i<Nx-1)and(j==0)):  #Distancias de volumenes y centros en frontera sur (s)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=x_centers[i]-x_centers[i-1] 
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=deltay_n #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[1]-y_grid[0]
                 
                 kx_e=kx; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W 
                 ky_n=ky; ky_s=0.0  #Permeabilidad en y evaluada en las caras N y S (Frontera cerrada sur)
                 
            if((i>0)and(i<Nx-1)and(j==Ny-1)):  #Distancias de volumenes y centros en frontera norte (n)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=x_centers[i]-x_centers[i-1] 
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W 
                 ky_n=0.0; ky_s=ky #Permeabilidad en y evaluada en las caras N y S (Frontera cerrada norte)
                 
            if((i>0)and(i<Nx-1)and(j==Ny-1)):  #Distancias de volumenes y centros en frontera norte (n)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=x_centers[i]-x_centers[i-1] 
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s #Se cliona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W 
                 ky_n=0.0; ky_s=ky #Permeabilidad en y evaluada en las caras N y S (Frontera cerrada norte)
           
            if((i==0)and(j==0)):  #Nodo de esquina (0,0)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=deltax_e
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=deltay_n #Se clona distancia hacia nodo fantasma y se inckuyte para no tener una division por cero. 
                 
                 delta_x=x_grid[1]-x_grid[0]
                 delta_y=y_grid[1]-y_grid[0]
                 
                 kx_e=kx; kx_w=0.0  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=ky; ky_s=0.0  #Permeabilidad en y evaluada en las caras N y S 
                 
            if((i==Nx-1)and(j==0)):  #Nodo de esquina (Nx-1,0)
                 deltax_w=x_centers[i]-x_centers[i-1]
                 deltax_e=deltax_w
                 deltay_n=y_centers[j+1]-y_centers[j]
                 deltay_s=deltay_n 
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[1]-y_grid[0]
                 
                 kx_e=0.0; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=ky; ky_s=0.0  #Permeabilidad en y evaluada en las caras N y S 
                 
            if((i==0)and(j==Ny-1)):  #Nodo de esquina (0,Ny-1)
                 deltax_e=x_centers[i+1]-x_centers[i]
                 deltax_w=deltax_e
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s
                 
                 delta_x=x_grid[1]-x_grid[0]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=kx; kx_w=0.0  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=0.0; ky_s=ky  #Permeabilidad en y evaluada en las caras N y S 
                 
            if((i==Nx-1)and(j==Ny-1)):  #Nodo de esquina (Nx-1,Ny-1)
                 deltax_w=x_centers[i]-x_centers[i-1]
                 deltax_e=deltax_w
                 deltay_s=y_centers[j]-y_centers[j-1]
                 deltay_n=deltay_s
                 
                 delta_x=x_grid[i]-x_grid[i-1]
                 delta_y=y_grid[j]-y_grid[j-1]
                 
                 kx_e=0.0; kx_w=kx  #Permeabilidad en x evaluada en las caras E y W (Caso especial Frontera cerrada)
                 ky_n=0.0; ky_s=ky  #Permeabilidad en y evaluada en las caras N y S 
            delta_v = delta_x*delta_y*delta_z
            
            Sw_e, Sw_w , Sw_n , Sw_s = fP.upwind_2D(Nx, Ny, i, j, P, Sw)
            krw_e,krw_w,krw_n,krw_s = fP.krw_2D(Sw_e,Sw_w,Sw_n,Sw_s,Srw,Sro,ommega)
            
            lamdaw_e=krw_e/Muw;  lamdaw_w=krw_w/Muw;  lamdaw_n=krw_n/Muw; lamdaw_s=krw_s/Muw;
            
            
            lamdaw = fP.krw_point(Sw_old[i][j],Srw,Sro,ommega)/Muw
            re = 0.14*((ky/kx)**0.5*delta_x**2.0+(kx/ky)**0.5*delta_y**2.0)**0.5/(0.5*((ky/kx)**0.25+(kx/ky)**0.25))
            
            if((coordInj_i==i)and(coordInj_j==j)):
                WI_inj = (2*np.pi*delta_z*(kx*ky**0.5))/(np.log(re/rwell_inj))
            else:
                WI_inj = 0
            if((coordProd_i==i)and(coordProd_j==j)):
                WI_prod = (2*np.pi*delta_z*(kx*ky**0.5))/(np.log(re/rwell_prod))
            else:
                WI_prod = 0
                
            
            AP[i][j] = ((kx_e*lamdaw_e)/deltax_e+(kx_w*lamdaw_w)/deltax_w)*(delta_t/(delta_x*phi))+((ky_n*lamdaw_n)/deltay_n+(ky_s*lamdaw_s)/deltay_s)*(delta_t/(delta_y*phi))+lamdaw*WI_inj*(delta_t/(phi))+lamdaw*WI_prod*(delta_t/phi)
            AE[i][j] = ((kx_e*lamdaw_e)/deltax_e)*(delta_t/(delta_x*phi))
            AW[i][j] = ((kx_w*lamdaw_w)/deltax_w)*(delta_t/(delta_x*phi))
            AN[i][j] = ((ky_n*lamdaw_n)/deltay_n)*(delta_t/(delta_y*phi))
            AS[i][j] = ((ky_s*lamdaw_s)/deltay_s)*(delta_t/(delta_y*phi))
            
            B[i][j] = Sw_old[i][j] + lamdaw*WI_inj*pbh_inj*(delta_t/(phi))+lamdaw*WI_prod*pbh_prod *(delta_t/(phi))
            
            
            
