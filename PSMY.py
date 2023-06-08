import matplotlib.pyplot as plt     #importa la biblioteca matplotlib con el sobrenombre plt (gráficar) 
import SolverIn as SL              #importa el archivo propio solverLin con el sobrenombre SL (resolver sistemas tridiagonales) 
import numpy as np                  #importa la biblioteca numpy para utilizar matrices y vectores (numerical python)
import CoeficientesPyS as coeff     #importa el archivo propio coefficients_BL con el sobrenombre  coeff (coeficientes de las ecuaciones de presión y saturación AP, AE, AW, B)
from PIL import Image              #De la biblioteca PIL importa el objeto Image (Nos permite abrir imagenes para manipularlas)
import pandas as pd                #Importa la biblioteca Pandas para utilizar bases de datos  
from openpyxl import load_workbook   #Del motor de bases de datos "openpyxl" importa load worbook para cargar o leer archivos de datos existentes 


"""
AUTOR: Dr. Victor Leonardo Teja Juárez
INSTITUCIÓN: UNAM, Facultad de Ingeniería, División de Ciencias de la Tierra, Ingeniería Petrolera.
DESCRIPCIÓN:
En este código se resuelve el problema de flujo bifasico incompresible en medios porosos en     
dos dimensión, tomando como referencia el problema propuesto por Buckley-Leverett.

Se utilizó la formulación de flujo fraccional para simplificar la solución del problema en una 
ecuación para la presión y otra para la saturación de agua. Cabe aclarar que el problema se resuelve
con el algoritmo IMPES, es decir se toma una discretización ímplicita en el tiempo para la presión
y una discretización explicita para la saturación de agua. (VER DIAPOSITIVAS SMY 10 y 11)
"""
#-------------------------------------------------------------------------#
#----DEFINICIÓN DE LAS PROPIEDAES FÍSICAS DEL DOMINIO Y DE LOS FLUIDOS----# 
#-------------------------------------------------------------------------#

Lx=300  #Longitud del dominio (m)
Ly=300
kx=1E-15 # Permeabilidad absoluta(m^2) 
ky=2E-15


phi=0.2  #Porosidad
Muw=1E-03 #Viscosidad del agua (Pa*s)
Muo=1E-03 #Viscosidad del aceite (Pa*s)
Srw=0.0   #Saturación residual de agua
Sro=0.2   #Saturación residual del aceite
Sw_ini=0.0  #Saturación inicial de agua en el medio poroso
vel_inj= 3.4722E-07   #Velocidad de inyección m/s
Sw_inj=0.8          #Saturación de inyección
Sw_prod=0.0          #Saturación de inyección
press_prod=10E06    #presión de producción en Pa (10 MPa) 
press_ini=11E06    #presión de inicial en Pa (10 MPa) 
ommega=1.0      # Exponente para calcular la permeabilidad relativa
#theta=2.0
delta_y=1.0   #El Área de la sección transversal se considera unitaria por falta de datos del paper
delta_z=1.0

#-------------------------------------------------------------------------#
#--------- DEFINICIÓN DEL TIEMPO DE SIMULACIÓN Y PASO  DE TIEMPO --------# 
#-------------------------------------------------------------------------#4

tiempo_total=2.0   #Días 
delta_t=0.1      #Días 

delta_save=1.0   #Paso de tiempo para salvar imagenes en días
num_save=int(tiempo_total/delta_save)  #Calculo del numero de veces que salvamos o guardamos imagenes 


#-------------------------------------------------------------------------#
#------------------------- FACTORES DE CONVERSIÓN --------# 
#-------------------------------------------------------------------------#
dia_a_seg=86400.0
delta_t=delta_t*dia_a_seg

#-------------------------------------------------------------------------#
#------ DEFINICIÓN DEL LA MALLA Y SUS CENTROS (CELDAS CENTRADAS)----------# 
#-------------------------------------------------------------------------#
Nx=10                               #Numero de celdas en x
Ny=10


x_grid=np.linspace(0,Lx,Nx+1)       #Creación de la malla vertices
delta_x=Lx/Nx                       #Ancho de la celda
x_centers=np.linspace(delta_x/2.0,Lx-delta_x/2.0,Nx)  #Creación de los centros

y_grid=np.linspace(0,Ly,Ny+1)       #Creación de la malla vertices
delta_y=Ly/Ny                       #Ancho de la celda
y_centers=np.linspace(delta_y/2.0,Ly-delta_y/2.0,Ny)  #Creación de los centros

#definicion de los datos de los pozos de inyeccion y de produccion
rwell_inj=0.1016  #radio del pozo de inyeccion metros 
coordInj_i=0; coordInj_j=0 #coordenadas del pozo de inyeccion
pbh_inj=15e06 #presion de inyeccion MPA

rwell_prod=0.1016  #radio del pozo de produccion metros 
coordProd_i=Nx-1; coordProd_j=Ny-1 #coordenadas del pozo de produccion
pbh_prod=1e06 #presion de produccion MPA


#----------------------------------------------------------------------------#
#-DEFINICIÓN DE LOS VALORES INICIALES PARA LAS VARIABLES PRIMARIAS P Y Sw----# 
#----------------------------------------------------------------------------#
P=np.ones((Nx,Ny))*press_ini    #Presión del tiempo actual
Sw=np.ones((Nx,Ny))*Sw_ini      #Saturación del tiempo actual

Sw[coordInj_i-1][coordInj_j-1]=Sw_inj

P_old=np.ones((Nx,Ny))*press_ini   # P tiempo anterior
Sw_old=np.ones((Nx,Ny))*Sw_ini     # Sw tiempo anteior

#----------------------------------------------------------------------------#
# DEFINICIÓN DE LOS COEFICIENTES PARA RESOLVER LAS ECUACIONES DE P Y Sw -----# 
#----------------------------------------------------------------------------#
AP=np.zeros((Nx,Ny))
AE=np.zeros((Nx,Ny))
AW=np.zeros((Nx,Ny))
AS=np.zeros((Nx,Ny))
AN=np.zeros((Nx,Ny))
B=np.zeros((Nx,Ny))

#----------------------------------------------------------------------------#
#- DEFINICIÓN DE LISTA Y CONTADOR PARA SALVAR LAS IMAGENES O FRAMES  --------# 
#----------------------------------------------------------------------------#

#Lista vacia para salvar las imagenes (se hará tan grande como tantas imagenes guardemos en ella)
framesSat=[]     #Salvar las imagenes del perfil de saturación
framesPress=[]  #Salvar las imagenes del perfil de presión
count_save=1    #Contador para salvar imagenes
Excel_name='ResultsBL_2D.xlsx'    #Nombre del archivo de excel a crear

#x_centV=np.zeros(Nx,Ny);x_centV=np.zeros(Nx,Ny);

"""
|*****************************************************************************|
|---------------------------INICIA EL ALGORITMO IMPES-------------------------|
|*****************************************************************************|
"""
tiempo=0.0
eps=1e-7
iteramax = 2000

while (tiempo<tiempo_total):  #Paso No 2 del algoritmo
    
    ##Paso No 3: del algoritmo Calculo de coeficientes para la ecuación de la presión (ecu 33 a 36 de Diapositivas)
    AP=np.zeros((Nx,Ny)); AE=np.zeros((Nx,Ny)); AW=np.zeros((Nx,Ny)); AS=np.zeros((Nx,Ny)); AN=np.zeros((Nx,Ny)) ; B=np.zeros((Nx,Ny)) #Se mandan a cero para ser reutilizados para la ecuación de presión
    #coeff.coeffPress_1DBL(Nx,Ny,x_centers,y_centers,x_grid,y_grid, delta_y, delta_z,rwell_inj,coordInj_i,coordInj_j,pbh_inj, vel_inj, Sw_inj,rwell_prod,coordProd_i,coordProd_j,pbh_prod, press_prod, Sw_prod, kx,ky, Muw, Muo, Srw, Sro, ommega, Sw, P, AP, AE, AW,AN,AS, B)
    coeff.coefPres_2D(Nx, Ny, x_centers, y_centers, x_grid, y_grid, delta_z, rwell_inj, coordInj_i, coordInj_j, pbh_inj, Sw_inj, rwell_prod, coordProd_i, coordProd_j, pbh_prod, Sw_prod, Srw, Sro, ommega, Muw, Muo, kx, ky, P, Sw, AP, AE, AW, AN, AS, B)
    P=SL.TDMA2Dadi(Nx,Ny,AP,AE,AW,AN,AS,B,P,iteramax,eps)
    
    ##Paso No 5: Calculo de coeficientes para la ecuación de saturación de agua (ecu 40 a 43 de Diapositivas)    
    AP=np.zeros((Nx,Ny)); AE=np.zeros((Nx,Ny)); AW=np.zeros((Nx,Ny)); AS=np.zeros((Nx,Ny)); AN=np.zeros((Nx,Ny)) ; B=np.zeros((Nx,Ny)) #Se mandan a cero para ser reutilizados para la ecuación de presión
   #coeff.coeffSat_1DBL(Nx,Ny,delta_t, delta_y, delta_z, x_grid,y_grid, x_centers,y_centers, phi, vel_inj, Sw_inj, press_prod, Sw_prod, kx, Muw, Muo, Srw, Sro, ommega, Sw, P, Sw_old, AP, AE, AW, B)
    coeff.coefSat_2D(Nx, Ny, delta_t, x_centers, y_centers, x_grid, y_grid, delta_z, rwell_inj, coordInj_i, coordInj_j, pbh_inj, Sw_inj, rwell_prod, coordProd_i, coordProd_j, pbh_prod, Sw_prod, phi, Srw, Sro, ommega, Muw, Muo, kx, ky, P, Sw, Sw_old, AP, AE, AW, AN, AS, B)
    #Paso No 6: Calculo de la saturación de agua de manera explicita 
    
    for j in range (0,Ny):
        for i in range(0,Nx):
            
            if ((i>0)and(i<Nx-1)and
                (j>0)and(j<Ny-1)):
                
                Sw[i][j]=-AP[i][j]*P[i][j]+AE[i][j]*P[i+1][j]+AW[i][j]*P[i-1][j]+AN[i][j]*P[i][j+1]+AS[i][j]*P[i][j-1]+B[i][j]
            
            if((i==0)and(j>0)and(j<Ny-1)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AE[i][j]*P[i+1][j]+AN[i][j]*P[i][j+1]+AS[i][j]*P[i][j-1]+B[i][j]
            if((i==Nx-1)and(j>0)and(j<Ny-1)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AW[i][j]*P[i-1][j]+AN[i][j]*P[i][j+1]+AS[i][j]*P[i][j-1]+B[i][j]
            if((i>0)and(i<Nx-1)and(j==0)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AE[i][j]*P[i+1][j]+AW[i][j]*P[i-1][j]+AN[i][j]*P[i][j+1]+B[i][j]
            if((i>0)and(i<Nx-1)and(j==Ny-1)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AE[i][j]*P[i+1][j]+AW[i][j]*P[i-1][j]+AS[i][j]*P[i][j-1]+B[i][j]
            
            if((i==0)and(j==0)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AE[i][j]*P[i+1][j]+AN[i][j]*P[i][j+1]+B[i][j]
            if((i==Nx-1)and(j==0)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AW[i][j]*P[i-1][j]+AN[i][j]*P[i][j+1]+B[i][j]
            if((i==0)and(j==Ny-1)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AE[i][j]*P[i+1][j]+AS[i][j]*P[i][j-1]+B[i][j]
            if((i==Nx-1)and(j==Ny-1)):
                Sw[i][j]=-AP[i][j]*P[i][j]+AW[i][j]*P[i-1][j]+AS[i][j]*P[i][j-1]+B[i][j]
                
            if((i==coordInj_i)and(j==coordInj_j)):
                Sw[i][j]=Sw_inj
                
            # x_centV[i+j*Nx]=x_centers[i]
            # y_centV[i+j*Nx]=y_centers[j]
            # P_V[i+j*Nx]=P[i][j]
            # Sw_V[i+j*Nx]=Sw[i][j]


    P_old=np.copy(P)   #Actualización de Presión    
    Sw_old=np.copy(Sw)   #Actualización de la saturación del tiempo anterior a la saturación actual (IMPORTANTE EN EL ALGORITMO)


    tiempo=tiempo+delta_t/dia_a_seg    #Paso 7 del algortimo    
    print("Días de simulación:", tiempo)
        
    
    #--------------------------------------------------------------------#
    #--------- GUARDADO DE ARCHIVOS Y GRAFICAS CADA "DELTA_SAVE"---------#
    #--------------------------------------------------------------------#


    if((abs(tiempo-delta_save*count_save)<1e-05)and(count_save<=num_save)):    #Condicional para el salvado de graficas (se puede utilizar también para salvar los resultados en tablas)
        #Variables tipo cadena para salvar las graficas de saturación
        stringTimes=round(tiempo,2)      #Variable que redondea el tiempo actual en dos digitos despues del punto
        stringSat = 'Saturacion %s dias' %str(stringTimes)   #Variable tipo cadena que indica la saturación y los días de simulación se le agrega stringtimes 
        stringSatSave = 'Saturacion %s dias.png' %str(stringTimes) #Variable tipo cadena que indica la saturación y los días de simulación se le agrega stringtimes servirá para salvar las imagenes en formato png
        
        plt.figure('Perfil de Sw', figsize=(8,8))   #Instrucciones para graficar el perfil de saturación vs centros de las celdas 
        plt.title(stringSat)
        plt.contourf(x_centers,y_centers,Sw.transpose(),10,alpha=0.75) 
        #plt.plot(x_centers,1-Sro-Sw) #Ambos frentes Sw y So 
        plt.xlabel("$ x(m) $");plt.ylabel("$ y(m) $");plt.grid()
        #plt.ylabel("Saturación de agua $S_{w}$")
        plt.grid()
        plt.savefig(stringSatSave, dpi=300)   #Se guarda el archivo o frame de imagene de la saturación
        #plt.show()   #Descomentar para ver en tiempo de ejecución
        
        imgSat = Image.open(stringSatSave)    #Se abre el archivo de imagen de la saturación recien guardado
        framesSat.append(imgSat)            #El archivo abierto de imagen de anexa a la lista de frames

        #Variables tipo cadena para salvar las graficas de presión
        stringPress = 'Presion %s dias' %str(stringTimes)   #Variable tipo cadena que indica la presión y los días de simulación se le agrega stringtimes 
        stringPressSave = 'Presion %s dias.png' %str(stringTimes) #Variable tipo cadena que indica la presión y los días de simulación se le agrega stringtimes servirá para salvar las imagenes en formato png
        
        plt.figure('Perfil de Po', figsize=(8,8))   #Instrucciones para graficar el perfil de presión vs centros de las celdas
        plt.title(stringPress)
        plt.contourf(x_centers,y_centers,P.transpose(),10,alpha=0.75)
        plt.colorbar()
        plt.xlabel("$ x(m) $");plt.ylabel("$ y(m) $");plt.grid()
        plt.savefig(stringPressSave, dpi=300)  #Se guarda el archivo o frame de imagen de la presion
        plt.show()   #Descomentar para ver en tiempo de ejecución
        
        imgPress = Image.open(stringPressSave)   #Se abre el archivo de imagen de la presion recien guardado
        framesPress.append(imgPress)   #El archivo abierto de imagen de anexa a la lista de frames

        #Variables tipo cadena para salvar los datos en un archivo de excel
        stringResSave='Resultados %s dias' %str(stringTimes)
        data = {'x (m)': x_centers[::], 'saturacion (Sw)': Sw[::], 'Presion P (Pa)': P[::]}  #Se crea el arreglo de datos con los centros de los volumenes, saturación de agua y presión
        df = pd.DataFrame(data, columns = ['x (m)', 'saturacion (Sw)', 'Presion P (Pa)'])   #Se le asigna el nombre a las columnas del dataFrame
        
        if (count_save==1):    
            with pd.ExcelWriter(Excel_name, engine='xlsxwriter') as writer:    
                df.to_excel(writer, stringResSave)  #Escribe en el archivo de excel y pone la primera hoja   
                writer.save() #Salva los datos         
       
        if (count_save>1):
            book = load_workbook(Excel_name)
            with pd.ExcelWriter(Excel_name, engine='openpyxl') as writer:
                writer.book = book  #Escribir en el libro de nombre Excel_name
                writer.sheets = dict((ws.title, ws) for ws in book.worksheets)  #Auxiliar para escribir las hojas subsecuentes y no sobreescribirlas     
                df.to_excel(writer, stringResSave)    #Escribe los datos del tiempo y simulación actual 
                writer.save() #Salva los datos
            
        count_save+=1
        
        
        
"""
|*****************************************************************************|
|---------------------------FIN DEL ALGORITMO IMPES-------------------------|
|*****************************************************************************|
"""
# #writer.close()  #Cierra el archivo de excel abierto para guardar los datos 

# framesPress[0].save('presionBL.gif', format='GIF', append_images=framesPress[1:], save_all=True, duration=200, loop=0)  #Crea el gif de la presión con las imagenes guardadas en la lista de frames
# framesSat[0].save('saturacionBL.gif', format='GIF', append_images=framesSat[1:], save_all=True, duration=200, loop=0)   #Crea el gif de la saturación  con las imagenes guardadas en la lista de frames
# SMY_15_BL.py
# Mostrando SMY_15_BL.py

