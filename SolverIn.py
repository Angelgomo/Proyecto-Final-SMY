import numpy as np



"""
El algortimo de Thomas es un método directo para la solución de matrices tridiagonales, es el más 
utilizado en este rubro, además se puede extender con facilidad a matrices pentadiagonales y heptadiagonales
las cuales resultan de aplicar la discretización de diferencias finitas o volumen finito a la ecuaciones diferenciales elipticas.  
"""

def thomas1D(AP, AE, AW, B):
    n = len(B)
   #print (n)
    x = np.zeros(n)
    P = np.zeros(n)                             # Se crean los vectores P y Q para las funciones de recurrencia
    Q = np.zeros(n)

    P[0] = AE[0]/AP[0]    #formula 9 
    Q[0] = B[0]/AP[0]     #formula 10
    
    
    for i in range(1,n,1):	#El error estaba en el "for" donde barre debe ser "n" para que barra los nodos 1 hasta n-1  
        P[i] = AE[i]/(AP[i]-AW[i]*P[i-1])					# Se calcula la función de recurrencia P, formula (7)
        Q[i] = (B[i]+AW[i]*Q[i-1])/(AP[i]-AW[i]*P[i-1])		# Se calcula la función de recurrencia Q, formula (8)
        #print(i)	
    
    x[n-1] = Q[n-1]   #formula 12
	
    i = n-2
    while i !=-1:
        x[i]=P[i]*x[i+1]+Q[i]	# Sustitución regresiva (formula 3)
        i-=1	
    return (x)   




def TDMA2Dadi(NX, NY, AP, AE, AW, AN, AS, B, FHI, iteramax, eps):
	
	Px = np.zeros(NX)
	Qx = np.zeros(NX)
	Py = np.zeros(NY)
	Qy = np.zeros(NY)
	FHIOLD = np.zeros((NX,NY))
	Residual = np.zeros((NX,NY))
	maxres = 1.0
	iteracion = 1

	
	while (iteracion <= iteramax)and(maxres > eps): 	
		for i in range(0,NX):			
			for j in range(0,NY):
				FHIOLD[i][j] = FHI[i][j]
#		print FHIOLD
# Barrido en dirección X
	
		for j in range(0,NY):											
			if (j > 0)and(j<NY-1):											# Punto inicial
				beta = AN[0][j]*FHIOLD[0][j+1]+AS[0][j]*FHIOLD[0][j-1]+B[0][j]
				Px[0] = AE[0][j]/AP[0][j]			
				Qx[0] = beta/AP[0][j]
			if(j == 0):
				beta = AN[0][j]*FHIOLD[0][j+1]+B[0][j]
				Px[0] = AE[0][j]/AP[0][j]			
				Qx[0] = beta/AP[0][j]
			if(j == NY-1):
				beta = AS[0][j]*FHIOLD[0][j-1]+B[0][j]
				Px[0] = AE[0][j]/AP[0][j]			
				Qx[0] = beta/AP[0][j]
			
			for i in range(1,NX):
				if (j > 0)and(j<NY-1):											
					beta = AN[i][j]*FHIOLD[i][j+1]+AS[i][j]*FHIOLD[i][j-1]+B[i][j]
					Px[i] = AE[i][j]/(AP[i][j]-AW[i][j]*Px[i-1])
					Qx[i] = (beta+AW[i][j]*Qx[i-1])/(AP[i][j]-AW[i][j]*Px[i-1])		
				if(j == 0):
					beta = AN[i][j]*FHIOLD[i][j+1]+B[i][j]
					Px[i] = AE[i][j]/(AP[i][j]-AW[i][j]*Px[i-1])
					Qx[i] = (beta+AW[i][j]*Qx[i-1])/(AP[i][j]-AW[i][j]*Px[i-1])
				if(j == NY-1):
					beta = AS[i][j]*FHIOLD[i][j-1]+B[i][j]
					Px[i] = AE[i][j]/(AP[i][j]-AW[i][j]*Px[i-1])
					Qx[i] = (beta+AW[i][j]*Qx[i-1])/(AP[i][j]-AW[i][j]*Px[i-1])

			FHI[NX-1][j] = Qx[NX-1]    # Asignación del ultimo punto    
			
			i = NX-2
			while i !=-1:
				FHI[i][j]=Px[i]*FHI[i+1][j]+Qx[i]
				i-=1

		for i in range(0,NX):			
			for j in range(0,NY):
				FHIOLD[i][j] = FHI[i][j]

# Barrido en dirección Y
	
		for i in range(0,NX):											
			if (i > 0)and(i < NX-1):											# Punto inicial
				beta = AE[i][0]*FHIOLD[i+1][0]+AW[i][0]*FHIOLD[i-1][0]+B[i][0]
				Py[0] = AN[i][0]/AP[i][0]			
				Qy[0] = beta/AP[i][0]
			if(i == 0):
				beta = AE[i][0]*FHIOLD[i+1][0]+B[i][0]
				Py[0] = AN[i][0]/AP[i][0]			
				Qy[0] = beta/AP[i][0]
			if(i == NX-1):
				beta = AW[i][0]*FHIOLD[i-1][0]+B[i][0]
				Py[0] = AN[i][0]/AP[i][0]			
				Qy[0] = beta/AP[i][0]
			
			for j in range(1,NY):
				if (i > 0)and(i < NX-1):										
					beta = AE[i][j]*FHIOLD[i+1][j]+AW[i][j]*FHIOLD[i-1][j]+B[i][j]
					Py[j] = AN[i][j]/(AP[i][j]-AS[i][j]*Py[j-1])
					Qy[j] = (beta+AS[i][j]*Qy[j-1])/(AP[i][j]-AS[i][j]*Py[j-1])		
				if(i == 0):
					beta = AE[i][j]*FHIOLD[i+1][j]+B[i][j]
					Py[j] = AN[i][j]/(AP[i][j]-AS[i][j]*Py[j-1])
					Qy[j] = (beta+AS[i][j]*Qy[j-1])/(AP[i][j]-AS[i][j]*Py[j-1])		
				if(i == NX-1):
					beta = AW[i][j]*FHIOLD[i-1][j]+B[i][j]
					Py[j] = AN[i][j]/(AP[i][j]-AS[i][j]*Py[j-1])
					Qy[j] = (beta+AS[i][j]*Qy[j-1])/(AP[i][j]-AS[i][j]*Py[j-1])		
	
			FHI[i][NY-1] = Qy[NY-1]    # Asignación del ultimo punto    
			
			j = NY-2
			while j !=-1:
				FHI[i][j]=Py[j]*FHI[i][j+1]+Qy[j]
				j-=1	

#----------------------- Calculo del residual --------------------------#
		maxresold = 0.0

		for i in range(1,NX-1):			
			for j in range(1,NY-1):
				Residual[i][j]=abs(AP[i][j]*FHI[i][j]-(AE[i][j]*FHI[i+1][j]+AW[i][j]*FHI[i-1][j]+AN[i][j]*FHI[i][j+1]+AS[i][j]*FHI[i][j-1])-B[i][j])
				maxres = max(Residual[i][j],maxresold)
				maxresold = maxres
#		maxres = max(Residual)

#		print 'iteracion = %d, Residual = %e' % (iteracion, maxres)		
		iteracion+=1
	return (FHI)
