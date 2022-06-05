#Calculating energy eigenvalues

import numpy as npy
import matplotlib.pyplot as mplt

N=1000
I=300
m=1
g=1
a=1
hbar=1

#kinetic energy
T=npy.zeros((2*N+1)**2).reshape(2*N+1,2*N+1)
for i in range(2*N+1):
    for j in range(2*N+1):        #Define your Kinetic energy matrix here
        if i==j:
            T[i,j]=(i-N)**2*((hbar**2/I)/2)
        else:
            T[i,j]=0
            
#potential energy
V = npy.zeros((2*N+1)**2).reshape(2*N+1,2*N+1)
for i in range(2*N+1):
    for j in range(2*N+1):
        if i==j:
            V[i,j]=m*g*a
        elif npy.abs(i-j)==1:        #Define your Potential matrix here
            V[i,j]=-m*g*a/2
        else:
            V[i,j]=0

#Hamiltonian
H=T+V

#eigenvalues
val,vec=npy.linalg.eig(H)
z=npy.argsort(val)
z=z[0:300]

#output
print("1st   eigenvalue:\t",val[z][0])
print("50th  eigenvalue:\t",val[z][49])
print("200th eigenvalue:\t",val[z][199])

#plot
x=npy.linspace(0,299,300)
y=[]
y=npy.append(y,val[z])

mplt.figure(1)
mplt.plot(x,y,"r",lw=2)
mplt.title("Energy Eigenvalues",size=18)
mplt.xlabel("Energy Level",size=18)
mplt.ylabel("Energy",size=18)
mplt.axis([0,300,0,40])
mplt.grid(True)
mplt.show()

#difference between consequiutive degenerate states
steps=20
aloo=[]
k=41
k0=k
l=42

#for i in range(steps):
#    aloo=npy.append(aloo,0)
A=npy.zeros(steps)
for i in range(steps):
    A[i]=(val[z[l]]-val[z[k]])
    k=k+2
    l=l+2

x1=npy.linspace(k0,k0+steps-1,steps)
y1=[]
y1=npy.append(y1,A[:])
mplt.figure(2)
mplt.plot(x1,y1,"m",lw=2)
mplt.title("Energy Difference Between Two Consecutive States",size=18)
mplt.xlabel("Energy Level",size=18)
mplt.ylabel("Energy Difference $\Delta$E",size=18)
mplt.axis([40,60,0.0,0.03])
mplt.grid(True)
mplt.show()
