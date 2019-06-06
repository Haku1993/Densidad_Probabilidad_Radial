import scipy as sc 
import math  as mt
import sympy as sp

n=3 # numero cuantico n
l=2 # numero cuantico l

r=sp.symbols('r')
x=sp.symbols('x')

a0=1
a=2/(n*a0) #alpha
j=2*n*(n+l)
m=n-l-1
k=2*l+1

def L(m): #polinomio de Laguerre
  f=sp.exp(-x)*pow(x,m+k) 
  ff=sp.diff(f,x,m) #Derivada m-sima
  L=sp.simplify(sp.exp(x)*pow(x,-k)/mt.factorial(m)*ff) #simplifica la expresion   
  return L.subs(x,a*r)

def F(m,a):
  F=pow(a**3*mt.factorial(m)/mt.factorial(j),1/2)*sp.exp(-a*r/2)*pow(a,l)*pow(r,l+1)
  return F

def Ur(r):
  ur=F(m,a)*L(m)
  ur=sp.simplify(ur)
  return ur

err2s=5*10**(3)
err2p=0.43*10**(7)
err3s=0.15*10**(15)
err3p=0.04*10**(23)
err3d=0.03*10**(31)

for i in range(500):
    r=0.04*i
    print(r," ",err3d*mt.fabs(Ur(r))**2)




