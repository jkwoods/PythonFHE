import numpy as np
import random
import sympy
from functools import reduce 
import math
from fractions import Fraction
from itertools import islice

def QuotientNear(a,b):
  "Gives the nearest integer to a/b"
  return (2*a+b)//(2*b)

def modNear(a,b):
  "Computes a mod b with a \in ]-b/2,b/2]"
  return a-b*QuotientNear(a,b)

def mod(c, p):
  return c % p

def random_prime(a,b):
  return sympy.randprime(a,b)

def random_element(a,b): #TODO hook up with seed??
  return random.randint(a,b)

def set_random_seed(seed): #if seed = 0, use rand and return it, else use seed
  s = seed
  if (seed == 0):
    random.seed()
    s = random.randint(2, 2**30) #TODO
    random.seed(s)
  else:
    random.seed(s)
  return s

def sumBinary(a,b):
  "Computes the sum of the binary vectors a and b, modulo 2^n where n is the length of the vectors a and b"
  c=[a[0]+b[0]]
  carry=a[0]*b[0]

  for i in range(1,len(a)-1):
    carry2=(a[i]+b[i])*carry+a[i]*b[i]    
    c.append(a[i]+b[i]+carry)             
    carry=carry2
  
  c.append(a[-1]+b[-1]+carry)
  return c

def toBinary(x,l):
  "Converts a positive integer x into binary with l digits"
  return digits(x+2**l)[:-1]

def digits(x): #always binary
  le = list('{0:0b}'.format(x))
  le.reverse()
  return le


def mul_inv(a, b):
    b0 = b
    x0, x1 = 0, 1
    if b == 1: return 1
    while a > 1:
        q = a // b
        a, b = b, a%b
        x0, x1 = x1 - q * x0, x0
    if x1 < 0: x1 += b0
    return x1

def CRT(n, a): #chinese remiander thm
    sum = 0
    prod = reduce(lambda a, b: a*b, n)
    for n_i, a_i in zip(n, a):
        p = prod // n_i
        sum += a_i * mul_inv(p, n_i) * p
    return sum % prod

def kd(i,j):
  if (i == j):
    return 1
  else:
    return 0

def arraymult(c,a):
  return [c*int(xi) for xi in a]

class Ciphertext():
  def __init__(self,val_,pk_):
    self.val,self.pk=val_,pk_

  @staticmethod
  def encrypt(pk,m):
    return Ciphertext(pk.encrypt(m),pk)

  def decrypt(self):
    m=self.pk.decrypt(self.val)
    return m

  def __add__(self,x):
    return self.__class__(self.pk.add(self.val,x.val),self.pk)

  def __mul__(self,x):
    return self.__class__(self.pk.mult(self.val,x.val),self.pk)

  def scalmult(self,x):
    if isinstance(x,list):
      print("multing list")

      return [self.__class__(self.val*int(xi),self.pk) for xi in x]
    else:
      return self.__class__(self.val*x,self.pk)

  def recrypt(self):
    return Ciphertext(self.pk.recrypt(self.val),self.pk)

def make_pri(x0,ell,seed): #generates X
    set_random_seed(seed)
    chi = [random_element(0,x0) for i in range(ell)]
    set_random_seed(0)
    return chi

def make_deltas(pk,lenv,rho,seed,cr):
  pr=make_pri(pk.x0,lenv,seed)
 
  r=[[random_element(-2**rho+1,2**rho) for i in range(pk.l)] for j in range(lenv)]
  E=[random_element(0,(2**(pk.lam+pk.log+(pk.l*pk.eta)))//pk.pi) for i in range(lenv)] #added from paper
  delta=[0 for i in range(lenv)]

  if (cr == 0):#x
    crts = [CRT(pk.p,[2*ri for ri in r[j]]) for j in range(lenv)]
  elif (cr == 1):#xi
    crts = [CRT(pk.p,[2*ri+kd(i,j) for ri,i in zip(r[j],range(pk.l))]) for j in range(lenv)]
  elif (cr == 2):#ii
    crts = [CRT(pk.p,[2*ri+(kd(i,j)*(2**(pk.rhoi+1))) for ri,i in zip(r[j],range(pk.l))]) for j in range(lenv)]
  else: #o
    crts = [CRT(pk.p,[2*ri+si for ri,si in zip(r[j],pk.verts[j])]) for j in range(lenv)]

  temp=[mod(Xi,pk.pi) for Xi in pr]
  delta=[te+(ei*pk.pi)-crti for te,ei,crti in zip(temp,E,crts)] #changed from paper
  return delta

def make_u_front(pk,seed):
  pr=make_pri(2**(pk.kap+1),pk.Theta,seed) #u draft
  u = pr

  n=0
  for j in range(pk.l):
    xpj = (2**pk.kap)//pk.p[j]

    su = [0 for i in range(pk.Theta)]
    for i in range(pk.Theta):
      su[i] = pk.s[j][i]*u[i]
    
    sumt = sum(su)
    sumt = mod(sumt, 2**(pk.kap+1))

    v = n
    n = n+1

    #change corresponding u
    su[v] = 0
    sumv = sum(su)
    k1 = 2**(pk.kap+1)
    nu = k1 - sumv + xpj
    while (nu < 0) or (nu >= k1):
      if (nu < 0):
        nu = nu+k1
      elif ():
        nu = nu-k1

    u[v] = nu

  return u[:pk.l]

class Pk(object):
  def __init__(self, key_size):
    self.lam = 12
    self.rho = 26 #p
    self.eta = 1988 #(n)
    self.gam = 147456 #y
    self.Theta = 150 #O
    self.alpha = 936
    self.tau = 188
    self.l = 10
    if (key_size==0):
        print("Making toy key")
        self.lam = 42
        self.rho = 26
        self.eta = 988
        self.gam = 290000
        self.Theta = 150
        self.alpha = 936
        self.tau = 188
        self.l = 10
    elif(key_size==1):
        print("making small key")
        self.lam = 52
        self.rho = 41
        self.eta = 1558
        self.gam = 1600000
        self.Theta = 555
        self.alpha = 1476
        self.tau = 661
        self.l = 37
    elif (key_size==2):
        print("making medium key")
        self.lam = 62
        self.rho = 56
        self.eta = 2128
        self.gam = 8500000
        self.Theta = 2070
        self.alpha = 2016
        self.tau = 2410
        self. l = 138
    elif (size == 3):
        self.lam = 72;
        self.rho = 71;
        self.eta = 2698;
        self.gam = 39000000;
        self.Theta = 7965;
        self.alpha = 2556;
        self.tau = 8713;
        self.l = 531;

 
    self.alphai = self.lam + self.alpha
    self.rhoi = self.lam + self.alpha
    self.n = 4;
    self.kap = 64*(self.gam//64+1)-1
    self.log = round(math.log2(self.l))
    self.theta = self.Theta//self.l

    self.p = [random_prime(2**(self.eta-1), 2**self.eta) for i in range(self.l)] #fix TODO ?????
    self.pi = reduce((lambda x, y: x * y), self.p) #product of p

    self.q0 = (2**self.gam)
    while (self.q0 > (2**self.gam)//self.pi):
      self.q0prime1 = random_prime(0, 2**(self.lam**2))
      self.q0prime2 = random_prime(0, 2**(self.lam**2))
      self.q0 = self.q0prime1*self.q0prime2
    self.x0=self.pi*self.q0
    
    self.x_seed = random.randint(2, 2**30)
    self.xi_seed = random.randint(2, 2**30)
    self.ii_seed = random.randint(2, 2**30)

    self.x_deltas = make_deltas(self,self.tau,self.rhoi-1,self.x_seed,0)
    self.xi_deltas = make_deltas(self,self.l,self.rho,self.xi_seed,1)
    self.ii_deltas = make_deltas(self,self.l,self.rho,self.ii_seed,2)
   
    self.B=self.Theta//self.theta

    self.s = [[0 for j in range(self.Theta)] for k in range(self.l)]

    for j in range(self.l):
      sj = []
      
      for t in range(self.theta):
        if (t==0): #if s[j][j] is in it
          fill = [0 for i in range(self.B)]
          fill[j] = 1
          sj = sj + fill
        else:
          fill = [0 for i in range(self.B)]
          sj = sj + fill
      self.s[j] = sj

    for t in range(1,self.theta):
      sri = random.sample(range(0, self.B), self.l)
      for j in range(self.l):
        k = (self.B*t)+sri[j]
        self.s[j][k] = 1

    self.verts = [[0 for j in range(self.l)] for k in range(self.Theta)]

    for i in range(self.Theta):
      for j in range(self.l):
        self.verts[i][j] = self.s[j][i]

    self.u_seed = random.randint(2, 2**30)
    self.o_seed = random.randint(2, 2**30)

    self.u_front = make_u_front(self, self.u_seed)

    self.o_deltas = make_deltas(self,self.Theta,self.rho,self.o_seed,3)

  def encrypt(self,m): #vector in {0,1}^l
    b = [random_element(-2**self.alpha,2**self.alpha) for i in range(self.tau)]
    bi= [random_element(-2**self.alphai,2**self.alphai) for i in range(self.l)]

    x = [c-d for c,d in zip(make_pri(self.x0,self.tau,self.x_seed),self.x_deltas)]
    xi= [c-d for c,d in zip(make_pri(self.x0,self.l,self.xi_seed),self.xi_deltas)]
    ii= [c-d for c,d in zip(make_pri(self.x0,self.l,self.ii_seed),self.ii_deltas)]

    sums=sum([mj*xij for mj,xij in zip(m,xi)])+sum([bij*iij for bij,iij in zip(bi,ii)])+sum([bj*xj for bj,xj in zip(b,x)])

    return modNear(sums,self.x0)

  def decrypt(self,c):
    return [mod(modNear(c,self.p[i]),2) for i in range(self.l)]

  def decrypt_squashed(self,c):
    y = [Fraction(ui)/Fraction(2**self.kap) for ui in self.u]
    z = [modNear(round((Fraction(c)*yi),4),2) for yi in y]

    m = [0 for i in range(self.l)]
    for j in range(self.l):
      sums = sum([sji*zi for sji,zi in zip(self.s[j],z)])
      sumsb = sum([mod(sji*zi,2) for sji,zi in zip(self.s[j],z)])
      
      m[j] = mod(int(round(sums)),2) ^ mod(c,2)
    return (m)

  def add(self,c1,c2):
    return mod(c1+c2,self.x0)

  def sub(self,c1,c2):
    return mod(c1-c2,self.x0)

  def mult(self,c1,c2):
    return mod(c1*c2,self.x0)

  def __repr__(self):
    return "<Pk with rho=%d, eta=%d, gam=%d>" % (self.rho,self.eta,self.gam) 

  def recrypt(self,c):
    #get u
    u_draft = make_pri(2**(self.kap+1),self.Theta,self.u_seed)
    u = self.u_front+u_draft[self.l:]

    #"expand"
    y = [Fraction(ui)/Fraction(2**self.kap) for ui in u]
    z = [mod(round((Fraction(c)*yi),4),2) for yi in y]#adjust bits of precision
    
    adjz = [int(round(zi*16)) for zi in z]

    #put z in binary arrays
    zbin = [toBinary(zi,self.n+1) for zi in adjz]
   

    #get o
    o = [c-d for c,d in zip(make_pri(self.x0,self.Theta,self.o_seed),self.o_deltas)]
   
    li = [arraymult(ski,cei) for ski,cei in zip(o,zbin)]

    Q_adds = [0 for i in range(self.n+1)]

    for t in range(self.Theta):
      Q_adds = sumBinary(Q_adds,li[t])
      Q_adds = [mod(qa,self.x0) for qa in Q_adds]

    rounded = Q_adds[-1] + Q_adds[-2] #"round"
 
    final = rounded + (c & 1)

    return final
