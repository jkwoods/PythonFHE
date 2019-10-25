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

def sum2(li):
  empty=True
  for x in li:
    if empty:
      res=x
      empty=False
    else:
      res=res+x
  if empty: 
    return 0
  else:
    return res
   

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

def xorBinary(a,b):
  return [ai+bi for ai,bi in zip(a,b)]

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
  def __init__(self,val_,pk_,degree_=1):
    self.val,self.pk,self.degree=val_,pk_,degree_

  @staticmethod
  def encrypt(pk,m):
    return Ciphertext(pk.encrypt(m),pk)

  def noise(self):
    return self.pk.noise(self.val)

  def decrypt(self):
    m=self.pk.decrypt(self.val)
    return m

  def decrypt_squashed(self):
    return self.pk.decrypt_squashed(self.val)

  def __add__(self,x):
    return self.__class__(self.pk.add(self.val,x.val),self.pk,max(self.degree,x.degree))

  def __mul__(self,x):
    return self.__class__(self.pk.mult(self.val,x.val),self.pk,self.degree+x.degree)

  def scalmult(self,x):
    if isinstance(x,list):
      print("multing list")

      return [self.__class__(self.val*int(xi),self.pk,self.degree) for xi in x]
    else:
      return self.__class__(self.val*x,self.pk,self.degree)

  def recrypt(self):
    return Ciphertext(self.pk.recrypt(self.val),self.pk)

class PRIntegers: #generates X
  "A list of pseudo-random integers."
  def __init__(self,x0,ell):
    self.x0,self.ell=x0,ell
    self.li=[None for i in range(self.ell)]
    self.se=set_random_seed(0)

  def __getitem__(self,i):
    return self.li[i]
    
  def __setitem__(self,i,val):
    self.li[i]=val

  def __iter__(self):
    set_random_seed(self.se)
    for i in range(self.ell):
      a=random_element(0,self.x0)
      if self.li[i]!=None:
        yield self.li[i]
      else:
        yield a
    set_random_seed(0)

class PRIntegersDelta(PRIntegers):
  """A list of pseudo-random integers, with their delta corrections"""
  def __iter__(self):
    return (c-d for c,d in zip(PRIntegers.__iter__(self),self.delta)) #changed from paper

  def ciphertexts(self,pk):
    return (Ciphertext(cd,pk) for cd in self)

  @staticmethod
  def encrypt(pk,v,rho,cr):
    pr=PRIntegersDelta(pk.x0,len(v))
    pr.r=[[random_element(-2**rho+1,2**rho) for i in range(pk.l)] for j in range(len(v))]
    E=[random_element(0,(2**(pk.lam+pk.log+(pk.l*pk.eta)))//pk.pi) for i in range(len(v))] #added from paper
    pr.delta=[0 for i in range(len(v))]

    if (cr == 0):#x
      crts = [CRT(pk.p,[2*ri for ri in pr.r[j]]) for j in range(len(v))]
    elif (cr == 1):#xi
      crts = [CRT(pk.p,[2*ri+kd(i,j) for ri,i in zip(pr.r[j],range(pk.l))]) for j in range(len(v))]
    elif (cr == 2):#ii
      crts = [CRT(pk.p,[2*ri+(kd(i,j)*(2**(pk.rhoi+1))) for ri,i in zip(pr.r[j],range(pk.l))]) for j in range(len(v))]
    else: #o
      
      crts = [CRT(pk.p,[2*ri+si for ri,si in zip(pr.r[j],pk.verts[j])]) for j in range(len(v))]

    temp=[mod(Xi,pk.pi) for Xi in PRIntegers.__iter__(pr)]
    pr.delta=[te+(ei*pk.pi)-crti for te,ei,crti in zip(temp,E,crts)] #changed from paper
    return pr

class PRIntegersU(PRIntegers):
  """A list of pseudo-random integers, frankencode crap"""
  def __iter__(self):
    return (PRIntegers.__iter__(self)) #changed from paper

  def ciphertexts(self,pk):
    return (Ciphertext(cd,pk) for cd in self)

  @staticmethod
  def encrypt(pk):
    pr=PRIntegersU(2**(pk.kap+1),pk.Theta) #u draft
    u = [ui for ui in PRIntegers.__iter__(pr)]

    for j in range(pk.l):
      s1indexs = []
      xpj = (2**pk.kap)//pk.p[j]

      su = [0 for i in range(pk.Theta)]
      for i in range(pk.Theta):
        su[i] = pk.s[j][i]*u[i]
        if (pk.s[j][i] == 1): s1indexs.append(i)

      sumt = sum2(su)
      sumt = mod(sumt, 2**(pk.kap+1))

      while(sumt != xpj):
        #pick random 1 in s
        v = random.choice(s1indexs)

        #change corresponding u
        su[v] = 0
        sumv = sum2(su)
        k1 = 2**(pk.kap+1)
        nu = k1 - sumv + xpj
        while (nu < 0) or (nu >= k1):
          if (nu < 0):
            nu = nu+k1
          elif ():
            nu = nu-k1

        u[v] = nu
        #check
        for i in range(pk.Theta):
          su[i] = pk.s[j][i]*u[i]
        
        sumt = sum2(su)
        sumt = mod(sumt, 2**(pk.kap+1))

    return u

class Pk(object):
  def __init__(self):
    self.lam = 12
    self.rho = 26 #p
    self.rhoi = self.rho # + self.lam
    self.eta = 1988 #(n)
    self.gam = 147456 #y
    self.Theta = 150 #O
    self.theta = 15 #0
    self.n = 4 #ceil(log2(theta+1))
    self.kap = self.gam + self.eta + 2
    self.alpha = 936
    self.alphai = self.alpha#??
    self.tau = 188
    self.l = 10
    self.log = 3 #math.log2(l)

    self.p = [random_prime(2**(self.eta-1), 2**self.eta) for i in range(self.l)] #fix TODO ?????
    self.pi = reduce((lambda x, y: x * y), self.p) #product of p

    self.q0 = (2**self.gam)
    while (self.q0 > (2**self.gam)//self.pi):
      self.q0prime1 = random_prime(0, 2**(self.lam**2))
      self.q0prime2 = random_prime(0, 2**(self.lam**2))
      self.q0 = self.q0prime1*self.q0prime2
    self.x0=self.pi*self.q0
    
    self.x = PRIntegersDelta.encrypt(self,[0 for i in range(self.tau)],self.rhoi-1,0)
    self.xi = PRIntegersDelta.encrypt(self,[0 for i in range(self.l)],self.rho,1)
    self.ii = PRIntegersDelta.encrypt(self,[0 for i in range(self.l)],self.rho,2)
   
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

    self.u = PRIntegersU.encrypt(self)

    self.o = PRIntegersDelta.encrypt(self,[0 for i in range(self.Theta)],self.rho,3)


    self.olist = list(self.o)

  def encrypt(self,m): #vector in {0,1}^l
    b = [random_element(-2**self.alphai,2**self.alpha) for i in range(self.tau)]
    bi= [random_element(-2**self.alphai,2**self.alphai) for i in range(self.l)] #shud be alpha i instead

    sums=sum2([mj*xij for mj,xij in zip(m,self.xi)])+sum2([bij*iij for bij,iij in zip(bi,self.ii)])+sum2([bj*xj for bj,xj in zip(b,self.x)])

    return modNear(sums,self.x0)

  def decrypt(self,c):
    return [mod(modNear(c,self.p[i]),2) for i in range(self.l)]

  def decrypt_squashed(self,c):
    y = [Fraction(ui)/Fraction(2**self.kap) for ui in self.u]
    z = [modNear(round((Fraction(c)*yi),4),2) for yi in y]

    m = [0 for i in range(self.l)]
    for j in range(self.l):
      sums = sum2([sji*zi for sji,zi in zip(self.s[j],z)])
      sumsb = sum2([mod(sji*zi,2) for sji,zi in zip(self.s[j],z)])
      
      m[j] = mod(int(round(sums)),2) ^ mod(c,2)
    return (m)

  def noise(self,c): #maybe?? shud prolly check this shit
    return modNear(c,self.pi)

  def add(self,c1,c2):
    return mod(c1+c2,self.x0)

  def sub(self,c1,c2):
    return mod(c1-c2,self.x0)

  def mult(self,c1,c2):
    return mod(c1*c2,self.x0)

  def __repr__(self):
    return "<Pk with rho=%d, eta=%d, gam=%d>" % (self.rho,self.eta,self.gam) 

  def recrypt(self,c):
    #"expand"
    y = [Fraction(ui)/Fraction(2**self.kap) for ui in self.u]
    z = [mod(round((Fraction(c)*yi),4),2) for yi in y]#adjust bits of precision
    
    adjz = [int(round(zi*16)) for zi in z]

    #put z in binary arrays
    zbin = [toBinary(zi,self.n+1) for zi in adjz]
   
    sec = self.olist
   
    li = [arraymult(ski,cei) for ski,cei in zip(sec,zbin)]

    Q_adds = [0 for i in range(self.n+1)]
    Q_adds2 = [0 for i in range(self.n+1)]

    print(len(li))
    for t in range(self.Theta):
      Q_adds = sumBinary(Q_adds,li[t])
      Q_adds = [mod(qa,self.x0) for qa in Q_adds]

    rounded = Q_adds[-1] + Q_adds[-2] #"round"
 
    final = rounded + (c & 1)

    return final


#fuck around
print("check")

pk = Pk()

#ii = list(pk.ii.__iter__())
#xi = list(pk.xi.__iter__())
#x = list(pk.x.__iter__())

#print(pk.eta >= (pk.alphai+pk.rhoi+1+pk.log))
#print(pk.eta >= (pk.lam*(math.log(pk.lam)**2))*pk.rho)

#for i in ii:
#  if (i//pk.pi > pk.q0):
#    print("alert")

#for i in xi:
#  if (i//pk.pi > pk.q0):
#    print("alert")

#for i in x:
#  if (i//pk.pi > pk.q0):
#    print("alert")

#print("all clear")



for i in range(1):
  print(i)
  #c1 = Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,1])
  #c0 = Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,0])
  c11 = Ciphertext.encrypt(pk, [1,1,1,1,1,1,1,1,1,1])
  c00 = Ciphertext.encrypt(pk, [0,1,0,1,0,1,0,1,0,1])
  #print(c00.noise())

  cd = [0,1,0,1,0,1,0,1,0,1]
  c00 = c11*c00
  cd = c00.decrypt()
  print(cd)
  
  r00 = c00.recrypt()

  print("recrypted")
  print(r00.decrypt())

  n00 = r00*c11
  print(n00.decrypt())





  #ca = c1 + c0
  #ca0 = c0 + c00
  #ca1 = c11 + c00
  #cm = c1 * c0
  #cm0 = c0 * c00
  #cm1 = c1 * c11
  #print(c11.decrypt_squashed()) #0011
  #print(c00.decrypt_squashed()) #01010101
  #print(c1.decrypt_squashed()) # 00001
  #print(c0.decrypt_squashed()) # 00000
  #print(ca0.decrypt()) # 0
  #print(ca1.decrypt()) # 10101010
  #print(cm.decrypt()) # 0
  #print(cm0.decrypt()) # 0
  #print(cm1.decrypt()) # 1
  #c1r = c1.recrypt()
  #c0r = c0.recrypt()
  #c11r = c11.recrypt()
  #c00r = c00.recrypt()
  #print("recrypting")
  #print(ca.recrypt()) #1111111
  #print(cb.recrypt()) #01010101
  #print(c11r.decrypt_squashed()) #111111
  #print(c00r.decrypt_squashed()) #01010101
  #print(c1r.decrypt_squashed()) # 00001
  #print(c0r.decrypt_squashed()) # 00000
  #print(c11r.decrypt()) #1111111
  #print(c00r.decrypt()) #01010101
  #print(c1r.decrypt()) # 00001
  #print(c0r.decrypt()) # 00000

  
