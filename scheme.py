import numpy as np
import operator
import random
import sympy
from functools import reduce 
import math
from fractions import Fraction
from itertools import islice
import gmpy2
from gmpy2 import mpz,mpq,c_div,random_state,mpz_random

class Rand(object):
  def __init__(self):
    self.state = random_state(42)

  def random_element(self,lb,ub,rs=None):
    if rs==None: rs=self.state
    r = mpz_random(rs,ub-lb)
    return r+lb

  def make_pri(self, x0,ell,seed):
    rsn = random_state(seed)
    chi = [self.random_element(0,x0,rsn) for i in range(ell)]
    return chi

def QuotientNear(a,b):
  "Gives the nearest integer to a/b"
  return (2*a+b)//(2*b)

def modNear(a,b):
  "Computes a mod b with a \in ]-b/2,b/2]"
  return a-b*QuotientNear(a,b)

def mod(c, p):
  return c % p

def random_prime(a,b):
  return mpz(sympy.randprime(a,b))

def sumBinary(a,b):
  c=[a[0]+b[0]]
  carry=a[0]*b[0]

  for i in range(1,len(a)-1):
    carry2=(a[i]+b[i])*carry+a[i]*b[i]    
    c.append(a[i]+b[i]+carry)             
    carry=carry2
  
  c.append(a[-1]+b[-1]+carry)
  return c

def arraymult(c,a):
  return [c*int(xi) for xi in a]

def frac1(u,k,c):
  return mpq(u,2**k)*c

def frac2(y):
  return mpz((y%2)*32)

def round3(z):
  return c_div(z,mpz(2))

def toBinary(x,l):
  if (x==32): return np.array([0]*l)
  return np.array(digits(x+2**l)[:-1])

def digits(x):
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

def CRT(n, a, prod, prod_div_array): #chinese remiander thm   #n = 1 d - array, a = 2 d rand array 2b split, should return 1-day array

  #mul_array = n.map_blocks((lambda ni: mul_inv((prod // ni), ni) * (prod // ni)), dtype=object) #1d
  mul_array = []
  for i in range(len(prod_div_array)):
    mul_array.append(mul_inv(prod_div_array[i], n[i]) * prod_div_array[i])

  a_p = [[a[i][j]*mul_array[j] for j in range(len(a[0]))] for i in range(len(a))] #TODO CHECK CORRECTNESS

  summed = [sum(a_p[i]) for i in range(len(a))] #1d #TODO axis correct??

  return [(s % prod) for s in summed]


def make_deltas(pk,lenv,rho,seed,cr):
  pr=pk.rgen.make_pri(pk.x0,lenv,seed)
  e_help = 2**(pk.lam+pk.log+(pk.l*pk.eta))//pk.pi

  r=[[pk.rgen.random_element(-2**rho+1,2**rho) for i in range(pk.l)] for j in range(lenv)]
  E=[pk.rgen.random_element(0,e_help) for i in range(lenv)] #added from paper

  kd_array = [[kd(i,j) for i in range(pk.l)] for j in range(lenv)]
  r2 = [[2*r[j][i] for i in range(pk.l)] for j in range(lenv)]
  #delta=[0 for i in range(lenv)]

  if (cr == 0):#x
    #crts = [CRT(pk.p,[2*ri for ri in r[j]]) for j in range(lenv)]
    crts = CRT(pk.p,r2,pk.pi,pk.pdiv)

  elif (cr == 1):#xi
    #crts = [CRT(pk.p,[2*ri+kd(i,j) for ri,i in zip(r[j],range(pk.l))]) for j in range(lenv)]
    r3 = [[kd_array[j][i]+r2[j][i] for i in range(pk.l)] for j in range(lenv)]
    
    crts = CRT(pk.p,r3,pk.pi,pk.pdiv)
  elif (cr == 2):#ii
    print("cr2")
    #crts = [CRT(pk.p,[2*ri+(kd(i,j)*(2**(pk.rhoi+1))) for ri,i in zip(r[j],range(pk.l))]) for j in range(lenv)]
    r3 = [[(kd_array[j][i]*(2**(pk.rhoi+1)))+r2[j][i] for i in range(pk.l)] for j in range(lenv)]

    crts = CRT(pk.p,r3,pk.pi,pk.pdiv)
  else: #o
    print("cr3")
    #crts = [CRT(pk.p,[2*ri+si for ri,si in zip(r[j],pk.verts[j])]) for j in range(lenv)]
    r3 = [[pk.verts[j][i]+r2[j][i] for i in range(pk.l)] for j in range(lenv)] #TODO check correctness

    crts = CRT(pk.p,r3,pk.pi,pk.pdiv)

  #temp=[mod(Xi,pk.pi) for Xi in pr]
  #delta=[te+(ei*pk.pi)-crti for te,ei,crti in zip(temp,E,crts)] #changed from paper
  
  temp = [(r % pk.pi) for r in pr]
 
  delta=[te+(ei*pk.pi)-crti for te,ei,crti in zip(temp,E,crts)]

  return delta


def make_u_front(pk,seed):
  pr=pk.rgen.make_pri(2**(pk.kap+1),pk.Theta,seed) #u draft
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


def kd(i,j):
  if (i == j):
    return 1
  else:
    return 0



class Ciphertext(object):
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


def read_primes(key,el):
  plist = []
  if (el==0): #eta
    s_type = "eta"
    
  elif (el==1):  #lam
    s_type = "lam"

  with open((s_type + str(key)), "r") as prime_file:
    lines = prime_file.readlines()
    for p in lines:
      plist.append(int(p))

  return plist


class Pk(object):
  def __init__(self, key_size):
    self.lam = 42
    self.rho = 26
    self.eta = 988
    self.gam = 290000
    self.Theta = 150
    self.alpha = 936
    self.tau = 188
    self.l = 10

    if (key_size==-1):
      print("correctness key test")
      self.lam = mpz(12)
      self.rho = mpz(26) #p
      self.eta = mpz(1988) #(n)
      self.gam = mpz(147456) #y
      self.Theta = 150 #O
      self.alpha = mpz(936)
      self.tau = 188
      self.l = 10
    elif (key_size==0):
      print("Making toy key")
      self.lam = mpz(42)
      self.rho = mpz(26)
      self.eta = mpz(988)
      self.gam = mpz(290000)
      self.Theta = 150
      self.alpha = mpz(936)
      self.tau = 188
      self.l = 10
    elif(key_size==1):
      print("making small key")
      self.lam = mpz(52)
      self.rho = mpz(41)
      self.eta = mpz(1558)
      self.gam = mpz(1600000)
      self.Theta = 555
      self.alpha = mpz(1476)
      self.tau = 661
      self.l = 37
    elif (key_size==2):
      print("making medium key")
      self.lam = mpz(62)
      self.rho = mpz(56)
      self.eta = mpz(2128)
      self.gam = mpz(8500000)
      self.Theta = 2070
      self.alpha = mpz(2016)
      self.tau = 2410
      self. l = 138
    elif (size == 3):
      self.lam = mpz(72)
      self.rho = mpz(71)
      self.eta = mpz(2698)
      self.gam = mpz(39000000)
      self.Theta = 7965
      self.alpha = mpz(2556)
      self.tau = 8713
      self.l = 531

 
    self.alphai = self.lam + self.alpha
    self.rhoi = self.lam + self.alpha
    self.n = 4;
    self.kap = 64*(self.gam//64+1)-1
    self.log = round(math.log2(self.l))
    self.theta = self.Theta//self.l

    self.rhoi = self.rho # + self.lam
    self.alphai = self.alpha#??
    self.rgen = Rand()

    print("rs created")
    #self.p_unwrapped = read_primes(key_size,0)
    #lam_prime_list = read_primes(key_size,1)

    self.p = [random_prime(2**(self.eta-1), 2**self.eta) for i in range(self.l)]
    self.pi = reduce(operator.mul, self.p)

    self.pdiv = [self.pi // p for p in self.p]

    self.q0 = (2**self.gam)

    i = 0
    while (self.q0 > (2**self.gam)//self.pi):
      q0prime1 = random_prime(0, 2**(self.lam**2))
      q0prime2 = random_prime(0, 2**(self.lam**2))
      #q0prime1 = lam_prime_list[i]
      #q0prime2 = lam_prime_list[i+1]
      i = i + 2
      self.q0 = q0prime1*q0prime2

    self.x0=self.pi*self.q0

    self.x_seed = random.randint(2, 2**30)
    self.xi_seed = random.randint(2, 2**30)
    self.ii_seed = random.randint(2, 2**30)

    #TODO MOVE pre CRT calcualtions

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

    self.u_front = make_u_front(self, self.u_seed) #make future TODO

    self.o_deltas = make_deltas(self,self.Theta,self.rho,self.o_seed,3)

  def encrypt(self,m): #vector in {0,1}^l
    b = [self.rgen.random_element(-2**self.alpha,2**self.alpha) for i in range(self.tau)]
    bi= [self.rgen.random_element(-2**self.alphai,2**self.alphai) for i in range(self.l)]

    x = [c-d for c,d in zip(self.rgen.make_pri(self.x0,self.tau,self.x_seed),self.x_deltas)]
    xi= [c-d for c,d in zip(self.rgen.make_pri(self.x0,self.l,self.xi_seed),self.xi_deltas)]
    ii= [c-d for c,d in zip(self.rgen.make_pri(self.x0,self.l,self.ii_seed),self.ii_deltas)]
   
    m_xi = [mj*xij for mj,xij in zip(m,xi)]
    bi_ii = [bij*iij for bij,iij in zip(bi,ii)]
    b_x = [bj*xj for bj,xj in zip(b,x)]

    big = sum(m_xi) + sum(bi_ii) + sum(b_x)

    final = modNear(big,self.x0)

    return final

  def decrypt(self,c):
    return [int(mod(modNear(c,self.p[i]),2)) for i in range(self.l)]

  def add(self,c1,c2):
    return mod(c1+c2,self.x0)

  def sub(self,c1,c2):
    return mod(c1-c2,self.x0)

  def mult(self,c1,c2):
    return mod(c1*c2,self.x0)

  def recrypt(self,c):
    #get u
    u_draft = self.rgen.make_pri(2**(self.kap+1),self.Theta,self.u_seed)
    u = self.u_front+u_draft[self.l:]

    #"expand"
    y = [frac1(ui,self.kap,c) for ui in u]
    z = [frac2(yi) for yi in y]
    z1 = [round3(zi) for zi in z]
    zbin = [toBinary(zi,self.n+1) for zi in z1]

    #get o
    o = [c-d for c,d in zip(self.rgen.make_pri(self.x0,self.Theta,self.o_seed),self.o_deltas)]

    li = [arraymult(ski,cei) for ski,cei in zip(o,zbin)]

    Q_adds = [0 for i in range(self.n+1)]

    for t in range(self.Theta):
      Q_adds = sumBinary(Q_adds,li[t])

    rounded = Q_adds[-1] + Q_adds[-2] #"round"
 
    final = rounded + (c & 1)

    return final
