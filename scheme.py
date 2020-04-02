from dask import delayed
from dask import bag
import dask
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

def CRTsub(ai,ni,prod):
  p = (prod // ni)
  return ai * mul_inv(p, ni) * p


def CRT(n, a, pi): #chinese remiander thm; both inputs = bag
    sum = 0
    #prod = n.fold(lambda x, y: x*y)

    to_sum = (bag.zip(a,n)).starmap(lambda ai,ni: CRTsub(ai,ni,pi))
    summed = to_sum.fold(lambda x,y: x+y)

    return bag.compute(summed)[0] % pi

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

def make_pri(x0,ell,seed): #generates X, CANNOT DELAY, order matters
    set_random_seed(seed)
    chi = [random_element(0,x0) for i in range(ell)]
    set_random_seed(0)
    return chi

def make_deltas(pk,lenv,rho,seed,cr,partitions):
  print("enter make_deltas; partition arg =")
  print(partitions)

  pr        = bag.from_sequence(make_pri(pk.x0,lenv,seed), npartitions=partitions)
  print("pr partitions=")
  print(pr.npartitions)

  r         = [[random_element(-2**rho+1,2**rho) for i in range(pk.l)] for j in range(lenv)]

  E         = bag.from_sequence([random_element(0,(2**(pk.lam+pk.log+(pk.l*pk.eta)))//pk.pi) for i in range(lenv)], npartitions=partitions) #added from paper


  if (cr == 0):#x
  
    crts = bag.from_sequence([CRT(pk.p,bag.from_sequence([ri for ri in r[j]], npartitions=partitions).map(lambda ri: 2*ri),pk.pi) for j in range(lenv)], npartitions=partitions)
  elif (cr == 1):#xi
    crts = bag.from_sequence([CRT(pk.p,bag.from_sequence([2*ri+kd(i,j) for ri,i in zip(r[j],range(pk.l))], npartitions=partitions),pk.pi) for j in range(lenv)], npartitions=partitions)
  elif (cr == 2):#ii
    crts= bag.from_sequence([CRT(pk.p,bag.from_sequence([2*ri+(kd(i,j)*(2**(pk.rhoi+1))) for ri,i in zip(r[j],range(pk.l))], npartitions=partitions),pk.pi) for j in range(lenv)], npartitions=partitions)
  else: #o
    crts = bag.from_sequence([CRT(pk.p,bag.from_sequence([2*ri+si for ri,si in zip(r[j],pk.verts[j])], npartitions=partitions),pk.pi) for j in range(lenv)], npartitions=partitions)



  temp= pr.map(lambda xi: xi % pk.pi)
  delta = (bag.zip(temp,E,crts)).starmap(lambda te,ei,crti: te+(ei*pk.pi)-crti)
  
  print("crts partitions=")
  print(crts.npartitions)
  print("delta partitions=")
  print(delta.npartitions)


  print("leave make_deltas")
  return delta

def make_u_front(pk,seed):
  u = make_pri(2**(pk.kap+1),pk.Theta,seed) #not bag

  n=0
  x_p = pk.p.map(lambda i: (2**pk.kap)//i).compute()

  for j in range(pk.l):
    xpj = x_p[j]

    su = [delayed(lambda x,y: x*y)(pk.s[j][i],u[i]) for i in range(pk.Theta)]
    
    v = n
    n = n+1

    #change corresponding u
    su[v] = 0
    sumv = sum(su).compute()
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
  def __init__(self, key_size, partitions):
    self.lam = 42
    self.rho = 26
    self.eta = 988
    self.gam = 290000
    self.Theta = 150
    self.alpha = 936
    self.tau = 188
    self.l = 10

    self.partitions = partitions

    if (key_size==-1):
      print("correctness key test")
      self.lam = 12
      self.rho = 26 #p
      self.eta = 1988 #(n)
      self.gam = 147456 #y
      self.Theta = 150 #O
      self.alpha = 936
      self.tau = 188
      self.l = 10
    elif (key_size==0):
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

    self.rhoi = self.rho # + self.lam
    self.alphai = self.alpha#??
    
    self.p = bag.from_sequence([random_prime(2**(self.eta-1), 2**self.eta) for i in range(self.l)], npartitions=self.partitions)
    self.pi = self.p.fold(lambda x, y: x * y).compute()

    self.q0 = (2**self.gam)
    while (self.q0 > (2**self.gam)//self.pi):
      q0prime1 = delayed(random_prime)(0, 2**(self.lam**2))#delayed
      q0prime2 = delayed(random_prime)(0, 2**(self.lam**2))#delayed
      q0_temp = q0prime1*q0prime2
      self.q0 = q0_temp.compute()
    self.x0=self.pi*self.q0
    
    self.x_seed = random.randint(2, 2**30)
    self.xi_seed = random.randint(2, 2**30)
    self.ii_seed = random.randint(2, 2**30)

    self.x_deltas = make_deltas(self,self.tau,self.rhoi-1,self.x_seed,0,self.partitions) #returns bag
    self.xi_deltas = make_deltas(self,self.l,self.rho,self.xi_seed,1,self.partitions)
    self.ii_deltas = make_deltas(self,self.l,self.rho,self.ii_seed,2,self.partitions)

    print("xi_deltas partitions=")
    print(self.xi_deltas.npartitions)

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

    self.u_front = bag.from_sequence(make_u_front(self, self.u_seed), npartitions=self.partitions)

    self.o_deltas = make_deltas(self,self.Theta,self.rho,self.o_seed,3,self.partitions)
   
  def encrypt(self,m_array): #vector in {0,1}^l
    b   = bag.from_sequence([random_element(-2**self.alpha,2**self.alpha) for i in range(self.tau)], npartitions=self.partitions)
    bi  = bag.from_sequence([random_element(-2**self.alphai,2**self.alphai) for i in range(self.l)], npartitions=self.partitions)

    x   = bag.zip(bag.from_sequence(make_pri(self.x0,self.tau,self.x_seed), npartitions=self.partitions),self.x_deltas).starmap(lambda c,d: c-d)
    xi  = bag.zip(bag.from_sequence(make_pri(self.x0,self.l,self.xi_seed), npartitions=self.partitions),self.xi_deltas).starmap(lambda c,d: c-d)
    ii  = bag.zip(bag.from_sequence(make_pri(self.x0,self.l,self.ii_seed), npartitions=self.partitions),self.ii_deltas).starmap(lambda c,d: c-d)

    m = bag.from_sequence(m_array, npartitions=self.partitions)

    print("m partitions=")
    print(m.npartitions)
    print("xi partitions=")
    print(xi.npartitions)

    m_xi = (bag.zip(m,xi)).starmap(lambda x,y: x*y)
    bi_ii = (bag.zip(bi,ii)).starmap(lambda x,y: x*y)
    b_x = (bag.zip(b,x)).starmap(lambda x,y: x*y)

    sums= (bag.concat([m_xi,bi_ii,b_x])).fold(lambda x,y: x+y)
    return modNear(sums.compute(),self.x0)

  def decrypt(self,c):
    pbag = bag.from_sequence(self.p, npartitions=self.partitions)
    m = pbag.map(lambda pi: mod(modNear(c,pi),2))
    return m.compute()
    

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
    u_end = bag.from_sequence(u_draft[self.l:])
    u = bag.concat([self.u_front,u_end])
    u = u.repartition(self.partitions)
    print("u part")
    print(u.npartitions)


    #"expand" #TODO CHANGE FRACTION TO RATIONAL
    y = u.map(lambda ui: Fraction(ui)/Fraction(2**self.kap))
    print("y part")
    print(y.npartitions)

    z = y.map(lambda yi: mod(round((Fraction(c)*yi),4),2))
    print("z part")
    print(z.npartitions)

    adjz = z.map(lambda zi: int(round(zi*16)))
    print("adjz part")
    print(adjz.npartitions)

    #put z in binary arrays
    zbin = adjz.map(lambda zi: toBinary(zi,self.n+1))
    print("zbin part")
    print(zbin.npartitions)

    #get o
    o = bag.zip(bag.from_sequence(make_pri(self.x0,self.Theta,self.o_seed), npartitions=self.partitions),self.o_deltas).starmap(lambda c,d: c-d)

    o_z = (bag.zip(o,zbin)).starmap(lambda oi,zi: arraymult(oi,zi))
   
    Q_adds = o_z.fold(lambda q,r: sumBinary(q,r))

    Q_sum = Q_adds.compute()

    rounded = Q_sum[-1] + Q_sum[-2] #"round"
 
    final = rounded + (c & 1)

    '''
    li = bag.compute(o_z)[0]
    Q_adds = [0 for i in range(self.n+1)]

    for t in range(self.Theta):
      Q_adds = sumBinary(Q_adds,li[t])
      #Q_adds = [mod(qa,self.x0) for qa in Q_adds]

    rounded = Q_adds[-1] + Q_adds[-2] #"round"
 
    final = rounded + (c & 1)
    '''

    return final
