from dask import delayed
import dask.array as da
import numpy as np
import dask
import operator
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

@dask.delayed
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

@dask.delayed
def frac1(u,k):
  return Fraction(u)/Fraction(2**k)

@dask.delayed
def frac2(c,y):
  return mod(round((Fraction(c)*y),4),2)

@dask.delayed
def round3(zi):
  return int(round(zi*16))

@dask.delayed
def toBinary(x,l):
  "Converts a positive integer x into binary with l digits"
  if (x==32): return np.array([0]*l)

  return np.array(digits(x+2**l)[:-1])

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

def CRT(n, a, prod, prod_div_array): #chinese remiander thm   #n = 1 d - array, a = 2 d rand array 2b split, should return 1-day array

  #mul_array = n.map_blocks((lambda ni: mul_inv((prod // ni), ni) * (prod // ni)), dtype=object) #1d
  mul_array0 = []
  for i in range(len(prod_div_array)):
    mul_array0.append(mul_inv(prod_div_array[i], n[i]) * prod_div_array[i])

  mul_array = da.from_array(mul_array0, chunks=1)
  pre = mul_array.reshape(mul_array.shape[0],1)
  pre2 = da.transpose(pre)

  #a_p = da.multiply(pre2,a) #2d
  a_p = da.blockwise(operator.mul, 'ij', pre2, 'ij', a, 'ij', dtype=object)

  summed = da.sum(a_p, axis=1) #1d #TODO axis correct??

  return summed.map_blocks((lambda s: s % prod), dtype=object)

def kd(i,j):
  if (i == j):
    return 1
  else:
    return 0

@dask.delayed
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

def make_deltas(pk,lenv,rho,seed,cr):
  pr=da.from_array(make_pri(pk.x0,lenv,seed), chunks=1)
 
  e_help = (2**(pk.lam+pk.log+(pk.l*pk.eta)))//pk.pi
  
  r0=[[delayed(random_element)(-2**rho+1,2**rho) for i in range(pk.l)] for j in range(lenv)]
  E0=[delayed(random_element)(0,e_help) for i in range(lenv)] #added from paper
  
  r1 = dask.compute(*r0)
  E1 = dask.compute(*E0)

  r = da.from_array(r1, chunks=1)
  E = da.from_array(E1, chunks=1)
  
  kd_array = da.from_array([[kd(i,j) for i in range(pk.l)] for j in range(lenv)])

  #delta=[0 for i in range(lenv)]

  if (cr == 0):#x
    #crts = [CRT(pk.p,[2*ri for ri in r[j]]) for j in range(lenv)]
    r2 = r.map_blocks((lambda x: 2*x), dtype=object)
    crts = CRT(pk.p_unwrapped,r2,pk.pi,pk.pdiv)
  elif (cr == 1):#xi
    #crts = [CRT(pk.p,[2*ri+kd(i,j) for ri,i in zip(r[j],range(pk.l))]) for j in range(lenv)]
    r2 = r.map_blocks((lambda x: 2*x), dtype=object)
    
    #r3 = da.add(r2,kd_array)
    r3 = da.blockwise(operator.add, 'ij', r2, 'ij', kd_array, 'ij', dtype=object)
    
    crts = CRT(pk.p_unwrapped,r3,pk.pi,pk.pdiv)
  elif (cr == 2):#ii
    #crts = [CRT(pk.p,[2*ri+(kd(i,j)*(2**(pk.rhoi+1))) for ri,i in zip(r[j],range(pk.l))]) for j in range(lenv)]
    r2 = r.map_blocks((lambda x: 2*x), dtype=object)
    kd2 = kd_array.map_blocks((lambda x: (2**(pk.rhoi+1)*x)), dtype=object)
    
    #r3 = da.add(r2,kd2)
    r3 = da.blockwise(operator.add, 'ij', r2, 'ij', kd2, 'ij', dtype=object)

    crts = CRT(pk.p_unwrapped,r3,pk.pi,pk.pdiv)
  else: #o
    #crts = [CRT(pk.p,[2*ri+si for ri,si in zip(r[j],pk.verts[j])]) for j in range(lenv)]
    r2 = r.map_blocks((lambda x: 2*x), dtype=object)
    
    #r3 = da.add(r2,pk.verts)
    r3 = da.blockwise(operator.add, 'ij', r2, 'ij', pk.verts, 'ij', dtype=object)

    crts = CRT(pk.p_unwrapped,r3,pk.pi,pk.pdiv)

  #temp=[mod(Xi,pk.pi) for Xi in pr]
  #delta=[te+(ei*pk.pi)-crti for te,ei,crti in zip(temp,E,crts)] #changed from paper

  
  temp = pr.map_blocks((lambda r: r % pk.pi), dtype=object)
  delta1 = E.map_blocks((lambda ei: pk.pi*ei), dtype=object)

  #delta2 = da.add(temp,delta1)
  #delta3 = da.subtract(delta2,crts)

  delta2 = da.blockwise(operator.add, 'i', temp, 'i', delta1, 'i', dtype=object) #TODO check?
  delta3 = da.blockwise(operator.sub, 'i', delta2, 'i', crts, 'i', dtype=object) #TODO check?

  #delta.visualize(filename=(str(cr)+'deltapart.svg'))

  return delta3

def make_u_front(pk,seed):
  pr=make_pri(2**(pk.kap+1),pk.Theta,seed) #u draft
  u = pr

  n=0
  for j in range(pk.l):
    xpj = (2**pk.kap)//pk.p_unwrapped[j]

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

    #self.p_unwrapped = read_primes(key_size,0)
    #lam_prime_list = read_primes(key_size,1)

    self.p_unwrapped = [random_prime(2**(self.eta-1), 2**self.eta) for i in range(self.l)]
    self.p = da.from_array(self.p_unwrapped, chunks=1).persist() #fix TODO ?????
    #self.pi = reduce((lambda x, y: x * y), self.p) #product of p
    pi_wrapped = da.prod(self.p)
    self.pi = pi_wrapped.compute()

    prod_array0 = da.broadcast_to(pi_wrapped, (self.l,))#expand pi
    prod_array1 = da.blockwise(operator.floordiv, 'i', prod_array0, 'i', self.p, 'i', dtype=object) #1d
    self.pdiv = prod_array1.compute()

    self.q0 = (2**self.gam)
    i = 0
    while (self.q0 > (2**self.gam)//self.pi):
      q0prime1 = delayed(random_prime)(0, 2**(self.lam**2))
      q0prime2 = delayed(random_prime)(0, 2**(self.lam**2))
      #q0prime1 = lam_prime_list[i]
      #q0prime2 = lam_prime_list[i+1]
      i = i + 2
      self.q0 = (q0prime1*q0prime2).compute()

    self.x0=self.pi*self.q0

    self.x_seed = random.randint(2, 2**30)
    self.xi_seed = random.randint(2, 2**30)
    self.ii_seed = random.randint(2, 2**30)

    #TODO MOVE pre CRT calcualtions

    self.x_deltas = make_deltas(self,self.tau,self.rhoi-1,self.x_seed,0).persist() #dask array 
    self.xi_deltas = make_deltas(self,self.l,self.rho,self.xi_seed,1).persist()
    self.ii_deltas = make_deltas(self,self.l,self.rho,self.ii_seed,2).persist()

    self.B=self.Theta//self.theta

    self.s = [[0 for j in range(self.Theta)] for k in range(self.l)]

    for j in range(self.l): #TODO DASK
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

    self.verts = [[0 for j in range(self.l)] for k in range(self.Theta)] #TODO DASK TRANSPOSE

    for i in range(self.Theta):
      for j in range(self.l):
        self.verts[i][j] = self.s[j][i]

    self.verts = da.from_array(self.verts)

    self.u_seed = random.randint(2, 2**30)
    self.o_seed = random.randint(2, 2**30)

    self.u_front = make_u_front(self, self.u_seed) #make future TODO

    self.o_deltas = make_deltas(self,self.Theta,self.rho,self.o_seed,3).persist()
   
  def encrypt(self,m): #vector in {0,1}^l
    b0 = [delayed(random_element)(-2**self.alpha,2**self.alpha) for i in range(self.tau)]
    bi0= [delayed(random_element)(-2**self.alphai,2**self.alphai) for i in range(self.l)]

    b1 = dask.compute(*b0)
    bi1 = dask.compute(*bi0)

    b = da.from_array(b1, chunks=1)
    bi = da.from_array(bi1, chunks=1)

    #x = [c-d for c,d in zip(make_pri(self.x0,self.tau,self.x_seed),self.x_deltas)]
    #xi= [c-d for c,d in zip(make_pri(self.x0,self.l,self.xi_seed),self.xi_deltas)]
    #ii= [c-d for c,d in zip(make_pri(self.x0,self.l,self.ii_seed),self.ii_deltas)]

    #x = da.subtract((da.from_array(make_pri(self.x0,self.tau,self.x_seed), chunks=1)),self.x_deltas)
    #xi = da.subtract((da.from_array(make_pri(self.x0,self.l,self.xi_seed), chunks=1)),self.xi_deltas)
    #ii = da.subtract((da.from_array(make_pri(self.x0,self.l,self.ii_seed), chunks=1)),self.ii_deltas)

    x = da.blockwise(operator.sub, 'i', (da.from_array(make_pri(self.x0,self.tau,self.x_seed), chunks=1)), 'i', self.x_deltas, 'i', dtype=object) #TODO check?
    xi = da.blockwise(operator.sub, 'i', (da.from_array(make_pri(self.x0,self.l,self.xi_seed), chunks=1)), 'i', self.xi_deltas, 'i', dtype=object) #TODO check?
    ii = da.blockwise(operator.sub, 'i', (da.from_array(make_pri(self.x0,self.l,self.ii_seed), chunks=1)), 'i', self.ii_deltas, 'i', dtype=object) #TODO check?

    #sums=sum([mj*xij for mj,xij in zip(m,xi)])+sum([bij*iij for bij,iij in zip(bi,ii)])+sum([bj*xj for bj,xj in zip(b,x)])
    #rmn = modNear(sums,self.x0)
    #rmn.visualize(filename='encrypt.svg')

    #m_xi = da.sum(da.multiply(da.from_array(m, chunks=1),xi))
    #bi_ii = da.sum(da.multiply(bi,ii))
    #b_x = da.sum(da.multiply(b,x))

    m_xi = da.blockwise(operator.mul, 'i', (da.from_array(m, chunks=1)), 'i', xi, 'i', dtype=object) #TODO check?
    bi_ii = da.blockwise(operator.mul, 'i', bi, 'i', ii, 'i', dtype=object) #TODO check?
    b_x = da.blockwise(operator.mul, 'i', b, 'i', x, 'i', dtype=object) #TODO check?

    #big_sum = da.sum(m_xi)+da.sum(bi_ii)+da.sum(b_x) DOESN"T RUN
    #cut x into chunks, add everything up
    half = da.blockwise(operator.add, 'i', m_xi, 'i', bi_ii, 'i', dtype=object) #TODO check?
    r = self.tau//self.l
    #for i in range(r+1):
    #  half = da.blockwise(operator.add, 'i', half, 'i', b_x[(i*self.l):((i*self.l)+self.l)], 'i', dtype=object) #TODO check?

    big = sum(half.compute()) + sum(b_x.compute())

    final = modNear(big,self.x0)

    return final

  def decrypt(self,c):
    #c-pi*((2*c+pi)//(2*pi))
    p_2 = self.p.map_blocks((lambda pi: pi * 2), dtype=object)
    c_2_plus_p = self.p.map_blocks((lambda pi: (2*c)+pi), dtype=object)
    qnear = da.blockwise(operator.floordiv, 'i', c_2_plus_p, 'i', p_2, 'i', dtype=object)
    qnear_p = da.blockwise(operator.mul, 'i', qnear, 'i', self.p, 'i', dtype=object)
    modresult = qnear_p.map_blocks((lambda pi: c-pi), dtype=object)
    #modresult % 2 #TODO change to &?
    result = modresult.map_blocks((lambda mi: mi % 2), dtype=object)
    return result.compute()

    
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
    y = [frac1(ui,self.kap) for ui in u]
    z = [frac2(c,yi) for yi in y]#adjust bits of precision
    adjz = [round3(zi) for zi in z]
    zbin = [toBinary(zi,self.n+1) for zi in adjz]
    z_comp = np.array(dask.compute(*zbin))
    print("zbin computed")
    #get o

    o = da.blockwise(operator.sub, 'i', (da.from_array(make_pri(self.x0,self.Theta,self.o_seed), chunks=1)), 'i', self.o_deltas, 'i', dtype=object) #TODO check?
    o1 = o.reshape(o.shape[0],1)
    #new_o = da.concatenate((o1,o1,o1,o1,o1), axis=1)
    #zarr = da.from_array(zcom, chunks=1)

    o_comp = dask.compute(*o1)

    #li = da.blockwise(operator.mul, 'ij', new_o, 'ij', zarr, 'ij', dtype=object)

    li = [arraymult(ski,cei) for ski,cei in zip(o_comp,z_comp)]

    Q_adds = [0 for i in range(self.n+1)]

    for t in range(self.Theta):
      Q_adds = sumBinary(Q_adds,li[t])

    rounded = Q_adds[-1] + Q_adds[-2] #"round"
 
    final = rounded + (c & 1)
    final.visualize(filename='2halfRecrypt.svg')

    return final.compute()
