import scheme
from timeit import default_timer as timer

#test

start = timer()
pk = scheme.Pk(-1)
end = timer()
print("Key Time = ") 
print(end - start)

#print((scheme.Ciphertext(1,pk)).decrypt())

c1 = scheme.Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,1])
start = timer()
c0 = scheme.Ciphertext.encrypt(pk, [0,1,0,1,0,1,0,1,0,1])
end = timer()
print("Encode Time = ")
print(end - start)

print(c1.decrypt())
start = timer()
print(c0.decrypt())
end = timer()
print("Decode Time = ")
print(end - start)


#print((c1.recrypt()).decrypt())

ca = c1*c0
ca.recrypt()
start = timer()
ca.recrypt()
end = timer()
print("Recode Time = ")
print(end - start)

cb = c1+c0
print(cb.decrypt())
  
print((ca.recrypt()).decrypt())

print((c1*c1).decrypt())
