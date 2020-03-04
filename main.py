import scheme
from timeit import default_timer as timer

#test

pk = scheme.Pk(2)

start = timer()
pk = scheme.Pk(2)
end = timer()
print("Key Time = ") 
print(end - start)

c1 = scheme.Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,1])
start = timer()
c0 = scheme.Ciphertext.encrypt(pk, [0,1,0,1,0,1,0,1,0,1])
end = timer()
print("Encode Time = ")
print(end - start)

c1.decrypt()
start = timer()
c0.decrypt()
end = timer()
print("Decode Time = ")
print(end - start)


print((c1.recrypt()).decrypt())

c11 = scheme.Ciphertext.encrypt(pk, [1,1,1,1,1,1,1,1,1,1])
c00 = scheme.Ciphertext.encrypt(pk, [0,1,0,1,0,1,0,1,0,1])

ca = c11*c00
ca.recrypt()
start = timer()
ca.recrypt()
end = timer()
print("Recode Time = ")
print(end - start)

cb = c11+c00
print(cb.decrypt())
  
print((ca.recrypt()).decrypt())

print((c1*c1).decrypt())

