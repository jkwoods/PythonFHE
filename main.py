import scheme


#fuck around
print("check")
pk = scheme.Pk()
print("pk made")


c1 = scheme.Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,1])
c0 = scheme.Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,0])

print(c1.decrypt())
  
print((c1.recrypt()).decrypt())

c11 = scheme.Ciphertext.encrypt(pk, [1,1,1,1,1,1,1,1,1,1])
c00 = scheme.Ciphertext.encrypt(pk, [0,1,0,1,0,1,0,1,0,1])

ca = c11*c00
print(ca.decrypt())

cb = c11+c00
print(cb.decrypt())
  
print((ca.recrypt()).decrypt())

print((c1*c1).decrypt())

