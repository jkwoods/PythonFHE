import scheme


#fuck around
print("check")
pk = scheme.Pk()
print("pk made")


c1 = scheme.Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,1])
c0 = scheme.Ciphertext.encrypt(pk, [0,0,0,0,0,0,0,0,0,0])

print(c1.decrypt())
  
print((c1.recrypt()).decrypt())






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
