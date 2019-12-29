import numpy as np
import sys
import FGT
print()
list_type=[np.int,np.cdouble,np.double,np.complex]
mytype=list_type[int(sys.argv[1])]
print("("+str(sys.argv[1])+"): Start python3 dtype exception test #"+str(mytype))
x=np.zeros((4,4),dtype=mytype)
#x=1j*np.zeros((4,4))
x+=np.random.randn(4,4)
x+=1j*np.random.randn(4,4)
#print("np.dtype(x[0][0])=")
#print(np.dtype(x[0][0]))
#print("x=")
#print(x)
#print("StdMapPythonC.fgt(x)=")
try:
  y=FGT.fgt(x)
except Exception as e:
  print("("+str(sys.argv[1])+"): Failed dtype exception check using dtype="+str(mytype))
  #print(e.what())
#print(y)
print("("+str(sys.argv[1])+"): Passed dtype exception check using dtype="+str(mytype))
