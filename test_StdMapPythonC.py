import numpy as np
import StdMapPythonC
x=np.zeros((4,4),dtype=np.cdouble)
x+=np.random.randn(4,4)
x+=1j*np.random.randn(4,4)
print("np.dtype(x[0][0])=")
print(np.dtype(x[0][0]))
print("x=")
print(x)
print("StdMapPythonC.fgt(x)=")
print(StdMapPythonC.fgt(x))
