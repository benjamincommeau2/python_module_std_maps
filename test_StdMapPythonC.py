import numpy as np
import StdMapPythonC
x=np.zeros((4,3),dtype=np.cdouble)
x+=np.random.randn(4,3)
x+=1j*np.random.randn(4,3)
print("x=")
print(x)
print("StdMapPythonC.fgt(x)=")
print(StdMapPythonC.fgt(x))
