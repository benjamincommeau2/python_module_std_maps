/*
This code was taken from Dan Foreman-Mackey and modified.
Source Code : https://gist.github.com/dfm/3247796
Source Documentation : http://dan.iel.fm/posts/python-c-extensions/
clear ; python setup_StdMapPythonC.py build_ext --inplace
clear ; python test_StdMapPythonC.py
*/
#define MOD_DEF(ob, name, doc, methods) \
    static struct PyModuleDef moduledef = { \
        PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
    ob = PyModule_Create(&moduledef);
#define NPY_NO_DEPRECATED_API NPY_1_17_API_VERSION
#include </usr/include/python3.6/Python.h>
#include </usr/local/lib/python3.6/dist-packages/numpy/core/include/numpy/arrayobject.h>
#include <stdio.h>
#include <complex.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include <cmath.h>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// c function declarations//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fht(std::complex<double> *x, int64_t n,int64_t NN, int64_t i1,int64_t j1);
void norm(std::complex<double> *x,std::complex<double> a,int64_t NN);
void g2m_coor(int64_t * x,int64_t* y,int64_t *mask);
void m2g_coor(int64_t * x,int64_t* y,int64_t *mask);
void matrix_pad_copy(std::complex<double> *x, std::complex<double> *y,int64_t *N,int64_t *NN);
int64_t CountOnesFromInteger(int64_t value);
void rear(std::complex<double> *x,int64_t NN);
void irear(std::complex<double> *x,int64_t NN);
int64_t init_arrays(PyObject *self, PyObject *args, PyObject **x_obj, PyObject **x_array,  PyObject **y_array, std::complex<double> **x,  std::complex<double> **y, int64_t *N, int64_t *NN);
void parity(std::complex<double> *z,int NN);
void iparity(std::complex<double> *z,int NN);
std::complex<double> imaginary_number=std::complex<double>{1,0};
/*
Available functions 
This is a declaration of a structure type.
This is our declaration of our 'fgt_StdMapPythonC' function for python.
It specifies all input and output objects of our function.
The name that we've given to the function (fgt_StdMapPythonC) is also a matter of convention and nessaray because {module_name}_{function_name}. 
The name for this c file given as '_StdMapPythonC' is a matter of convention and necessary.
Calling the function in python is _StdMapPythonC.fgt(), where _StdMapPythonC is the name of the module and fgt is the name of the function. 
the '_' for c is equivalent for the '.' in python in regards to class extensions
*/
static PyObject *rear_StdMapPythonC(PyObject *self, PyObject *args); 
static PyObject *irear_StdMapPythonC(PyObject *self, PyObject *args); 
static PyObject *fgt_StdMapPythonC(PyObject *self, PyObject *args);
static PyObject *ifgt_StdMapPythonC(PyObject *self, PyObject *args);
static PyObject *parity_StdMapPythonC(PyObject *self, PyObject *args);
static PyObject *iparity_StdMapPythonC(PyObject *self, PyObject *args);
static PyObject *auto_pad_StdMapPythonC(PyObject *self, PyObject *args);
static PyObject *fhgt_StdMapPythonC(PyObject *self, PyObject *args);
static PyObject *ifhgt_StdMapPythonC(PyObject *self, PyObject *args);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// python documentation ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Docstrings module_docstring
Below establishes documentation for our fgt function
In python executing 'fgt?' or 'help(fgt)' will display 'fgt_docstring' in the terminal
Malloc cannot be used, because the program has no specified shutdown time. The module is imported, not called.
 */
static char module_docstring[] =
  "Calculate the coefficients A[j_1...j_n], of a square matrix A "
  "with size N by N of complex double elements, whose weighted sum[1] "
  "with N^2 irreducible matrices ([j_1](X) ...(X)[j_n]) "
  "is equal to A. N=2^n, where 'n' is a natural number. [j] are the "
  "Pauli matrices[2] for j=(1,2,3) or identity for j=0. (X) "
  "is the Kronecker product[3]. "
  "\n"
  "\n"
  "Input: Square matrix A with size (2^n,2^n), where n is a natural number. "
  "\n"
  "\n"
  "Output: Weighted coefficients "
  "A[j_1...j_n], where j=(0,1,2,3). "
  "\n"
  "\n"  
  "Memory Management: "
  "The following functions perform memory in place calculations if no padding is performed: fgt(), ifgt(), rear(), irear(), and auto_pad()."
  "The folowing functions do not perform memory in place calculations: parity()."
  "\n"
  "\n"  
  "Details that might be relevant "
  "in identifying Missing Application: "
  "\n"
  "I wrote a program that computes the "
  "full decomposition of a N by N matrix into the sum of weighted "
  "gamma matrices[4] and their multi-linears[4] in O(N^2log(N)) "
  "steps. A is the matrix to decompose. N=2^n, where n is a natural number. "
  "A[k{j}] are the weighted elements to be computed. "
  "s^k{j}=(g^j_1)(g^j_2)(g^j_3)... "
  "are the multi-linears equal to the product of 'k' unique gamma matrices (g^j). "
  "k{j} is the jth set of uniquely chosen 'k' indices. For example in "
  "N=4 and k = 2: k{j} is pick the jth element from the list "
  "(01,02,03,04,12,13,14,23,24,34), where none of elements share the "
  "same indices. Another representation of A can be written as a sum of "
  "pauli matrices [j] expanded into each other multiple times using (X), where (X) is the Kronecker "
  "product, [1], [2], [3], are the Pauli matrices, [0] is identity, "
  "A[j_1...j_n]$ are the weighted elements to be computed in a different "
  "basis, and (j_1...j_n) can be treated as a number in base 4. One can "
  "show that multipliting two matrices ([j_1](X)...(X)[j_n]) and "
  "([l_1](X)...(X)[l_n]) to get ([g_1](X)...(X)[g_n]) "
  "is equivalent to (j_1...j_n) xor (l_1...l_n) = (g_1...g_n) "
  "but missing the sign. There exists a simple map between A[k{j}] "
  "<--> A[j_1...j_n]. "
  "\n"
  "\n"
  "References: "
  "\n"
  "[1]: https://en.wikipedia.org/wiki/Weight_function "
  "\n"
  "[2]: https://en.wikipedia.org/wiki/Pauli_matrices "
  "\n"
  "[3]: https://en.wikipedia.org/wiki/Kronecker_product "
  "\n"
  "[4]: https://en.wikipedia.org/wiki/Higher-dimensional_gamma_matrices#Symmetry_properties ";
static char rear_docstring[] =
  "Rearrange elements so each row's elements share the same irreducibles of whom who share the same non-zero elements."; 
static char irear_docstring[] =
  "Inverse of the Rearrange elements so each row's elements share the same irreducibles of whom who share the same non-zero elements."; 
static char fgt_docstring[] =
  "Computes the Fast Gamma Transform of a matrix"; 
static char ifgt_docstring[] =
  "Inverse of the Computes the Fast Gamma Transform of a matrix"; 
static char fhgt_docstring[] =
  "Computes the Fast Hermitian Gamma Transform of a matrix"; 
static char ifhgt_docstring[] =
  "Inverse of the Computes the Fast Hermitian Gamma Transform of a matrix"; 
static char parity_docstring[] =
  "Returns a coefficient which if multiplied into the gamma matrix cooresponding to the gamma cooridinate (x,y), makes it hermitian."; 
static char iparity_docstring[] =
  "Invere of the Returns a coefficient which if multiplied into the gamma matrix cooresponding to the gamma cooridinate (x,y), makes it hermitian."; 
static char auto_pad_docstring[] =
  "Returns new copy of matrix zero padded to the next power of 2 in matrix length."; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// python module initialization ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Module specification
This specifies all functional members of this python module
In our case there is only one, which is 'fgt_StdMapPythonC' as defined above
*/
static PyMethodDef module_methods[] = {
  {"rear", rear_StdMapPythonC, METH_VARARGS, rear_docstring},
  {"irear", irear_StdMapPythonC, METH_VARARGS, irear_docstring},
  {"fgt", fgt_StdMapPythonC, METH_VARARGS, fgt_docstring},
  {"ifgt", ifgt_StdMapPythonC, METH_VARARGS, ifgt_docstring},
  {"fhgt", fhgt_StdMapPythonC, METH_VARARGS, fhgt_docstring},
  {"ifhgt", ifhgt_StdMapPythonC, METH_VARARGS, ifhgt_docstring},
  {"parity", parity_StdMapPythonC, METH_VARARGS, parity_docstring},
  {"iparity", iparity_StdMapPythonC, METH_VARARGS, iparity_docstring},
  {"auto_pad", auto_pad_StdMapPythonC, METH_VARARGS, auto_pad_docstring},
  {NULL, NULL, 0, NULL} // tells the compiler there are no more methods
};
/*
Initialize the module.
This initializes are module containing our function.
This initializing module function must be named as 'init'+'ourmodulename',
which in our case will be 'init'+'_StdMapPythonC' = 'init_StdMapPythonC' .
*/
static PyObject * moduleinit(void)
{
  //PyObject *mm = PyModule_Create("_StdMapPythonC", module_methods, module_docstring);
  PyObject *mm;
  MOD_DEF(mm,"_StdMapPythonC",module_docstring,module_methods);
  if (mm == NULL)
  { // i think this checks if the initialization of our function fgt failed
    return NULL; // exit python function call?
  }
    /* Load `numpy` functionality. */
  import_array(); // has an implicit declaration warning, probabliy resolved when python compiles
  return mm;
}
PyMODINIT_FUNC PyInit_StdMapPythonC(void)
{
  return moduleinit();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// python methods //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
This is our functional code for fgt_StdMapPythonC,
which one can write any c programming functionality for our python function 'fgt'
*/
static PyObject *parity_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Computes the Coefficient of A , such that A*gamma[x,y] is hermitian*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  parity(y,NN);
  return y_array;
}
static PyObject *iparity_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Computes the Coefficient of A , such that A*gamma[x,y] is hermitian*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  iparity(y,NN);
  return y_array;
}
static PyObject *fhgt_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Returns Fast Hermitian Gamma Transform of Input Matrix
   * Input is one matrix 
   * If the input matrix is has length of power of two, then method returns reference to matrix
   * If the input matrix is not equal to power of two, then method returns larger matrix with length power of two*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  rear(y,NN);
  int64_t i;
  for(i=0;i<NN;i++)
  {
    fht(y, NN,NN, 0,i);
  }
  parity(y,NN);
  norm(y,NN,NN);
  //Py_XDECREF(x_array);//Py_XDECREF(x_obj);//Py_XDECREF(y_array);
  return y_array;
}
static PyObject *ifhgt_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Returns Inverse Fast Hermitian Gamma Transform of Input Matrix
   * Input is one matrix 
   * If the input matrix is has length of power of two, then method returns reference to matrix
   * If the input matrix is not equal to power of two, then method returns larger matrix with length power of two*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  iparity(y,NN);
  int64_t i;
  for(i=0;i<NN;i++)
  {
    fht(y, NN,NN, 0,i);
  }
  irear(y,NN);
//  Py_XDECREF(x_obj);Py_XDECREF(x_array);Py_XDECREF(y_array);
  return y_array;
}
static PyObject *fgt_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Returns Fast Gamma Transform of Input Matrix
   * Input is one matrix 
   * If the input matrix is has length of power of two, then method returns reference to matrix
   * If the input matrix is not equal to power of two, then method returns larger matrix with length power of two*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  rear(y,NN);
  int64_t i;
  for(i=0;i<NN;i++)
  {
    fht(y, NN,NN, 0,i);
  }
  norm(y,NN,NN);
  //Py_XDECREF(x_array);//Py_XDECREF(x_obj);//Py_XDECREF(y_array);
  return y_array;
}
static PyObject *ifgt_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Returns Inverse Fast Gamma Transform of Input Matrix
   * Input is one matrix 
   * If the input matrix is has length of power of two, then method returns reference to matrix
   * If the input matrix is not equal to power of two, then method returns larger matrix with length power of two*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  int64_t i;
  for(i=0;i<NN;i++)
  {
    fht(y, NN,NN, 0,i);
  }
  irear(y,NN);
//  Py_XDECREF(x_obj);Py_XDECREF(x_array);Py_XDECREF(y_array);
  return y_array;
}
static PyObject *auto_pad_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Returns new matrix with zero padding
   * Input is one matrix 
   * If the input matrix is has length of power of two, then method returns reference to matrix
   * If the input matrix is not equal to power of two, then method returns larger matrix with length power of two*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
//  Py_XDECREF(x_obj);Py_XDECREF(x_array);Py_XDECREF(y_array);
  return y_array;
}
static PyObject *rear_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Returns Fast Rearrange Transform of Input Matrix
   * Input is one matrix 
   * If the input matrix is has length of power of two, then method returns reference to matrix
   * If the input matrix is not equal to power of two, then method returns larger matrix with length power of two*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  rear(y,NN);
//  Py_XDECREF(x_obj);Py_XDECREF(x_array);Py_XDECREF(y_array);
  return y_array;
}
static PyObject *irear_StdMapPythonC(PyObject *self, PyObject *args)
{
  /* Returns Inverse Fast Rearrange Transform of Input Matrix
   * Input is one matrix 
   * If the input matrix is has length of power of two, then method returns reference to matrix
   * If the input matrix is not equal to power of two, then method returns larger matrix with length power of two*/
  PyObject *x_obj,*x_array,*y_array;
  std::complex<double> *x,*y;
  int64_t N,NN;
  if(init_arrays(self, args, &x_obj, &x_array, &y_array, &x, &y, &N, &NN)!=0)
  {
    return NULL;
  }
  irear(y,NN);
//  Py_XDECREF(x_obj);Py_XDECREF(x_array);Py_XDECREF(y_array);
  return y_array;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// c functions definitions /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fht(std::complex<double> *x, int64_t n,int64_t NN, int64_t i1,int64_t j1)
{
  /*
  Borrowed Fast Hadamard transform, modified to evaluate different rows at a time
  https://lcav.github.io/SparseFHT
  */
  /* stopping case */
  if (n == 1)
    return;
  /* recurse */
  fht(x, n/2,NN,i1,j1);
  fht(x, n/2,NN,i1+n/2,j1);
  /* butterflies */
  int64_t i;
  for (i = 0 ; i < n/2 ; i++)
  {
    std::complex<double> tmp = x[(i+i1)+NN*(j1)];
    x[(i+i1)+NN*(j1)] = tmp +  x[(n/2+i+i1)+NN*(j1)];
    x[(n/2+i+i1)+NN*(j1)] = tmp - x[(n/2+i+i1)+NN*(j1)];
  }
  return;
}
void norm(std::complex<double> *x,std::complex<double> a,int64_t NN)
{
  /* divide 'a' into all of x's elements */
  int64_t i,j;
  for(i=0;i<NN;i++)
  {
    for(j=0;j<NN;j++)
    {
      x[i+NN*j]=x[i+NN*j]/a;
    }
  }
  return;
}
void m2g_coor(int64_t * x,int64_t* y,int64_t *mask)
{
  /*translates matrix to gamma coordinates 
   * all variables are passed by reference for speed 
   * this algrothm uses bit wise addition of +1 into base 4 number xy */
  *x=(*x)^(*y);
  *y=(*y)^(*mask);
  return;
}
void g2m_coor(int64_t * x,int64_t* y,int64_t *mask)
{
  /* translates gamma to matrix coordinates
  all variables are passed by reference for speed 
  * this algrothm uses bit wise addition of +3 into base 4 number xy */
  *x=((*x)^(*mask))^(*y); 
  *y=(*y)^(*mask);
  return;
}
void rear(std::complex<double> *mat,int64_t NN)
{
  /* in place memory swap
   * at any given time, only four elements have be to simultaneously swapped, independent of matrix size
   * the last bit of x and y will only be toggled when stepping through all swap members in the same group
   * index to different swap groups by only stepping x's and y's bits except their last bits */
  std::complex<double> temp[4] ; 
  int64_t i1,i2,index,x0,y0;
  int64_t mask = NN-1;
  for(x0=0;x0<(NN/2);x0++)
  {
    for(y0=0;y0<(NN/2);y0++)
    {
      i1=x0;
      i2=y0;
      for(index=0;index<4;index++)
      {
        temp[index] = mat[NN*i1+i2];
        m2g_coor(&i1,&i2,&mask);
      }
      for(index=0;index<4;index++)
      {
        m2g_coor(&i1,&i2,&mask);
        mat[NN*i1+i2] = temp[index] ;
      }    
    }
  }
  return;
}
void irear(std::complex<double> *mat,int64_t NN)
{
  /* in place memory swap
   * at any given time, only four elements have to simultaneously swapped, independent of matrix size
   * the last bit of x and y will only be toggled when stepping through all swap members in the same group
   * index to different swap groups by only stepping x's and y's bits except their last bits */
  std::complex<double> temp[4] ; 
  int64_t i1,i2,index,x0,y0;
  int64_t mask = NN-1;
  for(x0=0;x0<(NN/2);x0++)
  {
    for(y0=0;y0<(NN/2);y0++)
    {
      i1=x0;
      i2=y0;
      for(index=0;index<4;index++)
      {
        temp[index] = mat[NN*i1+i2];
        g2m_coor(&i1,&i2,&mask);
      }
      for(index=0;index<4;index++)
      {
        g2m_coor(&i1,&i2,&mask);
        mat[NN*i1+i2] = temp[index] ;
      }    
    }
  }
  return;
}
void matrix_pad_copy(std::complex<double> *x, std::complex<double> *y,int64_t *N,int64_t *NN)
{
  /* copiies x into y
   * takes into account different sizes between x and y
   * extra y values are filled with zeros
   * */
  int64_t i1,i2;
  for(i1=0;i1<*NN;i1++)
  {
    for(i2=0;i2<*NN;i2++)
    {
      if((i1<*N)&&(i2<*N))
      {
        y[*NN*i1+i2]=x[*N*i1+i2];
        }
      else
      {
        y[*NN*i1+i2]=0;
      }
    }		
  }
  return;
}
int64_t CountOnesFromInteger(int64_t value)
{
  /* http://codereview.stackexchange.com/questions/38182/counting-number-of-1s-and-0s-from-integer-with-bitwise-operation 
   * returns sum of all non-zero bits in value*/
  int64_t count;
  for (count = 0; value != 0; count++, value &= value-1);
  return count;
}
int64_t init_arrays(PyObject *self, PyObject *args, PyObject **x_obj, PyObject **x_array,  PyObject **y_array, std::complex<double> **x,  std::complex<double> **y, int64_t *N, int64_t *NN)
{
  /* Redundent Object Creation in similar python methods
   * double stars for pointers are double layers of address-inception, reference of a reference */
  if (!PyArg_ParseTuple(args, "O", x_obj))
  { 
    printf("_StdMapPythonC error : only one input allowed\n");
    return 1;
  }
  try
  {
    *x_array = PyArray_FROM_OTF(*x_obj, NPY_COMPLEX128, NPY_IN_ARRAY);
  }
  catch(const std::exception& e)
  {
    std::cout << e.what() << std::endl;
    Py_XDECREF(*x_array); 
    return 1;
  }
  if ((*x_array == NULL))
  { // return 1 for unable to allocate memory to new python object
    printf("_StdMapPythonC error : unable to allocate memory\n");
    Py_XDECREF(*x_array); 
    return 1;
  }
  /* check dimension number */
  if(PyArray_NDIM(*x_array) !=2)
  {
    printf("_StdMapPythonC error : array's number of dimensions is not 2 \n");    
    return 2;
  }
  /* check square-ness */
  int64_t *x_dims = PyArray_DIMS(*x_array);
  if(x_dims[0]!=x_dims[1])
  {
    //printf("_StdMapPythonC error : array is not a square matrix : x_dims=(%"PRId64",%"PRId64")\n",x_dims[0],x_dims[1]);    
    std::cout << "_StdMapPythonC error : array is not a square matrix : x_dims=(" <<
      std::to_string(x_dims[0]) << "," << std::to_string(x_dims[1]) <<
      ")\n" << std::endl;    
    return 3;
  }
  /* assign array to c pointer */
  *x = (std::complex<double> *)PyArray_DATA(*x_array);
  *N = (int64_t)PyArray_DIM(*x_array, 0);
  int64_t n = (int64_t)ceil(log2(*N)); 
  *NN = 1<<n;
  if(*N!=*NN)
  { // create y as a new data object of greater size
    npy_intp dims[2]; //this could cause a problem ///////////////////////////////////////////////////////////////
    dims[0] = *NN;
    dims[1] = *NN;
    *y_array = PyArray_SimpleNew(2, dims, NPY_CDOUBLE);
    if ((*y_array == NULL))
    { 
      printf("_StdMapPythonC error : unable to allocate memory\n");
      Py_XDECREF(*y_array); 
      Py_XDECREF(*x_array); 
      return 1;
    }
    *y = (std::complex<double> *)PyArray_DATA(*y_array);
    matrix_pad_copy(*x, *y,N,NN);
  }
  else
  { // set y to equal x
    *y_array = *x_array;
    *y=*x;
  }
  return 0;
}
void parity(std::complex<double> *z,int NN)
{
  int64_t i1,i2,parity;
  for(i1=0;i1<NN;i1++)//vary irreducibles
  {
    for(i2=0;i2<NN;i2++)//vary groups
    {
      /*     i1 = |    A    (X)    D    (X)    A    (X)    D    | // rearrangement of matrix elements, D = diagonal pauli matrix, A=anti-diagonal pauli matrix
       *     i2 = | (+1 -1) (X) (+1 +1) (X) (+1 +1) (X) (+1 -1) | // hadamard transform
       *     i1 = |    A    (X)    D    (X)    A    (X)    D    | // rearrangement of matrix elements, D = diagonal pauli matrix, A=anti-diagonal pauli matrix
       *     i2 = | (+1 -1) (X) (+1 +1) (X) (+1 +1) (X) (+1 -1) | // hadamard transform
       * parity = |  (1)(1)  +   (0)(0)  +   (1)(0)  +   (0)(1) | = 1 : odd is anti-hermitian , even is hermitian              
* xy -->(+1)--> mn
* 00            01
* 01            10
* 10            11
* 11            00
* */
      parity=CountOnesFromInteger(i1&i2);
      parity=parity%4; // mod 4 keeps track of sign change, mod 2 does not
      if     (parity ==0)
      {
        z[NN*i1+i2]*=1;
      }
      else if(parity ==1)
      {
        z[NN*i1+i2]*=-imaginary_number;
      }
      else if     (parity ==2)
      {
        z[NN*i1+i2]*=-1;
      }
      else// (parity ==3)
      {
        z[NN*i1+i2]*=imaginary_number;
      }
    }		
  }
  return;
}
void iparity(std::complex<double> *z,int NN)
{
  int64_t i1,i2,parity;
  for(i1=0;i1<NN;i1++)//vary irreducibles
  {
    for(i2=0;i2<NN;i2++)//vary groups
    {
      /*     i1 = |    A    (X)    D    (X)    A    (X)    D    | // rearrangement of matrix elements, D = diagonal pauli matrix, A=anti-diagonal pauli matrix
       *     i2 = | (+1 -1) (X) (+1 +1) (X) (+1 +1) (X) (+1 -1) | // hadamard transform
       * parity = |  (1)(1)  +   (0)(0)  +   (1)(0)  +   (0)(1) | = 1 : odd is anti-hermitian , even is hermitian              */
      parity=CountOnesFromInteger(i1&i2);
      parity=parity%4; // mod 4 keeps track of sign change, mod 2 does not
      if     (parity ==0)
      {
        z[NN*i1+i2]*=1;
      }
      else if(parity ==1)
      {
        z[NN*i1+i2]*=imaginary_number;
      }
      else if(parity ==2)
      {
        z[NN*i1+i2]*=-1;
      }
      else// (parity ==3)
      {
        z[NN*i1+i2]*=-imaginary_number;
      }
    }		
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Things to Do for Completition ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * memory leakage check? Use proper dereferencing, create a dereferencing function
 * dereferening any of the python objects caused segmentation faults?
 * Debug Depriation Warning
 * transform from pauli base to gamma base?
 * raise python errors and their messages up function levels properly
 * Write a manaul for installing all required packages, provide installation instructions if import fails
 * * sudo apt-get purge python-numpy # if all ready installed and out of date
 * * sudo apt-get install python-pip # if not installed
 * * sudo apt-get install python-dev # if not installed
 * * sudo apt-get install g++ # if not installed
 * * sudo pip install numpy # if not installed or not update to date
 * 
 * */
