/*
This code was taken from Dan Foreman-Mackey and modified.
Source Code : https://gist.github.com/dfm/3247796
Source Documentation : http://dan.iel.fm/posts/python-c-extensions/
*/
#include <Python.h>
#include <stdio.h>
#include <vector>
void rfstoh(int num, int *total , char frompeg, char topeg, char auxpeg);
/*
Docstrings module_docstring
Below establishes documentation for our fstoh function
In python executing 'fstoh?' or 'help(fstoh)' will display 'fstoh_docstring' in the terminal
 */
static char module_docstring[] =
  "This module provides an interface for calculating the Tower of Hanoi.";
static char fstoh_docstring[] =
  "From wiki : The objective of the puzzle is to move the entire stack to another rod, obeying the following simple rules:\n(1) Only one disk can be moved at a time.\n(2) Each move consists of taking the upper disk from one of the stacks and placing it on top of another stack i.e. a disk can only be moved if it is the uppermost disk on a stack.\n(3)No disk may be placed on top of a smaller disk."; 
/*
Available functions 
This is a declaration of a structure type.
This is our declaration of our 'fstoh_fstoh' function for python. 
It specifies all input and output objects of our function. 
The name that we've given to the function (fstoh_fstoh) is also a matter of convention and nessaray because {module_name}_{function_name}. 
The name for this c file given as '_fstoh' is a matter of convention and necessary. 
Calling the function in python is _fstoh.fstoh(), where _fstoh is the name of the module and fstoh is the name of the function. 
the '_' for c is equivalent for the '.' in python in regards to class extensions 
*/
static PyObject *fstoh_fstoh(PyObject *self, PyObject *args); 
/*
Module specification 
This specifies all functional members of this python module
In our case there is only one, which is 'fstoh_fstoh' as defined above
*/
static PyMethodDef module_methods[] = {
  {"fstoh", fstoh_fstoh, METH_VARARGS, fstoh_docstring},
  {NULL, NULL, 0, NULL}
};
/*
Initialize the module. 
This initializes are module containing our function.
This initializing module function must be named as 'init'+'ourmodulename',
which in our case will be 'init'+'_fstoh' = 'init_fstoh' .
*/
PyMODINIT_FUNC init_fstoh(void)
{
  PyObject *mm = Py_InitModule3("_fstoh", module_methods, module_docstring);
  if (mm == NULL) // i think this checks if the initialization of our function fstoh failed
    return; // exit python function call?
}
/*
This is our functional code for fstoh_fstoh,
which one can write any c programming functionality for our python function 'fstoh' 
*/
static PyObject *fstoh_fstoh(PyObject *self, PyObject *args)
{
  double m ; // declares a c variable to hold our python input variable
  /*
  Parse the input tuple. From the python documentation :
  PyArg_ParseTuple will allocate a buffer of the needed size, copy the encoded data into this buffer and adjust *buffer to reference the newly allocated storage. 
  The caller is responsible for calling PyMem_Free() to free the allocated buffer after use (Note : only for pointers, m is not a pointer) . 
  I think it returns false if it was unsuccesful in allocating space for the variables
  */
  if (!PyArg_ParseTuple(args, "d", &m)) // passing address of python variables to c variables?
    return NULL;
  /*
  This is our main programing. It calls the c function 'rfstoh'.
  */
  int *total; // an internal c pointer to count the number steps in pur puzzel
	total = (int *)malloc(sizeof(int));    // Allocate an int pointee,
      *total = 0 ;// pointer to a block of memory to calculate the total number of steps in loop assigned a value of zero
  printf("The sequence of moves involved in the Tower of Hanoi are :\n");
  rfstoh(m,total, 'A', 'C', 'B');
  printf("\n Total number of steps : ");
  /*
  Build the output tuple 
  'PyObject *ret' tells python the pointer containing our output result
  'i' specifices the value of '*total' is an integer
  */
  PyObject *ret = Py_BuildValue("i", *total);
  free(total);//clean up memory leak
  return ret;
}
void rfstoh(int num, int *total , char frompeg, char topeg, char auxpeg)
{
  // Taken from : http://www.sanfoundry.com/c-program-tower-of-hanoi-using-recursion/
  if (num == 1)
  {
      printf("\n Move disk 1 from peg %c to peg %c", frompeg, topeg);
      return;
  }
	*total += 1 ; // computes total number of steps
  rfstoh(num - 1,total, frompeg, auxpeg, topeg);
  printf("\n Move disk %d from peg %c to peg %c", num, frompeg, topeg);
  rfstoh(num - 1,total, auxpeg, topeg, frompeg);
}
