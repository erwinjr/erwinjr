/* modulehello.c */

#include <Python.h>
#include <stdio.h>

/* The C definition of our function to print Hello World. By convention
the name is <module name>_<function name>, but this isn't necessary. */

static PyObject*
hello_sayHello(PyObject *self, PyObject *args) {

    printf("Hello World using printf");

    /* We wish for our function to return None, the Python
    equivalent of returning void. We have to do the following
    to return None. */

    Py_INCREF(Py_None);
    return Py_None; 
}

/* This array lists all of the methods we are putting into our module. Take
note of the sentinel value at the end to indicate the ending of the array. */

static PyMethodDef HelloMethods[] = {
    {"sayHello", hello_sayHello, METH_VARARGS, "Print Hello World"},
    {NULL, NULL, 0, NULL} /* The sentinel value. */
};


/* This declares set-up information for our module.*/

static struct PyModuleDef hellomodule = {

    PyModuleDef_HEAD_INIT,
    "hello",
    NULL, /*This is for documentation, which we won't use; so it is NULL. */
    -1,
    HelloMethods
};

/* Function to initialize the module. Note the necessary name format of
PyInit_<module name>. */

PyMODINIT_FUNC
PyInit_hello(void) {

    return PyModule_Create(&hellomodule);
}
