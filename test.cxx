
#include <Python.h> 

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include<vector>
#include <iostream>

double RosenBrock(const double *xx ) // real function
{
  const double x = xx[0];
  const double y = xx[1];
  const double tmp1 = y-x*x;
  const double tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

int main(int argc, char *argv[])
{

auto spam_system = [](PyObject *self, PyObject *args)->PyObject *
{
    
    std::cout<<PyTuple_Size(args)<<std::endl;
    PyArrayObject *arr = (PyArrayObject *)PyTuple_GetItem(args, 0);
    
    auto params = (double *)PyArray_DATA(arr);
    std::cout<<params[0]<<" "<<params[0]<<std::endl;
    auto r = RosenBrock(params);
    std::cout<<r<<std::endl;

    return PyFloat_FromDouble(r);
};

static PyObject *SpamError;
static PyMethodDef SpamMethods[] = {
    {"system",  spam_system, METH_VARARGS,
     "Execute a shell command."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef spammodule = {
    PyModuleDef_HEAD_INIT,
    "spam",   /* name of module */
    "spam_doc", /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    SpamMethods
};

    auto PyInit_spam = [](void)->PyObject *
    {
    PyObject *m;

    m = PyModule_Create(&spammodule);
    if (m == NULL)
        return NULL;

    SpamError = PyErr_NewException("spam.error", NULL, NULL);
    Py_XINCREF(SpamError);
    if (PyModule_AddObject(m, "error", SpamError) < 0) {
        Py_XDECREF(SpamError);
        Py_CLEAR(SpamError);
        Py_DECREF(m);
        return NULL;
    }

    return m;  
    };
    

    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }    
/* Add a built-in module, before Py_Initialize */
    if (PyImport_AppendInittab("spam", PyInit_spam) == -1) {
        fprintf(stderr, "Error: could not extend in-built modules table\n");
        exit(1);
    }
    Py_SetProgramName(program);

    Py_Initialize();
    _import_array();
    PyObject *pmodule = PyImport_ImportModule("spam");
    if (!pmodule) {
        PyErr_Print();
        fprintf(stderr, "Error: could not import module 'spam'\n");
    }

   PyObject *mpmodule = PyModule_GetDict(pmodule);
   PyObject *fsystem = PyDict_GetItemString(mpmodule, "system");
   if (fsystem) Py_INCREF(fsystem);

    double Array [] = {0.1,0.1};
    npy_intp dims[1] = { 2 };
    PyObject *py_array;
    py_array = PyArray_SimpleNewFromData(1, dims, NPY_CDOUBLE, Array);
   auto pycode = Py_BuildValue("(O)",py_array);
   PyObject *result = PyObject_CallObject(fsystem,pycode);
   Py_DECREF(pycode);



    PyObject *pName = PyUnicode_FromString("scipy.optimize");
   // Import the file as a Python module.
   // return object is a new reference !
   PyObject *fModuleOpt = PyImport_Import(pName);
   //Py_DECREF(pName);
   if (!fModuleOpt) {
    std::cerr<<"--- ERROR: Can't import scipy "<< std::endl;
   }
   PyObject *mDict = PyModule_GetDict(fModuleOpt);
   PyObject *fMinimize = PyDict_GetItemString(mDict, "minimize");
   if (fMinimize) Py_INCREF(fMinimize);
   pycode = Py_BuildValue("(OO)",fsystem,py_array);
   result = PyObject_CallObject(fMinimize,pycode);
   Py_DECREF(pycode);

  // PyObject *pycode = Py_BuildValue("()", "root_scipy.wrap_function");
   //PyObject *result = PyObject_CallObject(fMinimize, pycode);
   //Py_DECREF(pycode);
   //pName = PyUnicode_FromString("root_scipy");
   //PyObject *fRootScipy = PyImport_Import(pName);
   //PyRun_String("root_scipy.test('a')", Py_single_input, NULL, NULL);
   //Py_DECREF(pName);
   //if (!fModuleOpt) {
   // std::cerr<<"--- ERROR: Can't import root_scipy "<< std::endl;
   //}
   //mDict = PyModule_GetDict(fRootScipy);
   //PyObject *fTest = PyDict_GetItemString(mDict, "test");
   //pycode = Py_BuildValue("()");
   //PyObject *result = PyObject_CallObject(fTest,NULL);
   //Py_DECREF(pycode);

   if (fMinimize) Py_DECREF(fMinimize);
    PyMem_RawFree(program);
    Py_Finalize();
    return 0;
}