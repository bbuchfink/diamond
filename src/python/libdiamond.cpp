#include <Python.h>
#include "../run/main.h"

/**
 * Run diamond by cmd options
*/
static PyObject* method_main(PyObject* self, PyObject* args)
{
    try {
        int size = PyTuple_Size(args);
        int argc = size + 1;
        const char* argv[argc];
        argv[0] = "diamond";
        for (int i = 0; i < size; i++) {
            PyObject* item = PyTuple_GetItem(args, i);
            if (!PyArg_Parse(item, "s", &(argv[i + 1]))) {
                return NULL;
            }
        }
        int status = diamond(argc, argv);
        return Py_BuildValue("i", status);
	}
	catch(const std::bad_alloc &e) {
        PyErr_SetString(PyExc_MemoryError, e.what());
        return NULL;
	} catch(const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
    catch(...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception occurred from Cpp");
        return NULL;
    }
}

/**
 * Return diamond's version
*/
static PyObject* method_version(PyObject *self, PyObject*args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }
    char* version = (char *)Const::version_string;
    return Py_BuildValue("s", version);
}

static PyMethodDef libdiamond_methods[] = {
    {
        "main", 
        method_main, 
        METH_VARARGS, 
        "Run diamond by its command options. For option details, just pass a 'help' argument"},
    {
        "version", 
        method_version, 
        METH_VARARGS, 
        "Return the version of diamond."
    },
    // the last one used just to tell Python the end of method list.
    {
        NULL, 
        NULL, 
        0, 
        NULL}
};

static struct PyModuleDef libdiamond_module {
    PyModuleDef_HEAD_INIT,
    "diamondpy",
    "Diamond's python wrapper module",
    -1,
    libdiamond_methods
};

PyMODINIT_FUNC PyInit_libdiamond(void)
{
    return PyModule_Create(&libdiamond_module);
}


