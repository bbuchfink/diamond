#include <Python.h>
#include "../run/main.h"

static PyObject* method_run(PyObject* self, PyObject* args)
{
    const char* db, *query, *out, *command;
    // int status = diamond(c_argc, c_args);
    if (!PyArg_ParseTuple(args, "ssss", &db, &query, &out, &command)) {
        return NULL;
    }
    try {
        // std::cout << command << std::endl;
        const char * argv[] = {"diamond", command, "-d", db, "-q", query, "-o", out};
        config = Config(8, argv);
        switch (config.command) {
        case config.blastp:
        case config.blastx:
            Search::run();
            // printf("Diamond finished!\n");
            std::cout << "Finished!" << std::endl;
            break;
        default:
            PyErr_SetString(PyExc_NotImplementedError, "Command not supported yet.");
            return NULL;
        }
    }
    catch(const std::bad_alloc &e) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate sufficient memory. Please refer to the manual for instructions on memory usage.");
        return NULL;
    }
    catch(const std::exception& e) {
        string err = e.what();
        if (err.find("Invalid command", 0) == 0) {
            PyErr_SetString(PyExc_NotImplementedError, "Unsupported command.");
        }
        else if (err.find("Error calling stat on file ", 0) == 0) {
            string file = err.substr(27);
            PyErr_SetString(PyExc_FileNotFoundError, file.c_str());
        }
        else {
            PyErr_SetString(PyExc_RuntimeError, e.what());
            std::cout << 1 << std::endl;
        }
        return NULL;
    }
    catch (const std::string& e) {
        PyErr_SetString(PyExc_RuntimeError, e.c_str());
        std::cout << 2 << std::endl;
        return NULL;
    }
    catch (const int& e) {
        PyErr_SetString(PyExc_RuntimeError, "Error code returned.");
        return NULL;
    }
    catch(...) {
        std::cout << "Hello World" << std::endl;
        std::exception_ptr e = std::current_exception();
        std::clog << (e ? e.__cxa_exception_type()->name() : "null") << std::endl;
        PyErr_SetString(PyExc_RuntimeError, "Exception of unknown type!");
        // std::cerr << 
        return NULL;
    }
}

static PyObject* method_version(PyObject *self, PyObject*args)
{   
    // 无参数
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }
    char* version = (char *)Const::version_string;
    return Py_BuildValue("s", version);
}


static PyMethodDef libdiamond_methods[] = {
    {"run", method_run, METH_VARARGS, "Run blastp or blastx."},
    {"version", method_version, METH_VARARGS, "Return internal version of diamond"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef libdiamond_module {
    PyModuleDef_HEAD_INIT,
    "diamondpy",
    "Diamond module",
    -1,
    libdiamond_methods
};

PyMODINIT_FUNC PyInit_libdiamond(void)
{
    return PyModule_Create(&libdiamond_module);
}


