#include <Python.h>

static PyObject * phasogram(PyObject *self, PyObject *args)
{
    PyObject *chr_values;
    int window_size;
    int i;
    int j;
    
    if (!PyArg_ParseTuple(args, "iO", &window_size, &chr_values)){
	  PyErr_SetString(PyExc_TypeError, "phasogram parameters errors, expected parameters: int, [float]");
	  return NULL;
    }

    int nb = PyList_Size(chr_values);

    double *lPhases = malloc(window_size * sizeof(double));

    for (i=0;i<window_size;i++) {
        lPhases[i] = 0.0;
    }

    double *lvalues = malloc(nb * sizeof(double));

    for(i=0;i<nb;i++){
        lvalues[i] = PyFloat_AsDouble(PyList_GetItem(chr_values,i));
    }

    for (i=0; i<(nb-window_size); i++) {
        for (j=0; j<window_size; j++){
           lPhases[j] += lvalues[i] * lvalues[i+j];
	}
    }

    PyObject *lst = PyList_New(window_size);
    PyObject *num;

    for(i=0;i<window_size;i++){
        num=PyFloat_FromDouble(lPhases[i]);
	PyList_SET_ITEM(lst,i,num);
    }

    PyObject *result= Py_BuildValue("O",lst);
    free(lPhases);
    free(lvalues);
    return result;
}

static char phasogram_doc[] = "phasogram(), compute phasogram of input sequence\n"; 

static PyMethodDef nucleosome_funcs[] = {
    {"phasogram", (PyCFunction)phasogram, METH_VARARGS, phasogram_doc},
    {NULL, NULL, 0, NULL} /* sentinel */ 
};	



PyMODINIT_FUNC initCMSTS(void)
{
   Py_InitModule3("CMSTS", nucleosome_funcs, "CMSTS module"); 

}
