#include <Python.h>

#include<stdio.h>
#include<string.h>

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
    Py_DECREF(lst);
    free(lPhases);
    free(lvalues);
    return result;
}

static PyObject * phasogramraw(PyObject *self, PyObject *args)
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
           lPhases[j] += lvalues[i]; 
	}
    }

    PyObject *lst = PyList_New(window_size);
    PyObject *num;

    for(i=0;i<window_size;i++){
        num=PyFloat_FromDouble(lPhases[i]);
	PyList_SET_ITEM(lst,i,num);
    }

    PyObject *result= Py_BuildValue("O",lst);
    Py_DECREF(lst);
    free(lPhases);
    free(lvalues);
    return result;
}

static PyObject * dinucfrequency(PyObject *self, PyObject *args)
{
    int distance;
    PyObject *sequences;
    
    if (!PyArg_ParseTuple(args, "iO", &distance, &sequences)){
	  PyErr_SetString(PyExc_TypeError, "dinucfrequency parameters errors, expected parameters: int, [string]");
	  return NULL;
    }

   int i;
   int j;
   double *ATheap = malloc((2*distance+1) * sizeof(double));
   double *GCheap = malloc((2*distance+1) * sizeof(double));
   for(i=0;i<(2*distance+1);i++){
       ATheap[i] = 0.0;
       GCheap[i] = 0.0;
   }
   double *ATvalues = malloc((2*distance+1) * sizeof(double));
   double *GCvalues = malloc((2*distance+1) * sizeof(double));
   int nb = PyList_Size(sequences);

   for(j=0;j<nb;j++){
       
      char *sequence = PyString_AsString(PyList_GetItem(sequences,j));

      int size = strlen(sequence);
      int middle = size/2;
      char revsequence[size];
      for(i=0;i<size;i++) {
          revsequence[size-1-i]=sequence[i]; 
      }

      for(i=0;i<(2*distance+1);i++){
          ATvalues[i] = 0.0;
          GCvalues[i] = 0.0;
      }
      for(i = middle; i<size-1;i++) {
          char substr[3];
          char revsubstr[3];
          strncpy(substr, sequence + i,2);
          substr[sizeof(substr)-1] = '\0';
          strncpy(revsubstr, revsequence + i,2);
          revsubstr[sizeof(revsubstr)-1] = '\0';

          if ((strcmp(substr,"AA") == 0) ||  (strcmp(substr,"TT") == 0) || (strcmp(substr,"AT") == 0) || (strcmp(substr,"TA") == 0))  {
              ATvalues[i]++;
          }
          if ((strcmp(revsubstr,"AA") == 0) ||  (strcmp(revsubstr,"TT") == 0) || (strcmp(revsubstr,"AT") == 0) || (strcmp(revsubstr,"TA") == 0))  {
              ATvalues[middle-(i-middle)]++;
          } 
          if ((strcmp(substr,"GG") == 0) ||  (strcmp(substr,"CC") == 0) || (strcmp(substr,"GC") == 0) || (strcmp(substr,"CG") == 0))  {
              GCvalues[i]++;
          }
          if ((strcmp(revsubstr,"GG") == 0) ||  (strcmp(revsubstr,"CC") == 0) || (strcmp(revsubstr,"GC") == 0) || (strcmp(revsubstr,"CG") == 0))  {
              GCvalues[middle-(i-middle)]++;
          } 
      }

       ATvalues[middle] = ATvalues[middle]/2.0;
       GCvalues[middle] = GCvalues[middle]/2.0;

      for(i=0;i<(2*distance+1);i++){
          ATheap[i] += ATvalues[i];
          GCheap[i] += GCvalues[i];
      }
   }

    PyObject *ATlst = PyList_New((2*distance+1));
    PyObject *GClst = PyList_New((2*distance+1));
    PyObject *num;

    for(i=0;i<(2*distance+1);i++){
        num=PyFloat_FromDouble(ATheap[i]);
	PyList_SET_ITEM(ATlst,i,num);
        num=PyFloat_FromDouble(GCheap[i]);
        PyList_SET_ITEM(GClst,i,num);
    }

    PyObject *result= Py_BuildValue("OO",ATlst,GClst);
    Py_DECREF(ATlst);
    Py_DECREF(GClst);
    free(ATvalues);
    free(GCvalues);
    free(ATheap);
    free(GCheap);
    return result;
}

static char phasogram_doc[] = "phasogram(), compute phasogram of input sequence\n";

static char dinucfrequency_doc[] = "dinucfrequency(), doc\n"; 

static PyMethodDef nucleosome_funcs[] = {
    {"phasogram", (PyCFunction)phasogram, METH_VARARGS, phasogram_doc},
    {"phasogramraw", (PyCFunction)phasogramraw, METH_VARARGS, phasogram_doc},
    {"dinucfrequency", dinucfrequency, METH_VARARGS, dinucfrequency_doc},
    {NULL, NULL, 0, NULL} /* sentinel */ 
};	



PyMODINIT_FUNC initCMSTS(void)
{
   Py_InitModule3("CMSTS", nucleosome_funcs, "CMSTS module"); 

}
