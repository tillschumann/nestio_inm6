#include "Python.h"
#include "nestio_func.h"

static PyObject *
RandomGenerator_normal(PyObject *self, PyObject *args, PyObject *keywds)
{
    int n;
    double mean;
    double variance;

    static char *kwlist[] = {"number", "mean", "variance", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "idd", kwlist,&n, &mean, &variance))
        return NULL;
    
    PyObject *lst_values = PyList_New(n);
    if (!lst_values)
      return NULL;

    nestio::StandardDistribution dist(mean,variance);
    for (int i=0; i<n; i++) {
      double v = dist.getValue();
      PyObject *num = PyFloat_FromDouble(v);
      if (!num) {
        Py_DECREF(lst_values);
        return NULL;
      }
      PyList_SET_ITEM(lst_values, i, num);   // reference to num stolen
    }

    return lst_values;
}

static PyObject *
RandomGenerator_poisson(PyObject *self, PyObject *args, PyObject *keywds)
{
    int n;
    double lambda;

    static char *kwlist[] = {"number", "lambda", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "id", kwlist,&n, &lambda))
        return NULL;
    
    PyObject *lst_values = PyList_New(n);
    if (!lst_values)
      return NULL;

    nestio::PoissonDistribution dist(lambda);
    for (int i=0; i<n; i++) {
      double v = dist.getValue();
      PyObject *num = PyFloat_FromDouble(v);
      if (!num) {
        Py_DECREF(lst_values);
        return NULL;
      }
      PyList_SET_ITEM(lst_values, i, num);   // reference to num stolen
    }

    return lst_values;
}

static PyObject *
RandomGenerator_binominal(PyObject *self, PyObject *args, PyObject *keywds)
{
    int n;
    double t;
    double p;

    static char *kwlist[] = {"number", "t", "p", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "idd", kwlist,&n, &t, &p))
        return NULL;
    
    PyObject *lst_values = PyList_New(n);
    if (!lst_values)
      return NULL;

    nestio::BinominalDistribution dist(t,p);
    for (int i=0; i<n; i++) {
      double v = dist.getValue();
      PyObject *num = PyFloat_FromDouble(v);
      if (!num) {
        Py_DECREF(lst_values);
        return NULL;
      }
      PyList_SET_ITEM(lst_values, i, num);   // reference to num stolen
    }

    return lst_values;
}

static PyMethodDef RandomGenerator_methods[] = {
    /* The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and RandomGenerator_generate() takes
     * three.
     */
    {"normal", (PyCFunction)RandomGenerator_normal, METH_VARARGS | METH_KEYWORDS,
     "Generate random variables"},
    {"poisson", (PyCFunction)RandomGenerator_poisson, METH_VARARGS | METH_KEYWORDS,
     "Generate random variables"},
    {"binominal", (PyCFunction)RandomGenerator_binominal, METH_VARARGS | METH_KEYWORDS,
     "Generate random variables"},
    {NULL, NULL, 0, NULL}   /* sentinel */
};

PyMODINIT_FUNC initRandomGenerator(void)
{
    Py_InitModule("RandomGenerator", RandomGenerator_methods);
}