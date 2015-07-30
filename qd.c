#include <Python.h>
#include <qd/c_qd.h>
#include <qd/c_dd.h>
#include <qd/fpu.h>

// TODO: check whether it is orthodox to call PyType_GenericNew for creating
// e ne winstance of the type (by skipping the useless initialization step)

/* using isnan(x) from math.h but linking to the library is not needed */
#include <math.h>

/* ========================================================================== */
static PyTypeObject PyQDTypeObjectType;
static PyTypeObject PyDDTypeObjectType;
typedef struct {
    PyObject_HEAD;
    double content_data[2];
} PyDDTypeObject;
typedef struct {
    PyObject_HEAD;
    double content_data[4];
} PyQDTypeObject;
unsigned int fpu_state;
/* ========================================================================== */

static PyObject *fpu_init(PyObject *s, PyObject *a) { fpu_fix_start(&fpu_state); Py_RETURN_NONE; }
static PyObject *fpu_restore(PyObject *s, PyObject *a) { fpu_fix_end(&fpu_state); Py_RETURN_NONE; }

#define DD_GENERICWRAPPER1(f,self) \
            PyObject *o; \
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL); \
            f(self->content_data, ((PyDDTypeObject *) o)->content_data); \
            return o;

#define DD_GENERICWRAPPER2(f,self) \
            PyObject *a,*b; \
            a = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL); \
            b = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL); \
            f(self->content_data, ((PyDDTypeObject *) a)->content_data, \
                                  ((PyDDTypeObject *) b)->content_data); \
            return PyTuple_Pack(2, a,b);
#define DD_GENERICWRAPPER3(f,self, args) \
            PyObject *o; \
            int n; \
            if (!PyArg_ParseTuple(args, "i", &n)) return NULL; \
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL); \
            f(self->content_data, n, ((PyDDTypeObject *) o)->content_data); \
            return o;

static int DD_init(PyDDTypeObject *self, PyObject *args) {
    PyObject *a = NULL;
    PyObject *tmp;

    if (!PyArg_ParseTuple(args, "|O", &a)) return -1;

    if(a==NULL) return 0; /* 0 argument */

    if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) { /* qd.DD */
        self->content_data[0] = ((PyDDTypeObject *)a)->content_data[0];
        self->content_data[1] = ((PyDDTypeObject *)a)->content_data[1];
        return 0;
    }

    if(PyFloat_Check(a)) { /* Float */
        self->content_data[0] = PyFloat_AS_DOUBLE(a);
        self->content_data[1] = 0.0;
        return 0;
    }

    if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) { /* qd.QD */
        self->content_data[0] = ((PyQDTypeObject *)a)->content_data[0];
        self->content_data[1] = ((PyQDTypeObject *)a)->content_data[1];
        return 0;
    }

    if(PyInt_Check(a)||PyLong_Check(a)) { /* Int/Long */
        tmp = PyObject_Str(a);
        c_dd_read(PyString_AS_STRING(tmp), self->content_data);
        Py_DECREF(tmp);
        return 0;
    }

    if(PyString_Check(a)) { /* String */
        c_dd_read(PyString_AS_STRING(a), self->content_data);
        if( isnan(self->content_data[0] )) {
            PyErr_SetString(PyExc_ValueError,"Converting string to DD gave a NaN Value");
            return -1;
        }
        return 0;
    }

    tmp = PyNumber_Float(a);
    if (tmp != NULL) { /* __float__ */
        self->content_data[0] = PyFloat_AS_DOUBLE(tmp);
        self->content_data[1] = 0.0;
        Py_DECREF(tmp);
        return 0;
    }

    return -1;
}

static PyObject *DD_reset(PyDDTypeObject *self, PyObject *args) {
    (void) DD_init(self, args);
    Py_RETURN_NONE;
}

static PyObject *DD_repr(PyDDTypeObject *self) {
    char repr[40];
    c_dd_swrite( self->content_data, 32, repr, 40 );
    return PyString_FromString(repr);
}

static PyObject *DD_sqrt(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_sqrt, self) }
static PyObject *DD_sqr(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_sqr, self) }
static PyObject *DD_abs(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_abs, self) }
static PyObject *DD_neg(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_neg, self) }
static PyObject *DD_nint(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_nint, self) }
static PyObject *DD_aint(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_aint, self) }
static PyObject *DD_floor(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_floor, self) }
static PyObject *DD_ceil(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_ceil, self) }
static PyObject *DD_exp(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_exp, self) }
static PyObject *DD_log(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_log, self) }
static PyObject *DD_log10(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_log10, self) }
static PyObject *DD_sin(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_sin, self) }
static PyObject *DD_cos(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_cos, self) }
static PyObject *DD_tan(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_tan, self) }
static PyObject *DD_asin(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_asin, self) }
static PyObject *DD_acos(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_acos, self) }
static PyObject *DD_atan(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_atan, self) }
static PyObject *DD_sinh(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_sinh, self) }
static PyObject *DD_cosh(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_cosh, self) }
static PyObject *DD_tanh(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_tanh, self) }
static PyObject *DD_asinh(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_asinh, self) }
static PyObject *DD_acosh(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_acosh, self) }
static PyObject *DD_atanh(PyDDTypeObject *self) { DD_GENERICWRAPPER1(c_dd_atanh, self) }
static PyObject *DD_sincos(PyDDTypeObject *self) { DD_GENERICWRAPPER2(c_dd_sincos, self) }
static PyObject *DD_sincosh(PyDDTypeObject *self) { DD_GENERICWRAPPER2(c_dd_sincosh, self) }
static PyObject *DD_npwr(PyDDTypeObject *self, PyObject *args) {
    DD_GENERICWRAPPER3(c_dd_npwr, self, args) }
static PyObject *DD_nroot(PyDDTypeObject *self, PyObject *args) {
    DD_GENERICWRAPPER3(c_dd_nroot, self, args) }

static PyObject *DD_atan2(PyDDTypeObject *self, PyObject *args) {
    PyObject *x;
    PyObject *y;
    PyObject *o;

    if (!PyArg_ParseTuple(args, "OO", &x, &y)) return NULL;

    if(!PyObject_TypeCheck(x, &PyDDTypeObjectType)) {
        PyErr_SetString(PyExc_TypeError,"Wrong type for argument 1 of qd.DD.atan2()");
        return NULL;
    }
    if(!PyObject_TypeCheck(y, &PyDDTypeObjectType)) {
        PyErr_SetString(PyExc_TypeError,"Wrong type for argument 2 of qd.DD.atan2()");
        return NULL;
    }

    o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
    c_dd_atan2( ((PyDDTypeObject *) x)->content_data,
                ((PyDDTypeObject *) y)->content_data,
                ((PyDDTypeObject *) o)->content_data );

    return o;
}

static PyObject *DD_float(PyDDTypeObject *self) {
    return PyFloat_FromDouble( self->content_data[0] );
}

static PyObject *DD_set_random(PyDDTypeObject *self) {
    c_dd_rand( self->content_data );
    Py_RETURN_NONE;
}

static PyObject *DD_clone(PyDDTypeObject *self) {
    PyObject *o;
    o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
    ((PyDDTypeObject *) o)->content_data[0] = self->content_data[0];
    ((PyDDTypeObject *) o)->content_data[1] = self->content_data[1];
    return o;
}

static PyObject *DD_get_address(PyDDTypeObject *self) {
    return PyInt_FromLong( (long) self->content_data );
}

static PyObject *DD_get_content(PyDDTypeObject *self) {
    return PyTuple_Pack(2, 
           PyFloat_FromDouble( self->content_data[0] ),
           PyFloat_FromDouble( self->content_data[1] ) );
}

static PyObject *DD_pow(PyDDTypeObject *self, PyObject *a, PyObject *b) {
    PyObject *o;

    if (b != Py_None) {
        PyErr_SetString(PyExc_NotImplementedError,"Not implemented operation");
        return NULL;
    }

    if((!PyInt_Check(a))&&(!PyLong_Check(a))) {
        PyErr_SetString(PyExc_TypeError,"Wrong type for exponent");
        return NULL;
    }

    o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
    c_dd_npwr(self->content_data, PyInt_AS_LONG(a),
            ((PyDDTypeObject *) o)->content_data);
    return o;
}

static PyObject *DD_add(PyObject *a, PyObject *b) {
    PyDDTypeObject *left, *right;
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
        left = (PyDDTypeObject *) a;
        right = (PyDDTypeObject *) b;
    } else {
        left = (PyDDTypeObject *) b;
        right = (PyDDTypeObject *) a;
    }

    if(PyObject_TypeCheck(right, &PyDDTypeObjectType)) {
        o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
        c_dd_add(left->content_data, right->content_data,
                        ((PyDDTypeObject *) o)->content_data);
        return o;
    }

    if(PyFloat_Check(right)) {
        o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
        c_dd_add_dd_d(left->content_data, PyFloat_AS_DOUBLE(right),
                        ((PyDDTypeObject *) o)->content_data);
        return o;
    }

    if(PyObject_TypeCheck(right, &PyQDTypeObjectType)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_add_dd_qd(left->content_data, ((PyQDTypeObject *)right)->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyInt_Check(right)||PyLong_Check(right)) { /* Int/Long */
        o = PyTuple_Pack(1, right);
        right = (PyDDTypeObject *) PyObject_CallObject( 
                           (PyObject *)&PyDDTypeObjectType, o);
        Py_DECREF(o);

        o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
        c_dd_add(left->content_data, right->content_data,
                        ((PyDDTypeObject *) o)->content_data);
        Py_DECREF(right);
        return o;
    }

    return Py_NotImplemented;
}

static PyObject *DD_sub(PyObject *a, PyObject *b) {
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
        if(PyObject_TypeCheck(b, &PyDDTypeObjectType)) {
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_sub(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            return o;
        }

        if(PyFloat_Check(b)) {
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_sub_dd_d(((PyDDTypeObject *)a)->content_data, PyFloat_AS_DOUBLE(b),
                            ((PyDDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(b, &PyQDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub_dd_qd(((PyDDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(b)||PyLong_Check(b)) { /* Int/Long */
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyDDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_sub(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            Py_DECREF(b);
            return o;
        }
    } else {
        if(PyFloat_Check(a)) {
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_sub_d_dd(PyFloat_AS_DOUBLE(a), ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub_qd_dd(((PyQDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(a)||PyLong_Check(a)) { /* Int/Long */
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyDDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_sub(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            Py_DECREF(a);
            return o;
        }
    }

    return Py_NotImplemented;
}

static PyObject *DD_mul(PyObject *a, PyObject *b) {
    PyDDTypeObject *left, *right;
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
        left = (PyDDTypeObject *) a;
        right = (PyDDTypeObject *) b;
    } else {
        left = (PyDDTypeObject *) b;
        right = (PyDDTypeObject *) a;
    }

    if(PyObject_TypeCheck(right, &PyDDTypeObjectType)) {
        o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
        c_dd_mul(left->content_data, right->content_data,
                        ((PyDDTypeObject *) o)->content_data);
        return o;
    }

    if(PyFloat_Check(right)) {
        o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
        c_dd_mul_dd_d(left->content_data, PyFloat_AS_DOUBLE(right),
                        ((PyDDTypeObject *) o)->content_data);
        return o;
    }

    if(PyObject_TypeCheck(right, &PyQDTypeObjectType)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_mul_dd_qd(left->content_data, ((PyQDTypeObject *)right)->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyInt_Check(right)||PyLong_Check(right)) { /* Int/Long */
        o = PyTuple_Pack(1, right);
        right = (PyDDTypeObject *) PyObject_CallObject( 
                           (PyObject *)&PyDDTypeObjectType, o);
        Py_DECREF(o);

        o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
        c_dd_mul(left->content_data, right->content_data,
                        ((PyDDTypeObject *) o)->content_data);
        Py_DECREF(right);
        return o;
    }

    return Py_NotImplemented;
}

static PyObject *DD_div(PyObject *a, PyObject *b) {
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
        if(PyObject_TypeCheck(b, &PyDDTypeObjectType)) {
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_div(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            return o;
        }

        if(PyFloat_Check(b)) {
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_div_dd_d(((PyDDTypeObject *)a)->content_data, PyFloat_AS_DOUBLE(b),
                            ((PyDDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(b, &PyQDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div_dd_qd(((PyDDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(b)||PyLong_Check(b)) { /* Int/Long */
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyDDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_div(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            Py_DECREF(b);
            return o;
        }
    } else {
        if(PyFloat_Check(a)) {
            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_div_d_dd(PyFloat_AS_DOUBLE(a), ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div_qd_dd(((PyQDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(a)||PyLong_Check(a)) { /* Int/Long */
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyDDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyDDTypeObjectType, NULL, NULL);
            c_dd_div(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyDDTypeObject *) o)->content_data);
            Py_DECREF(a);
            return o;
        }
    }

    return Py_NotImplemented;
}

static PyObject *DD_cmp(PyObject *a, PyObject *b, int op) {
    int result;
    PyObject *o;
    if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
        if(PyObject_TypeCheck(b, &PyDDTypeObjectType)) {
            c_dd_comp( ((PyDDTypeObject *) a)->content_data,
                       ((PyDDTypeObject *) b)->content_data, &result);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyFloat_Check(b)) {
            c_dd_comp_dd_d( ((PyDDTypeObject *) a)->content_data,
                       PyFloat_AS_DOUBLE(b), &result);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyObject_TypeCheck(b, &PyQDTypeObjectType)) {
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            c_qd_comp( ((PyQDTypeObject *) a)->content_data,
                       ((PyQDTypeObject *) b)->content_data, &result);
            Py_DECREF(a);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyInt_Check(b)||PyLong_Check(b)) {
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyDDTypeObjectType, o);
            Py_DECREF(o);

            c_dd_comp( ((PyDDTypeObject *) a)->content_data,
                       ((PyDDTypeObject *) b)->content_data, &result);
            Py_DECREF(b);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
    } else {
        if(PyFloat_Check(a)) {
            c_dd_comp_d_dd( PyFloat_AS_DOUBLE(a),
                            ((PyDDTypeObject *) b)->content_data, &result);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            c_qd_comp( ((PyQDTypeObject *) a)->content_data,
                       ((PyQDTypeObject *) b)->content_data, &result);
            Py_DECREF(b);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyInt_Check(a)||PyLong_Check(a)) {
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyDDTypeObjectType, o);
            Py_DECREF(o);

            c_dd_comp( ((PyDDTypeObject *) a)->content_data,
                       ((PyDDTypeObject *) b)->content_data, &result);
            Py_DECREF(a);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
    }
    return Py_NotImplemented;
}

static PyMethodDef DD_methods[] = {
    {"__reset__",(PyCFunction) DD_reset, METH_VARARGS,"Reset the value of a DD instance."},
    {"setRandom",(PyCFunction) DD_set_random, METH_NOARGS,"Set the DD instance to a random value."},
    {"clone",(PyCFunction) DD_clone, METH_NOARGS,"Clone the DD instance."},
    {"getAddress",(PyCFunction) DD_get_address, METH_NOARGS,"Get the address of the internal components of the DD instance."},
    {"getContent",(PyCFunction) DD_get_content, METH_NOARGS,"Get the internal components of the DD instance."},
    {"sqrt",(PyCFunction) DD_sqrt, METH_NOARGS,"Compute sqrt(x) for a DD instance."},
    {"sqr",(PyCFunction) DD_sqr, METH_NOARGS,"Compute x**2 for a DD instance."},
    {"nint",(PyCFunction) DD_nint, METH_NOARGS,"Compute nint(x) for a DD instance."},
    {"aint",(PyCFunction) DD_aint, METH_NOARGS,"Compute aint(x) for a DD instance."},
    {"floor",(PyCFunction) DD_floor, METH_NOARGS,"Compute floor(x) for a DD instance."},
    {"ceil",(PyCFunction) DD_ceil, METH_NOARGS,"Compute ceil(x) for a DD instance."},
    {"exp",(PyCFunction) DD_exp, METH_NOARGS,"Compute exp(x) for a DD instance."},
    {"log",(PyCFunction) DD_log, METH_NOARGS,"Compute log(x) for a DD instance."},
    {"log10",(PyCFunction) DD_log10, METH_NOARGS,"Compute log10(x) for a DD instance."},
    {"sin",(PyCFunction) DD_sin, METH_NOARGS,"Compute sin(x) for a DD instance."},
    {"cos",(PyCFunction) DD_cos, METH_NOARGS,"Compute cos(x) for a DD instance."},
    {"tan",(PyCFunction) DD_tan, METH_NOARGS,"Compute tan(x) for a DD instance."},
    {"asin",(PyCFunction) DD_asin, METH_NOARGS,"Compute asin(x) for a DD instance."},
    {"acos",(PyCFunction) DD_acos, METH_NOARGS,"Compute acos(x) for a DD instance."},
    {"atan",(PyCFunction) DD_atan, METH_NOARGS,"Compute atan(x) for a DD instance."},
    {"sinh",(PyCFunction) DD_sinh, METH_NOARGS,"Compute sinh(x) for a DD instance."},
    {"cosh",(PyCFunction) DD_cosh, METH_NOARGS,"Compute cosh(x) for a DD instance."},
    {"tanh",(PyCFunction) DD_tanh, METH_NOARGS,"Compute tanh(x) for a DD instance."},
    {"asinh",(PyCFunction) DD_asinh, METH_NOARGS,"Compute asinh(x) for a DD instance."},
    {"acosh",(PyCFunction) DD_acosh, METH_NOARGS,"Compute acosh(x) for a DD instance."},
    {"atanh",(PyCFunction) DD_atanh, METH_NOARGS,"Compute atanh(x) for a DD instance."},
    {"sincos",(PyCFunction) DD_sincos, METH_NOARGS,"Compute (sin(x), cos(x)) for a DD instance."},
    {"sincosh",(PyCFunction) DD_sincosh, METH_NOARGS,"Compute (sinh(x), cosh(x)) for a DD instance."},
    {"npwr",(PyCFunction) DD_npwr, METH_VARARGS,"Compute x**n for a DD instance."},
    {"nroot",(PyCFunction) DD_nroot, METH_VARARGS,"Compute x**(1/n) for a DD instance."},
    {"atan2",(PyCFunction) DD_atan2, METH_STATIC,"Static function for computing atan2(x,y)."},
    {NULL}
};
static PyNumberMethods DD_number_methods = {
    (binaryfunc) DD_add, /*    binaryfunc nb_add; */
    (binaryfunc) DD_sub, /*    binaryfunc nb_subtract; */
    (binaryfunc) DD_mul, /*    binaryfunc nb_multiply; */
    (binaryfunc) DD_div, /*    binaryfunc nb_divide; */
    0, /*    binaryfunc nb_remainder; */
    0, /*    binaryfunc nb_divmod; */
    (ternaryfunc) DD_pow, /*    ternaryfunc nb_power; */
    (unaryfunc) DD_neg, /*    unaryfunc nb_negative; */
    (unaryfunc) DD_clone, /*    unaryfunc nb_positive; */
    (unaryfunc) DD_abs, /*    unaryfunc nb_absolute; */
    0, /*    inquiry nb_nonzero; */
    0, /*    unaryfunc nb_invert; */
    0, /*    binaryfunc nb_lshift; */
    0, /*    binaryfunc nb_rshift; */
    0, /*    binaryfunc nb_and; */
    0, /*    binaryfunc nb_xor; */
    0, /*    binaryfunc nb_or; */
    0, /*    coercion nb_coerce; */
    0, /*    unaryfunc nb_int; */
    0, /*    unaryfunc nb_long; */
    (unaryfunc) DD_float, /*    unaryfunc nb_float; */
    0, /*    unaryfunc nb_oct; */
    0, /*    unaryfunc nb_hex; */
    /* Added in release 2.0 */
    0, /*    binaryfunc nb_inplace_add; */
    0, /*    binaryfunc nb_inplace_subtract; */
    0, /*    binaryfunc nb_inplace_multiply; */
    0, /*    binaryfunc nb_inplace_divide; */
    0, /*    binaryfunc nb_inplace_remainder; */
    0, /*    ternaryfunc nb_inplace_power; */
    0, /*    binaryfunc nb_inplace_lshift; */
    0, /*    binaryfunc nb_inplace_rshift; */
    0, /*    binaryfunc nb_inplace_and; */
    0, /*    binaryfunc nb_inplace_xor; */
    0, /*    binaryfunc nb_inplace_or; */

    /* Added in release 2.2 */
    /* The following require the Py_TPFLAGS_HAVE_CLASS flag */
    0, /*    binaryfunc nb_floor_divide; */
    0, /*    binaryfunc nb_true_divide; */
    0, /*    binaryfunc nb_inplace_floor_divide; */
    0, /*    binaryfunc nb_inplace_true_divide; */

    /* Added in release 2.5 */
    0, /*    unaryfunc nb_index; */
};

static PyTypeObject PyDDTypeObjectType = {
        PyObject_HEAD_INIT(NULL)
        0,                              /* ob_size        */
        "qd.DD",                        /* tp_name        */
        sizeof(PyDDTypeObject),         /* tp_basicsize   */
        0,                              /* tp_itemsize    */
        0,                              /* tp_dealloc     */
        0,                              /* tp_print       */
        0,                              /* tp_getattr     */
        0,                              /* tp_setattr     */
        0,                              /* tp_compare     */
        (reprfunc)DD_repr,              /* tp_repr        */
        &DD_number_methods,             /* tp_as_number   */
        0,                              /* tp_as_sequence */
        0,                              /* tp_as_mapping  */
        0,                              /* tp_hash        */
        0,                              /* tp_call        */
        0,                              /* tp_str         */
        0,                              /* tp_getattro    */
        0,                              /* tp_setattro    */
        0,                              /* tp_as_buffer   */
        Py_TPFLAGS_DEFAULT|Py_TPFLAGS_CHECKTYPES,             /* tp_flags       */
        "Wrapper for the quad double type.",   /* tp_doc         */
        0,                              /* tp_traverse       */
        0,                              /* tp_clear          */
        (richcmpfunc)DD_cmp,            /* tp_richcompare    */
        0,                              /* tp_weaklistoffset */
        0,                              /* tp_iter           */
        0,                              /* tp_iternext       */
        DD_methods,                     /* tp_methods        */
        0,                              /* tp_members        */
        0,                              /* tp_getset         */
        0,                              /* tp_base           */
        0,                              /* tp_dict           */
        0,                              /* tp_descr_get      */
        0,                              /* tp_descr_set      */
        0,                              /* tp_dictoffset     */
        (initproc)DD_init,              /* tp_init           */
};

#define QD_GENERICWRAPPER1(f,self) \
            PyObject *o; \
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL); \
            f(self->content_data, ((PyQDTypeObject *) o)->content_data); \
            return o;

#define QD_GENERICWRAPPER2(f,self) \
            PyObject *a,*b; \
            a = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL); \
            b = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL); \
            f(self->content_data, ((PyQDTypeObject *) a)->content_data, \
                                  ((PyQDTypeObject *) b)->content_data); \
            return PyTuple_Pack(2, a,b);
#define QD_GENERICWRAPPER3(f,self, args) \
            PyObject *o; \
            int n; \
            if (!PyArg_ParseTuple(args, "i", &n)) return NULL; \
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL); \
            f(self->content_data, n, ((PyQDTypeObject *) o)->content_data); \
            return o;

static int QD_init(PyQDTypeObject *self, PyObject *args) {
    PyObject *a = NULL;
    PyObject *tmp;

    if (!PyArg_ParseTuple(args, "|O", &a)) return -1;

    if(a==NULL) return 0; /* 0 argument */

    if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) { /* qd.QD */
        self->content_data[0] = ((PyQDTypeObject *)a)->content_data[0];
        self->content_data[1] = ((PyQDTypeObject *)a)->content_data[1];
        self->content_data[2] = ((PyQDTypeObject *)a)->content_data[2];
        self->content_data[3] = ((PyQDTypeObject *)a)->content_data[3];
        return 0;
    }

    if(PyFloat_Check(a)) { /* Float */
        self->content_data[0] = PyFloat_AS_DOUBLE(a);
        self->content_data[1] = 0.0;
        self->content_data[2] = 0.0;
        self->content_data[3] = 0.0;
        return 0;
    }

    if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) { /* qd.DD */
        self->content_data[0] = ((PyDDTypeObject *)a)->content_data[0];
        self->content_data[1] = ((PyDDTypeObject *)a)->content_data[1];
        self->content_data[2] = 0.0;
        self->content_data[3] = 0.0;
        return 0;
    }

    if(PyInt_Check(a)||PyLong_Check(a)) { /* Int/Long */
        tmp = PyObject_Str(a);
        c_qd_read(PyString_AS_STRING(tmp), self->content_data);
        Py_DECREF(tmp);
        return 0;
    }

    if(PyString_Check(a)) { /* String */
        c_qd_read(PyString_AS_STRING(a), self->content_data);
        if( isnan(self->content_data[0] )) {
            PyErr_SetString(PyExc_ValueError,"Converting string to QD gave a NaN Value");
            return -1;
        }
        return 0;
    }

    tmp = PyNumber_Float(a);
    if (tmp != NULL) { /* __float__ */
        self->content_data[0] = PyFloat_AS_DOUBLE(tmp);
        self->content_data[1] = 0.0;
        self->content_data[2] = 0.0;
        self->content_data[3] = 0.0;
        Py_DECREF(tmp);
        return 0;
    }

    return -1;
}

static PyObject *QD_reset(PyQDTypeObject *self, PyObject *args) {
    (void) QD_init(self, args);
    Py_RETURN_NONE;
}

static PyObject *QD_repr(PyQDTypeObject *self) {
    char repr[72];
    c_qd_swrite( self->content_data, 64, repr, 72 );
    return PyString_FromString(repr);
}

static PyObject *QD_sqrt(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_sqrt, self) }
static PyObject *QD_sqr(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_sqr, self) }
static PyObject *QD_abs(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_abs, self) }
static PyObject *QD_neg(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_neg, self) }
static PyObject *QD_nint(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_nint, self) }
static PyObject *QD_aint(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_aint, self) }
static PyObject *QD_floor(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_floor, self) }
static PyObject *QD_ceil(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_ceil, self) }
static PyObject *QD_exp(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_exp, self) }
static PyObject *QD_log(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_log, self) }
static PyObject *QD_log10(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_log10, self) }
static PyObject *QD_sin(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_sin, self) }
static PyObject *QD_cos(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_cos, self) }
static PyObject *QD_tan(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_tan, self) }
static PyObject *QD_asin(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_asin, self) }
static PyObject *QD_acos(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_acos, self) }
static PyObject *QD_atan(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_atan, self) }
static PyObject *QD_sinh(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_sinh, self) }
static PyObject *QD_cosh(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_cosh, self) }
static PyObject *QD_tanh(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_tanh, self) }
static PyObject *QD_asinh(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_asinh, self) }
static PyObject *QD_acosh(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_acosh, self) }
static PyObject *QD_atanh(PyQDTypeObject *self) { QD_GENERICWRAPPER1(c_qd_atanh, self) }
static PyObject *QD_sincos(PyQDTypeObject *self) { QD_GENERICWRAPPER2(c_qd_sincos, self) }
static PyObject *QD_sincosh(PyQDTypeObject *self) { QD_GENERICWRAPPER2(c_qd_sincosh, self) }
static PyObject *QD_npwr(PyQDTypeObject *self, PyObject *args) {
    QD_GENERICWRAPPER3(c_qd_npwr, self, args) }
static PyObject *QD_nroot(PyQDTypeObject *self, PyObject *args) {
    QD_GENERICWRAPPER3(c_qd_nroot, self, args) }

static PyObject *QD_atan2(PyQDTypeObject *self, PyObject *args) {
    PyObject *x;
    PyObject *y;
    PyObject *o;

    if (!PyArg_ParseTuple(args, "OO", &x, &y)) return NULL;

    if(!PyObject_TypeCheck(x, &PyQDTypeObjectType)) {
        PyErr_SetString(PyExc_TypeError,"Wrong type for argument 1 of qd.QD.atan2()");
        return NULL;
    }
    if(!PyObject_TypeCheck(y, &PyQDTypeObjectType)) {
        PyErr_SetString(PyExc_TypeError,"Wrong type for argument 2 of qd.QD.atan2()");
        return NULL;
    }

    o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
    c_qd_atan2( ((PyQDTypeObject *) x)->content_data,
                ((PyQDTypeObject *) y)->content_data,
                ((PyQDTypeObject *) o)->content_data );

    return o;
}

static PyObject *QD_float(PyQDTypeObject *self) {
    return PyFloat_FromDouble( self->content_data[0] );
}

static PyObject *QD_set_random(PyQDTypeObject *self) {
    c_qd_rand( self->content_data );
    Py_RETURN_NONE;
}

static PyObject *QD_clone(PyQDTypeObject *self) {
    PyObject *o;
    o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
    ((PyQDTypeObject *) o)->content_data[0] = self->content_data[0];
    ((PyQDTypeObject *) o)->content_data[1] = self->content_data[1];
    ((PyQDTypeObject *) o)->content_data[2] = self->content_data[2];
    ((PyQDTypeObject *) o)->content_data[3] = self->content_data[3];
    return o;
}

static PyObject *QD_get_address(PyQDTypeObject *self) {
    return PyInt_FromLong( (long) self->content_data );
}

static PyObject *QD_get_content(PyQDTypeObject *self) {
    return PyTuple_Pack(4, 
           PyFloat_FromDouble( self->content_data[0] ),
           PyFloat_FromDouble( self->content_data[1] ),
           PyFloat_FromDouble( self->content_data[2] ),
           PyFloat_FromDouble( self->content_data[3] ) );
}

static PyObject *QD_pow(PyQDTypeObject *self, PyObject *a, PyObject *b) {
    PyObject *o;

    if (b != Py_None) {
        PyErr_SetString(PyExc_NotImplementedError,"Not implemented operation");
        return NULL;
    }

    if((!PyInt_Check(a))&&(!PyLong_Check(a))) {
        PyErr_SetString(PyExc_TypeError,"Wrong type for exponent");
        return NULL;
    }

    o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
    c_qd_npwr(self->content_data, PyInt_AS_LONG(a),
            ((PyQDTypeObject *) o)->content_data);
    return o;
}

static PyObject *QD_add(PyObject *a, PyObject *b) {
    PyQDTypeObject *left, *right;
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
        left = (PyQDTypeObject *) a;
        right = (PyQDTypeObject *) b;
    } else {
        left = (PyQDTypeObject *) b;
        right = (PyQDTypeObject *) a;
    }

    if(PyObject_TypeCheck(right, &PyQDTypeObjectType)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_add(left->content_data, right->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyFloat_Check(right)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_add_qd_d(left->content_data, PyFloat_AS_DOUBLE(right),
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyObject_TypeCheck(right, &PyDDTypeObjectType)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_add_qd_dd(left->content_data, right->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyInt_Check(right)||PyLong_Check(right)) { /* Int/Long */
        o = PyTuple_Pack(1, right);
        right = (PyQDTypeObject *) PyObject_CallObject( 
                           (PyObject *)&PyQDTypeObjectType, o);
        Py_DECREF(o);

        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_add(left->content_data, right->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        Py_DECREF(right);
        return o;
    }

    return Py_NotImplemented;
}

static PyObject *QD_sub(PyObject *a, PyObject *b) {
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
        if(PyObject_TypeCheck(b, &PyQDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub(((PyQDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyFloat_Check(b)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub_qd_d(((PyQDTypeObject *)a)->content_data, PyFloat_AS_DOUBLE(b),
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(b, &PyDDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub_qd_dd(((PyQDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(b)||PyLong_Check(b)) { /* Int/Long */
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub(((PyQDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            Py_DECREF(b);
            return o;
        }
    } else {
        if(PyFloat_Check(a)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub_d_qd(PyFloat_AS_DOUBLE(a), ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub_dd_qd(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(a)||PyLong_Check(a)) { /* Int/Long */
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_sub(((PyQDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            Py_DECREF(a);
            return o;
        }
    }

    return Py_NotImplemented;
}

static PyObject *QD_mul(PyObject *a, PyObject *b) {
    PyQDTypeObject *left, *right;
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
        left = (PyQDTypeObject *) a;
        right = (PyQDTypeObject *) b;
    } else {
        left = (PyQDTypeObject *) b;
        right = (PyQDTypeObject *) a;
    }

    if(PyObject_TypeCheck(right, &PyQDTypeObjectType)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_mul(left->content_data, right->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyFloat_Check(right)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_mul_qd_d(left->content_data, PyFloat_AS_DOUBLE(right),
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyObject_TypeCheck(right, &PyDDTypeObjectType)) {
        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_mul_qd_dd(left->content_data, right->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        return o;
    }

    if(PyInt_Check(right)||PyLong_Check(right)) { /* Int/Long */
        o = PyTuple_Pack(1, right);
        right = (PyQDTypeObject *) PyObject_CallObject( 
                           (PyObject *)&PyQDTypeObjectType, o);
        Py_DECREF(o);

        o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
        c_qd_mul(left->content_data, right->content_data,
                        ((PyQDTypeObject *) o)->content_data);
        Py_DECREF(right);
        return o;
    }

    return Py_NotImplemented;
}

static PyObject *QD_div(PyObject *a, PyObject *b) {
    PyObject *o;

    if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
        if(PyObject_TypeCheck(b, &PyQDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div(((PyQDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyFloat_Check(b)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div_qd_d(((PyQDTypeObject *)a)->content_data, PyFloat_AS_DOUBLE(b),
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(b, &PyDDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div_qd_dd(((PyQDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(b)||PyLong_Check(b)) { /* Int/Long */
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div(((PyQDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            Py_DECREF(b);
            return o;
        }
    } else {
        if(PyFloat_Check(a)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div_d_qd(PyFloat_AS_DOUBLE(a), ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div_dd_qd(((PyDDTypeObject *)a)->content_data,
                     ((PyDDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            return o;
        }

        if(PyInt_Check(a)||PyLong_Check(a)) { /* Int/Long */
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            o = PyType_GenericNew( &PyQDTypeObjectType, NULL, NULL);
            c_qd_div(((PyQDTypeObject *)a)->content_data,
                     ((PyQDTypeObject *)b)->content_data,
                            ((PyQDTypeObject *) o)->content_data);
            Py_DECREF(a);
            return o;
        }
    }

    return Py_NotImplemented;
}

static PyObject *QD_cmp(PyObject *a, PyObject *b, int op) {
    int result;
    PyObject *o;
    if(PyObject_TypeCheck(a, &PyQDTypeObjectType)) {
        if(PyObject_TypeCheck(b, &PyQDTypeObjectType)) {
            c_qd_comp( ((PyQDTypeObject *) a)->content_data,
                       ((PyQDTypeObject *) b)->content_data, &result);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyFloat_Check(b)) {
            c_qd_comp_qd_d( ((PyQDTypeObject *) a)->content_data,
                       PyFloat_AS_DOUBLE(b), &result);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyObject_TypeCheck(b, &PyDDTypeObjectType)) {
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            c_qd_comp( ((PyQDTypeObject *) a)->content_data,
                       ((PyQDTypeObject *) b)->content_data, &result);
            Py_DECREF(b);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyInt_Check(b)||PyLong_Check(b)) {
            o = PyTuple_Pack(1, b);
            b = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            c_qd_comp( ((PyQDTypeObject *) a)->content_data,
                       ((PyQDTypeObject *) b)->content_data, &result);
            Py_DECREF(b);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
    } else {
        if(PyFloat_Check(a)) {
            c_qd_comp_d_qd( PyFloat_AS_DOUBLE(a),
                            ((PyQDTypeObject *) b)->content_data, &result);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyObject_TypeCheck(a, &PyDDTypeObjectType)) {
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            c_qd_comp( ((PyQDTypeObject *) a)->content_data,
                       ((PyQDTypeObject *) b)->content_data, &result);
            Py_DECREF(a);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
        if(PyInt_Check(a)||PyLong_Check(a)) {
            o = PyTuple_Pack(1, a);
            a = PyObject_CallObject( (PyObject *)&PyQDTypeObjectType, o);
            Py_DECREF(o);

            c_qd_comp( ((PyQDTypeObject *) a)->content_data,
                       ((PyQDTypeObject *) b)->content_data, &result);
            Py_DECREF(a);
            return PyBool_FromLong(
                (op==0)?result==-1:(
                    (op==2)?result==0:(
                        (op==4)?result==1:(
                            (op==1)?result<=0:(
                                (op==3)?result!=0:result>=0)))));
        }
    }
    return Py_NotImplemented;
}

static PyMethodDef QD_methods[] = {
    {"__reset__",(PyCFunction) QD_reset, METH_VARARGS,"Reset the value of a QD instance."},
    {"setRandom",(PyCFunction) QD_set_random, METH_NOARGS,"Set the QD instance to a random value."},
    {"clone",(PyCFunction) QD_clone, METH_NOARGS,"Clone the QD instance."},
    {"getAddress",(PyCFunction) QD_get_address, METH_NOARGS,"Get the address of the internal components of the QD instance."},
    {"getContent",(PyCFunction) QD_get_content, METH_NOARGS,"Get the internal components of the QD instance."},
    {"sqrt",(PyCFunction) QD_sqrt, METH_NOARGS,"Compute sqrt(x) for a QD instance."},
    {"sqr",(PyCFunction) QD_sqr, METH_NOARGS,"Compute x**2 for a QD instance."},
    {"nint",(PyCFunction) QD_nint, METH_NOARGS,"Compute nint(x) for a QD instance."},
    {"aint",(PyCFunction) QD_aint, METH_NOARGS,"Compute aint(x) for a QD instance."},
    {"floor",(PyCFunction) QD_floor, METH_NOARGS,"Compute floor(x) for a QD instance."},
    {"ceil",(PyCFunction) QD_ceil, METH_NOARGS,"Compute ceil(x) for a QD instance."},
    {"exp",(PyCFunction) QD_exp, METH_NOARGS,"Compute exp(x) for a QD instance."},
    {"log",(PyCFunction) QD_log, METH_NOARGS,"Compute log(x) for a QD instance."},
    {"log10",(PyCFunction) QD_log10, METH_NOARGS,"Compute log10(x) for a QD instance."},
    {"sin",(PyCFunction) QD_sin, METH_NOARGS,"Compute sin(x) for a QD instance."},
    {"cos",(PyCFunction) QD_cos, METH_NOARGS,"Compute cos(x) for a QD instance."},
    {"tan",(PyCFunction) QD_tan, METH_NOARGS,"Compute tan(x) for a QD instance."},
    {"asin",(PyCFunction) QD_asin, METH_NOARGS,"Compute asin(x) for a QD instance."},
    {"acos",(PyCFunction) QD_acos, METH_NOARGS,"Compute acos(x) for a QD instance."},
    {"atan",(PyCFunction) QD_atan, METH_NOARGS,"Compute atan(x) for a QD instance."},
    {"sinh",(PyCFunction) QD_sinh, METH_NOARGS,"Compute sinh(x) for a QD instance."},
    {"cosh",(PyCFunction) QD_cosh, METH_NOARGS,"Compute cosh(x) for a QD instance."},
    {"tanh",(PyCFunction) QD_tanh, METH_NOARGS,"Compute tanh(x) for a QD instance."},
    {"asinh",(PyCFunction) QD_asinh, METH_NOARGS,"Compute asinh(x) for a QD instance."},
    {"acosh",(PyCFunction) QD_acosh, METH_NOARGS,"Compute acosh(x) for a QD instance."},
    {"atanh",(PyCFunction) QD_atanh, METH_NOARGS,"Compute atanh(x) for a QD instance."},
    {"sincos",(PyCFunction) QD_sincos, METH_NOARGS,"Compute (sin(x), cos(x)) for a QD instance."},
    {"sincosh",(PyCFunction) QD_sincosh, METH_NOARGS,"Compute (sinh(x), cosh(x)) for a QD instance."},
    {"npwr",(PyCFunction) QD_npwr, METH_VARARGS,"Compute x**n for a QD instance."},
    {"nroot",(PyCFunction) QD_nroot, METH_VARARGS,"Compute x**(1/n) for a QD instance."},
    {"atan2",(PyCFunction) QD_atan2, METH_STATIC,"Static function for computing atan2(x,y)."},
    {NULL}
};
static PyNumberMethods QD_number_methods = {
    (binaryfunc) QD_add, /*    binaryfunc nb_add; */
    (binaryfunc) QD_sub, /*    binaryfunc nb_subtract; */
    (binaryfunc) QD_mul, /*    binaryfunc nb_multiply; */
    (binaryfunc) QD_div, /*    binaryfunc nb_divide; */
    0, /*    binaryfunc nb_remainder; */
    0, /*    binaryfunc nb_divmod; */
    (ternaryfunc) QD_pow, /*    ternaryfunc nb_power; */
    (unaryfunc) QD_neg, /*    unaryfunc nb_negative; */
    (unaryfunc) QD_clone, /*    unaryfunc nb_positive; */
    (unaryfunc) QD_abs, /*    unaryfunc nb_absolute; */
    0, /*    inquiry nb_nonzero; */
    0, /*    unaryfunc nb_invert; */
    0, /*    binaryfunc nb_lshift; */
    0, /*    binaryfunc nb_rshift; */
    0, /*    binaryfunc nb_and; */
    0, /*    binaryfunc nb_xor; */
    0, /*    binaryfunc nb_or; */
    0, /*    coercion nb_coerce; */
    0, /*    unaryfunc nb_int; */
    0, /*    unaryfunc nb_long; */
    (unaryfunc) QD_float, /*    unaryfunc nb_float; */
    0, /*    unaryfunc nb_oct; */
    0, /*    unaryfunc nb_hex; */
    /* Added in release 2.0 */
    0, /*    binaryfunc nb_inplace_add; */
    0, /*    binaryfunc nb_inplace_subtract; */
    0, /*    binaryfunc nb_inplace_multiply; */
    0, /*    binaryfunc nb_inplace_divide; */
    0, /*    binaryfunc nb_inplace_remainder; */
    0, /*    ternaryfunc nb_inplace_power; */
    0, /*    binaryfunc nb_inplace_lshift; */
    0, /*    binaryfunc nb_inplace_rshift; */
    0, /*    binaryfunc nb_inplace_and; */
    0, /*    binaryfunc nb_inplace_xor; */
    0, /*    binaryfunc nb_inplace_or; */

    /* Added in release 2.2 */
    /* The following require the Py_TPFLAGS_HAVE_CLASS flag */
    0, /*    binaryfunc nb_floor_divide; */
    0, /*    binaryfunc nb_true_divide; */
    0, /*    binaryfunc nb_inplace_floor_divide; */
    0, /*    binaryfunc nb_inplace_true_divide; */

    /* Added in release 2.5 */
    0, /*    unaryfunc nb_index; */
};

static PyTypeObject PyQDTypeObjectType = {
        PyObject_HEAD_INIT(NULL)
        0,                              /* ob_size        */
        "qd.QD",                        /* tp_name        */
        sizeof(PyQDTypeObject),         /* tp_basicsize   */
        0,                              /* tp_itemsize    */
        0,                              /* tp_dealloc     */
        0,                              /* tp_print       */
        0,                              /* tp_getattr     */
        0,                              /* tp_setattr     */
        0,                              /* tp_compare     */
        (reprfunc)QD_repr,              /* tp_repr        */
        &QD_number_methods,             /* tp_as_number   */
        0,                              /* tp_as_sequence */
        0,                              /* tp_as_mapping  */
        0,                              /* tp_hash        */
        0,                              /* tp_call        */
        0,                              /* tp_str         */
        0,                              /* tp_getattro    */
        0,                              /* tp_setattro    */
        0,                              /* tp_as_buffer   */
        Py_TPFLAGS_DEFAULT|Py_TPFLAGS_CHECKTYPES,             /* tp_flags       */
        "Wrapper for the quad double type.",   /* tp_doc         */
        0,                              /* tp_traverse       */
        0,                              /* tp_clear          */
        (richcmpfunc)QD_cmp,            /* tp_richcompare    */
        0,                              /* tp_weaklistoffset */
        0,                              /* tp_iter           */
        0,                              /* tp_iternext       */
        QD_methods,                     /* tp_methods        */
        0,                              /* tp_members        */
        0,                              /* tp_getset         */
        0,                              /* tp_base           */
        0,                              /* tp_dict           */
        0,                              /* tp_descr_get      */
        0,                              /* tp_descr_set      */
        0,                              /* tp_dictoffset     */
        (initproc)QD_init,              /* tp_init           */
};

static PyMethodDef functions[] = {
        {"fpu_init", fpu_init, METH_NOARGS,"Initialize the FPU state for using the QD library. The function should be called before computing with DD or QD types. The initial state of the FPU can be reset later with the function fpu_restore()."},
        {"fpu_restore", fpu_restore, METH_NOARGS,"Restore the initial FPU state."},
        { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initqd(void)
{
    PyObject* m;

    PyQDTypeObjectType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PyQDTypeObjectType) < 0)
            return;
    PyDDTypeObjectType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PyDDTypeObjectType) < 0)
            return;

    m = Py_InitModule3("qd", functions,
                       "Interface to the libqd library.");
    if (m == NULL)
            return;

    Py_INCREF(&PyQDTypeObjectType);
    PyModule_AddObject(m, "QD", (PyObject *)&PyQDTypeObjectType);
    PyModule_AddObject(m, "DD", (PyObject *)&PyDDTypeObjectType);
}
