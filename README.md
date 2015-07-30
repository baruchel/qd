# qd
A python wrapper for the high precision QD library

### The QD module

The QD module for Python is a wrapper for the high-precision old and famous [QD](http://crd.lbl.gov/~dhbailey/mpdist/) library providing new types for ~32 or ~64 decimal digits of precision. Unlike other high-precision libraries, the idea was to store the number as a sum of two or four double numbers and to use relevant combination of arithmetic operations for computing with such numbers. This system provides a very high speed for computing many mathematical operations or functions.

The wrapper has been written in C with the API provided with CPython in order to be as quick as possible. It was written and tried with Python 2.

A new version will allow later to use these new types ad Numpy dtypes.

### Compilation of the module

First, install the excellent QD library; both the binary package and the development files have no dependancy and they are very easy to install. A package should exist for your distribution. For instance on Debian-based distributions, just type:

    sudo apt-get install libqd0 libqd-dev

Then compile the module by copying the files in some directory and typing:

    python setup.py build_ext --inplace

The resulting binary file `qd.so` is the compiled module; you can install it in some relevant location or just put it in the same directory as a project for trying it first.

    If someone can provide some explanations for compiling it and using it on
    other platforms, please let me know and I will copy them below.

### Using the module

Using the module is very easy; first import it with:

    import qd

Before starting to do dome computation with the DD or QD types provided by the module, you should set the FPU in the needed state with:

    qd.fpu_init()

When your code doesn't require computing with the new types any longer, restore the state of the FPU with:

    qd.fpu_restore()

Now you can start computing with both types:

    a = qd.DD(3.0)
    b = a+2
    c = b.sqrt()
    c < 3.0
        ... etc. ...

You can mix the DD and QD types with many other numerical types but try to mix them with float numbers rather than integer numbers (computation will be slower with integer numbers). If you want to mix them with other high-precision modules, convert your numbers to string before converting them to DD or QD types (otherwise they will first be converted to float numbers and you will lose precision).
