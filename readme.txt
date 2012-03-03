ABOUT
=====

ErwinJr is an open source design and simulation program for quantum semiconductor devices including quantum cascade lasers. It is useful for modeling intersubband semiconductor devices.  The code base was written at the NASA Jet Propulsion Laboratory.

ErwinJr is a multi-platform application that runs on most desktop operating systems including Windows, Mac OS X, and Linux.  This is because it is written in Python and uses the Qt graphical user interface framework.


== WINDOWS INSTALLATION ==

1) Do a FULL install of pythonxy from http://www.pythonxy.com/

2) Under Start|Run type cmd

3) type the following at the command prompt
   cd \erwinjr
   gcc -c cFunctions.c
   gcc -shared -o cFunctions.dll cFunctions.o
   
4) open (double click) on the file erwinjr.pyw


== Python Dependences ==

PyQt4, PyQwt5, numpy, scipy, matplotlib
psyco is also useful

Built and tested on Python 2.6.