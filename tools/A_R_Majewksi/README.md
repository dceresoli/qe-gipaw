# EFG extraction and analysis tool
This tool has been contributed by A. R. Majewski and it's available from [https://github.com/Altoidnerd/python-pwscf](https://github.com/Altoidnerd/python-pwscf).

Example of interactive usage:

	ceresoli@adamello:~/Codes/qe-gipaw/tools/A_R_Majewksi$ ipython
	Python 2.7.13 (default, Sep 26 2018, 18:42:22) 
	Type "copyright", "credits" or "license" for more information.
	
	IPython 5.8.0 -- An enhanced Interactive Python.
	?         -> Introduction and overview of IPython's features.
	%quickref -> Quick reference.
	help      -> Python's own help system.
	object?   -> Details about 'object', use 'object??' for extra details.
	
	In [1]: from efg import Efg
	
	In [2]: e = Efg('quartz-efg.out')
	
	In [3]: print e.Cqs
	[0.1876, 0.1876, 0.1876, 5.4829, 5.4828, 5.4829, 5.4829, 5.4828, 5.4829]
