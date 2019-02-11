#!/usr/bin/env python3
#
#
#
#################################################################################
#										#
# Copyright (c) 2016 Allen Majewski (altoidnerd)				#
# Permission is hereby granted, free of charge, to any person obtaining a 	#
# copy of this software and associated documentation files (the "Software"),	#
# to deal in the Software without restriction, including without limitation	#
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 	#
# and/or sell copies of the Software, and to permit persons to whom the		#	
# Software is furnished to do so, subject to the following conditions:		#
#										#
# The above copyright notice and this permission notice shall be included 	#
# in all copies or substantial portions of the Software.			#
#										#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 	#
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 	#
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 	#
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 	#
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,	#
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE	#
# THE SOFTWARE.									#                                
#										#
#################################################################################

#   python-pwscf was designed for python3; 
#   run on > 2.7 at your own risk.						

import sys
import numpy as np
from nqr_parser import f32
# import f1 for spin 1; f32 is for spin 3/2 (default)

def filtr(pattern, array):
  """
  filtr(str "parrtern", list array): -> list
  
    returns a new list containing the elements
    in array for which "pattern" is matched.
  """
  return list( filter( lambda x: pattern in x, array ) )

def indices(pattern, array):
  """
  indices(str "pattern", array): -> list
  
    returns a new list containing the indices 
    of the elements in array that match "pattern."
  """
  return [ array.index(thing) for thing in filtr(pattern, array) ]

def lmap(func, array):
  """
  lmap(function, array): -> list

  just like python map, but returns a list instead of a
 """
  return list( map( func, array ) ) 

def dict_to_object(dict_item):
  d = dict_item
  class _:
    pass
  for key, value in d.items():
    setattr(_, key, value)
  return _

asobject = dict_to_object


def is_right_handed(x, y, z):
  """
   determine if a set of axes
   is a right handed coordinate 
   system using the cross product
  """
  testz=np.cross(x, y)
  return np.all(np.isclose(z, testz, rtol=0.01))
    


class Efg(object):
  """
  Class for parsing gipaw.x output files
  for 'efg' calculation
  """
  
  def __init__(self, efgfile):
    self.efgfile = efgfile
    self.infile = efgfile  
    self.efgfile_array = open(efgfile, 'r').readlines()
    self.file_array = self.efgfile_array
    self.labels = self.atom_labels
    self.specie = set(self.atomic_species)
  #  self.data = "Attribute not loaded.  Use self.get_data() to assign self.data."  

  def __repr__(self):
    return "< efg object; efgfile: {}\tlength: {} lines >".format(self.efgfile, len(self.efgfile_array))

  def __str__(self):
    s = ""
    for line in self.file_array:
      s += line
    return s

  def prnt(self):
    for line in self.file_array:
      sys.stdout.write(line)


  @property
  def mdstep(self):
    def getnum(string):
      _ = ""
      for i in string:
        if i in '0123456789':
          _ += i
      return int(_)
    return getnum(self.infile)


  @property
  def nat(self):
    return len([ line for line in self.file_array if 'Q=' in line ])

    
  #@property
  #def WALL_TIME(self):
  #  return float([ line.split()[4] for line in filtr('WALL',self.file_array) ][-1].replace('s',''))

    
  #@property
  #def CPU_TIME(self):
  #  return float([ line.split()[2] for line in filtr('WALL',self.file_array) ][-1].replace('s',''))


  @property
  def atom_labels(self):
    a =  filtr("Q=", self.file_array)
    b =  [ line.strip()[:6] for line in a ]
    return b


  @property
  def atomic_species(self):
    a =  filtr("Q=",self.file_array)
    b =  [ line.split()[0] for line in a ]
    return b


  @property
  def total_efg(self, atom=None):
    begin_on_pattern= '----- total EFG -----'
    begin_on_index  = self.file_array.index(filtr(begin_on_pattern, self.file_array)[0])
    end_on_pattern  = '(symmetrized)'
    end_on_index    = self.file_array.index(filtr(end_on_pattern, self.file_array)[0])
    begin, end      = begin_on_index, end_on_index 
    file_slice      = self.file_array[begin+1:end]
    tensors         = [None]*self.nat
    labels          = self.atom_labels
    for i in range(len(labels)):
      tens = [ lmap( float, thing.split()[2:]) for thing in filtr(self.atom_labels[i], file_slice) ] 
      tensors[i] = tens
    return tensors
    

#  @property
#  def symmetrized_efg(self, atom=None):
#    begin_on_pattern= '(symmetrized)'
#    begin_on_index  = self.file_array.index(filtr(begin_on_pattern, self.file_array)[0])
#    end_on_pattern  = 'NQR/NMR SPECTROSCOPIC PARAMETERS'
#    end_on_index    = self.file_array.index(filtr(end_on_pattern, self.file_array)[0])
#    begin, end      = begin_on_index, end_on_index 
#    file_slice      = self.file_array[begin+1:end]
#    tensors         = [None]*self.nat
#    labels          = self.atom_labels
#    for i in range(len(labels)):
#      tens = [ lmap( float, thing.split()[2:]) for thing in filtr(self.atom_labels[i], file_slice) ] 
#      tensors[i] = tens
#    return tensors
#    print('Not implemented.')
#    pass

  #@property
  #def compd_eigvals(self, atom=None, sym=True, herm=False):
  #  if sym:
  #    efgs = self.symmetrized_efg
  #  else:
  #    efgs = self.total_efg
  #  if not herm:
  #    eigfunc = np.linalg.eig
  #  else:
  #    eigfunc = np.linalg.eigh
  #  return [ eigfunc(thing) for thing  in  efgs ]
  #  print('Not implemented.')
  #   pass


  @property
  def principle_axes(self):
    axes=[]
    for i in range(self.nat):
      label = self.atom_labels[i]
      vijs = filtr(label, self.file_array)[-4:][:3]
      Xaxis  = np.array( lmap(float, filtr('Vxx', vijs)[0].replace(')','').replace('(','').split()[5:]))
      Yaxis  = np.array( lmap(float, filtr('Vyy', vijs)[0].replace(')','').replace('(','').split()[5:]))
      Zaxis  = np.array( lmap(float, filtr('Vzz', vijs)[0].replace(')','').replace('(','').split()[5:]))
      if not(is_right_handed(Xaxis,Yaxis,Zaxis)):
        Xaxis,Yaxis,Zaxis = -Xaxis,-Yaxis,-Zaxis
      ax = []
      ax.append(Xaxis)
      ax.append(Yaxis)
      ax.append(Zaxis)
      axes.append(ax)
    return axes


  @property
  def Vii(self):
    viis = []
    for i in range(self.nat):
      label = self.atom_labels[i]
      vijs = filtr(label, self.file_array)[-4:][:3]
      vxx  = float(filtr('Vxx', vijs)[0].split()[3])
      vyy  = float(filtr('Vyy', vijs)[0].split()[3])
      vzz  = float(filtr('Vzz', vijs)[0].split()[3])
      vii  = [vxx,vyy,vzz]
      viis.append(vii)
    return viis
      
      
  @property
  def Qs(self, indices=None):
    a =  filtr("Q=",self.file_array)
    b =  [ line.split()[3] for line in a ]
    return lmap(float, b)
    

  @property
  def Cqs(self, indices=None):
    a =  filtr("Cq=",self.file_array)
    b =  [ line.split()[7] for line in a ]
    return lmap(float, b)


  @property
  def etas(self, indices=None):
    a =  filtr("eta=",self.file_array)
    b =  [ line.split()[10] for line in a ]
    return lmap(float, b)

  def tidy_zip(self, atom=None, verbose=False):
    help_message="""
		self.
		0	atom_labels
		1	Cqs
		2	etas
		3	Vii
	`	4	principle_axes
		5	Qs
		6	total_efg
		7	symmetrized_efg	
                 """

    if verbose:
      print(help_message)
   
    all_atoms = list( zip( 
			self.atom_labels, 
			self.Cqs, 
			self.etas, 
			self.Vii, 
			self.principle_axes, 
			self.Qs, 
			self.total_efg, 
#			self.symmetrized_efg
			)
		)
    if atom is None:
      return all_atoms
    else:
      return all_atoms[atom-1]

  
  def atomdict(self, specie=None):
    if specie is None:
      l=[]
      Nspecie=1
      while Nspecie <= self.nat:
        print(Nspecie,self.nat)
        l.append(self.atomdict(Nspecie))
        Nspecie+=1
      return l
    tz = self.tidy_zip(specie)
    vxx, vyy, vzz = tz[3][0], tz[3][1], tz[3][2]
    pax = tz[4]
    xax = np.array( pax[0] )
    yax = np.array( pax[1] )
    zax = np.array( pax[2] )
    freq = np.abs(f32(tz[1],tz[2]))
    label = tz[0]
    symbol = label[:2].strip()
    index = int( label.replace(symbol,'').strip() )
    totalefg = np.array(tz[6])
#   symmefg  = np.array(tz[7])
    return {
	'label': 	  tz[0],
	'name': 	  tz[0],
	'symbol':	  symbol,
        'index':	  index,
        'mdstep': 	  self.mdstep,
	'Cq':		  tz[1],
	'cq':		  tz[1],
	'Eta':		  tz[2],
	'eta':		  tz[2],
	'fq':		  freq,
	'Vii':		  tz[3],
	'vii':		  tz[3],
	'Vxx':		  vxx,
	'vxx':		  vxx,
	'Vyy':		  vyy,
	'vyy':		  vyy,
	'Vzz':		  vzz,
	'vzz':		  vzz,
	'principle_axes': pax,
	'axes':		  pax,
	'xaxis':	  xax,
        'x':		  xax,
	'yaxis':	  yax,
	'y':		  yax,
	'zaxis':	  zax,
	'z':		  zax,
	'Q':		  tz[5],
	'efg':         totalefg,
	'tensor':      totalefg,
	'total_efg':   totalefg,
#	'symmetrized_efg':symmefg,
        'is_right_handed': is_right_handed(xax,yax,zax)
}


  def atom(self, specie):
    return asobject(self.atomdict(specie))
		

  def compute_eta(self,atomic_specie_index):
    i = atomic_specie_index
    return (self.atom(i).vxx - self.atom(i).vyy)/self.atom(i).vzz


  @property
  def atoms(self):
    return [None] + [ self.atom(i) for i in range(1,self.nat + 1 ) ]    


  @property
  def data(self):
    return atoms 

