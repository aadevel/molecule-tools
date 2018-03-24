#pdb file to png output 
import sys,os

import Bio
import Bio.PDB
import numpy
from numpy import *

from Bio.SVDSuperimposer import SVDSuperimposer

pymolstr='run\ color_by_restype.py\;load\ \%s,prot\;hide\ all\;show\ cartoon\;bg_color\ white\;color_by_restype\ prot\;ray\ 800,800\;viewport\ 300,300\;png\ \%s\.png\;'
cent_out=sys.argv[1]
os.system('pymol -c -d ' + (pymolstr %(cent_out, cent_out)))
