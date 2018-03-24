# Align structures in pymol
import sys,os

import Bio
import Bio.PDB
import numpy
from numpy import *

from Bio.SVDSuperimposer import SVDSuperimposer

#pymol string
pymolstr='load\ \%s,prot\;load\ T1ubq_ref.pdb,prot2\;color\ green,prot2\;align\ prot,prot2\;hide\ all\;show\ cartoon\;spectrum\ b\;cartoon\ putty\;viewport\ 300,300\;ray\ 800,800\;png\ \%s\.png\;orient\;mset\ 1\ x60\;mplay\;util.mroll\ 1,60,180\;mstop\;set\ ray_trace_frames=1\;mpng\ frame'
cent_out=sys.argv[1]
os.system('pymol -c -d ' + (pymolstr %(cent_out, cent_out)))
