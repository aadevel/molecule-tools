# Get a segment of a protein structure
#!/usr/bin/env python

#usage: python GetResidueXtoY.py pdbfile lowrange highrange

import sys
sys.path.insert(0, '/home/aashish/protlib2/lib/python')

from Bio.PDB import *
from Bio.PDB.Fold import *

stru=PDBParser().get_structure('x',sys.argv[1])

i=int(sys.argv[2])
j=int(sys.argv[3])

for model in stru:
    for chain in model:
        for res in chain:
            id=res.id
            if id[1] < i or id[1] > j: 
                child=chain.child_dict[id]
                child.detach_parent()
                del chain.child_dict[id]
                chain.child_list=chain.child_dict.values()
                chain.child_list.sort(chain._sort)

w = PDBIO()
w.set_structure(stru)
renum=sys.argv[1]+ '-' + str(i) + 'to' + str(j) + '.pdb'
w.save(renum)
