# Track the torsional angles of the amino acids in protein dynamics trajectory

#!/usr/bin/env python

import sys

from Bio.PDB import *
#from Bio.PDB.Polypeptide import *
from Bio.PDB.Fold import *
import time
import StringIO

import math
def deg(rad):
    if rad is None:
        return None
    ang = rad * 180 / 3.14
    while ang > 180:
        ang = ang - 360
    while ang < -180:
        ang = ang + 360
    return ang


def basin(phi,psi):
    if (psi > 50 and psi < 180 and phi < -100 and phi > -180):
        return "E"
    elif (psi > -180 and psi < -100 and phi < -180 and phi > -100):
        return "E"
    elif phi < 0 and phi > -180 and psi > -100 and psi < 50:
        if phi > -71 and phi < -57 and psi > -49 and psi < -35:
            return "H"
        else:
            return "Hb"
    elif phi < 0 and phi > 100 and psi > 50 and psi < 180:
        return "PP2"
    elif phi < 180 and phi > 0 and psi > -50 and psi < 100:
        return "RH"
    else:
        return "C"


# E = 1, Hcore = 4, Hbasin = 3, PP2 = 2, RH= 5, rest = 6
def basinNumeric(phi,psi):
    if (psi > 50 and psi < 180 and phi < -100 and phi > -180):
        return 1
    #elif (psi > -180 and psi < -100 and phi < -180 and phi > -100):
    #    return 1
    elif phi < 0 and phi > -180 and psi > -100 and psi < 50:
        #if phi > -71 and phi < -57 and psi > -49 and psi < -35:
        if phi > -78 and phi < -50 and psi > -59 and psi < -25:
            return 4
        else:
            return 3
    elif phi < 0 and phi > 100 and psi > 50 and psi < 180:
        return 2
    elif phi < 180 and phi > 0 and psi > -50 and psi < 100:
        return 5
    else:
        return 6





parser = PDBParser()
protein = PDBParser().get_structure('xx',sys.argv[1])

ppb=PPBuilder()


import Bio.PDB
for model in protein:
    s=''
    for chain in model:
        pp=ppb.build_peptides(chain)
        for poly_index, poly in enumerate(pp):
            phi_psi=poly.get_phi_psi_list()
            for ri, res in enumerate(poly):
                phi= deg(phi_psi[ri][1])
                psi= deg(phi_psi[ri][2])
                s=s+str(basinNumeric(phi,psi))
                print ri,phi,psi,basin(phi,psi),basinNumeric(phi,psi)
