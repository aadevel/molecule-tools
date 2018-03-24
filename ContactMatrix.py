#Create a contact matrix of residues in a pdb structure

import sys
import Bio.PDB
import numpy

from Bio.PDB import *

#argv 1 is pdt or pdb file from which to get the Hbond statistics for next round, arg2 is output energy matrix, arg3 if "init" makes a
#zero matrix, else still need something as arg3 to print out the energy matrix
#


Hmax=20
Hmin=2


def isNHbonded(residue_one, residue_two):
    if residue_one.has_id('HN') and residue_two.has_id('O'):
        diff=residue_one["HN"].coord - residue_two["O"].coord
        dist=numpy.sqrt(numpy.sum(diff*diff))
        v1=residue_one["HN"].get_vector()
        v2=residue_one["N"].get_vector()
        v3=residue_two["O"].get_vector()
        angle=calc_angle(v2,v1,v3)
        sdist=numpy.abs(residue_one.id[1]-residue_two.id[1])>3

#        print residue_one.id[1],residue_one.resname, residue_two.id[1], dist, angle
        if dist < 2.5 and angle > 2.26 and sdist:
            return 1
        else:
            return 0
    else:
        return 0


def isOHbonded(residue_one, residue_two):
    if residue_one.has_id('O') and residue_two.has_id('HN'):
        diff=residue_one["O"].coord - residue_two["HN"].coord
        dist=numpy.sqrt(numpy.sum(diff*diff))
        v1=residue_one["O"].get_vector()
        v2=residue_two["HN"].get_vector()
        v3=residue_two["N"].get_vector()
        angle=calc_angle(v1,v2,v3)
        sdist=numpy.abs(residue_one.id[1]-residue_two.id[1])>2

#        print residue_one.id[1],residue_one.resname, residue_two.id[1], dist, angle
        if dist < 3.5 and angle > 1.9 and sdist:
            return 1
        else:
            return 0
    else:
        return 0


def isContact(residue_one, residue_two,threshold):
    if residue_one.has_id('CB') and residue_two.has_id('CB'):
        diff  = residue_one["CB"].coord - residue_two["CB"].coord
        dist=numpy.sqrt(numpy.sum(diff*diff))
        v1=residue_one["CA"].get_vector()
        v2=residue_one["CB"].get_vector()
        v3=residue_two["CA"].get_vector()
        v4=residue_two["CB"].get_vector()
        angle1=calc_angle(v2,v1,v3)
        angle2=calc_angle(v4,v3,v1)
        sdist=numpy.abs(residue_one.id[1]-residue_two.id[1])>3
            
    #    if dist < threshold and dist > 1.5 and (angle1 < 1.57 or angle2 < 1.57) and sdist:
        if dist < threshold and dist > 1.5  and sdist:
            if sdist is 4:
                if (isNHbonded(residue_one,residue_two) or  isOHbonded(residue_one,residue_two)):
                    return 1
                else:
                    return 0
            else:
                return 1
        else:
            return 0
    else:
        return 0



def calc_residue_dist(residue_one, residue_two) :
	"""Returns the C-alpha distance between two residues"""
	diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
	return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

"""Returns a matrix of C-alpha distances between two chains"""
def calc_dist_matrix(chain_one, chain_two) :
	answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
	for row, residue_one in enumerate(chain_one) :
		for col, residue_two in enumerate(chain_two) :
			answer[row, col] = calc_residue_dist(residue_one, residue_two)
	return answer

'''
#returns integer contact matrix
def calc_contact_matrix(chain_one, chain_two,threshold) :
	answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
	for row, residue_one in enumerate(chain_one) :
		for col, residue_two in enumerate(chain_two) :
			#answer[row, col] = calc_residue_dist(residue_one, residue_two)
            if(isHbonded(residue_one, residue_two)) : answer[row,col]=1
            #if (answer[row,col]) < threshold: answer[row,col]=1
			else: answer[row,col]=0
	return answer
'''


#returns Hbond matrix
def calc_Cmatrix(chain_one, chain_two) :
        answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
        for row, residue_one in enumerate(chain_one) :
            for col, residue_two in enumerate(chain_two) :   
                #if(isNHbonded(residue_one, residue_two) or isOHbonded(residue_one,residue_two)) : answer[row,col]=1
                if(isContact(residue_one, residue_two, 7.0 )) : answer[row,col]=1
                else: answer[row,col]=0

   #Remove comment to normalize
        #if two H bonds, each half - ie Normalize Hbonds
        #for i in range(len(chain_one)):
        #    s=answer[i].sum()
        #    for j in range(len(chain_one)):
        #        if answer[i][j] > 0:
         #           answer[i][j]=answer[i][j]/s
        
        print answer
        return answer


structure=Bio.PDB.PDBParser().get_structure("1ubq", sys.argv[1])
model=structure[0]
cm = numpy.zeros((len(model["A"]), len(model["A"])), numpy.float)


count=0
for model in structure.get_list():
#	dm=dm+calc_dist_matrix(model["A"],model["A"])
	cm=cm + calc_Cmatrix(model["A"],model["A"]) 
	count += 1
#	print count



x=len(model["A"])
out_path = sys.argv[2]
out_file = open(out_path, 'w')
for j in range(x): 
    cmjsum=cm[j].sum()
    pHb=cmjsum/count
    sigma=1-pHb
    print cmjsum
    #cmjnzsum=(cm[j]!=0).sum()
    for k in range(x):
        if sys.argv[3]=="init": val = 0.0
        else:
            val=0.0
            if pHb > 0:
		print j,k,cm[j][k],count,cm[j][k]/count
                if (cm[j][k]/count) > 0.1:
                    val = 1.0 * cm[j][k]/count
        out_file.write(str(val) + ' ')
        #out_file.write(str(val) + ' ' + str(cm[j][k]) + ' ' + str(cmjsum) + "===")
    out_file.write('\n')


