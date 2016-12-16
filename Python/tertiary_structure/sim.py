import random, math, os, imp, sys, argparse
from pyrosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, Pose, MoveMap, \
                    toolbox
from rosetta import protocols, core

#protocols.simple_moves.MinMover()
#protocols.simple_moves.RotamerTrialsMover()

#sys.path.append('/Users/WRShoemaker/PyRosetta4.Release.python27.mac.release-123/setup/pyrosetta/toolbox')

cleaning = imp.load_source('cleaning', '/Users/WRShoemaker/PyRosetta4.' \
    'Release.python27.mac.release-123/setup/pyrosetta/toolbox/cleaning.py')

#from toolbox import mutate_residue
#from toolbox import cleanATOM

#foo = imp.load_source('module.name', '/path/to/file.py')
#foo.MyClass()

args =  sys.argv
#pdb_file = '/Users/WRShoemaker/GitHub/Evol_16S/Evolution/data/3j28'
os.system('python make_rna_rosetta_ready.py 3j28.pdb')
pdb_file = '3j28_RNA'
#pdb_file = pdb_file.split('.')[0]
init(extra_options='-mute basic -mute core')
# Constants
PACK_RADIUS = 10.0
#Nucleic acids, notice there is no C
NAs = ('A','C','G','U')
#Number of mutations to accept
max_accept_mut = 1500
#Population size
N = 100
#Beta (temp term)
beta = 1
 #Prepare data headers
data = ['Variant,Rosetta Score,"delta-delta-G",Probability,Generation\n']

#Load and clean up pdb file
name=pdb_file + '.pdb'
cleaning.cleanATOM(name)
clean_name=pdb_file+".clean.pdb"
#rna_set = core.chemical.ChemicalManager.get_instance().residue_type_set("coarse_rna").get()
#rna_set = core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")

initial_pose = pose_from_pdb(clean_name)
print initial_pose
print toolbox.get_secstruct(initial_pose)
#Set up ScoreFunction
sf = get_fa_scorefxn()

#Set up MoveMap.
mm = MoveMap()
#change these for more or less flexability
mm.set_bb(True)
mm.set_chi(True)

#Pack and minimize initial pose to remove clashes.
pre_pre_packing_score = sf(initial_pose)
print "1"
task = standard_packer_task(initial_pose)
task.restrict_to_repacking()
task.or_include_current(True)
print "2"
pack_rotamers_mover = protocols.simple_moves.RotamerTrialsMover(sf, task)
pack_rotamers_mover.apply(initial_pose)
print "3"

min_mover = protocols.simple_moves.MinMover()
print "4"

min_mover.movemap(mm)
print "5"
min_mover.score_function(sf)
print "6"
min_mover.min_type('dfpmin_armijo_nonmonotone')
print "7"
min_mover.apply(initial_pose)
print "8"
post_pre_packing_score = sf(initial_pose)
print "9"

print post_pre_packing_score
#445792
