'''

LSMauto is a python code for chaining multiple functions to simplify the use of LSM in LAMMPS.

Anass Harmal 
Worcester Polytechnic Institute

MIT license 
Copyright (c) 7th of February 2023

Please begin by reading the README.md file to prepare your system for this code.

'''

from datetime import datetime
start_time = datetime.now()

import os
import argparse
import LSMfunc

from LSMfunc import BnMGen
from LSMfunc import Img2Particle
from LSMfunc import SimscriptGen


##### argument parsing #####


parser = argparse.ArgumentParser(prog='LSMauto: A one command executable to run a Lattice Spring Modeling simulation')

# Arguments to generate image

parser.add_argument('-b', dest='b',type = float, default = 300, # number of pixels (with 1 pixel = 0.01mm)
                                                help='Brick size')       
parser.add_argument('-a', dest='a',type = float, default = 3,
                                                help='Aspect ratio of bricks')
parser.add_argument('-ov', dest='ov',type = float, default = 2,
                                                help='Brick layers overlap')
parser.add_argument('-an', dest='an',type = float, default = 0,
                                                help='angle in degrees')
parser.add_argument('-sp', dest='sp',type = float, default = 50, # number of pixels (with 1 pixel = 0.01mm)
                                                help='Soft phase thickness')                                                
parser.add_argument('--folder', default='output',
                                                help='output file name')

# Arguments to generate LAMMPS data file

parser.add_argument('-o', dest='outfile', metavar='data file', type=str) # is this needed?
parser.add_argument('-lx', dest='lx', metavar='Lx', type=float, default=1536, #needs to be 1536
                    help='model size in x dimension\t')
parser.add_argument('-ly', dest='ly', metavar='Ly',type=float, default=496,   #needs to be 496
                    help='model size in y dimension\t')
parser.add_argument('-lz', dest='lz', metavar='Lz',type=float, default=1,
                    help='model size in z dimension\t(default = 1)')
parser.add_argument('-s', dest='s', metavar='lc', type=float, default=4.0,
                    help='lattice constant = particle spacing \t') 
parser.add_argument('-tt', dest='types', metavar='types', type=int, default=3,
                    help='number of particle types')
parser.add_argument('-n',dest='nsize', metavar='len', nargs=2, type=float, default=[16, 400],
                    help='create notch with lengths in x and y dimension\t')

# Arguments to generate LAMMPS script
parser.add_argument('-stg1', dest='stg1', type=float, default=0.00046, 
                    help='Strength of the soft phase\t')
parser.add_argument('-stg2', dest='stg2', type=float, default=0.00046, 
                    help='Strength of the interface\t')
parser.add_argument('-stg3', dest='stg3', type=float, default=0.0004101, 
                    help='Strength of the hard phase\t')
parser.add_argument('-stf1', dest='stf1', type=float, default=0.000495, 
                    help='Stiffness of the soft phase\t')
parser.add_argument('-stf2', dest='stf2', type=float, default=0.000495, 
                    help='Stiffness of the interface\t')
parser.add_argument('-stf3', dest='stf3', type=float, default=0.0009, 
                    help='Stiffness of the hard phase\t')

args = parser.parse_args()

##### Preparing simulation #####

#producing an image of the structure

BnMGen(ar=args.a,bh=args.b,ov=args.ov,an=args.an,sp=args.sp)

#converting that image into a LAMMPS data file

filename = 'Bnm_sz{}_ar{}_an{}_ov{}_sp{}.png'.format(args.b,args.a,args.an,args.ov,args.sp)
    
# Create folder for image

finalname = f"{args.b}_{args.a}_{args.ov}_{args.an}_{args.sp}"

now = datetime.now()

pathtoinput = '{}/{}'.format(finalname,filename)

Img2Particle(input_file = pathtoinput,l_x = args.lx ,l_y = args.ly ,l_z = args.lz ,s_ = args.s ,t_ = 3, c1 = 8, c2 =100)

#Writing the simulation script
SimscriptGen(stiff1 = args.stf1,stiff2 = args.stf2,stiff3 = args.stf3,strength1 = args.stg1,strength2 = args.stg2,strength3 = args.stg3, directory = finalname)

scriptname =  f"bend_lmps_{args.stf1}_{args.stf2}_{args.stf3}_{args.stg1}_{args.stg2}_{args.stg3}.in"

# submiting simulation script using LAMMPS (make sure that LAMMPS is installed as lmp in your system)

os.chdir(f"{finalname}")

os.system(f"lmp -in {scriptname}") #replace script.in by a variable text

# renaming file 

os.chdir("..")

now = datetime.now()

new_name = finalname + "_" + now.strftime("%Y-%m-%d_%H-%M-%S")

os.rename(finalname,new_name)

#Printing run time

end_time = datetime.now()
print('Run time: {}'.format(end_time - start_time))
