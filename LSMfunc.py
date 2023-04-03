'''

LSMfunc is a module that includes functions used to automate Lattice Spring Modeling

Anass Harmal 
Worcester Polytechnic Institute

MIT license 
Copyright (c) 7th of February 2023

'''


def BnMGen(ar,bh,ov,an,sp): #aspect ratio, brick height, width, height, overlap, angle, soft phase thickness
    '''

    BnMGen generates a greyscale image of a brick and mortar composite

    Anass Harmal 
    Worcester Polytechnic Institute

    MIT license 
    Copyright (c) 6th of January 2023

    '''
    import os
    import argparse
    import math
    import PIL
    from PIL import Image, ImageChops, ImageDraw
    
    # number of pixels for image (represent 0.01mm/pixel for our samples)
    w=15400
    h=5000
    
    
    # inputs 
    aspect_ratio = ar
    brick_height = bh # add to arg parse
    brick_width = brick_height*aspect_ratio
    width = w
    width = width + brick_width
    height = h # add to arg parse (10 pixels/mm) (pixels added in the x direction to eliminate artifacts later) 
    # overlap (overlap in only 1/2 for bowtie geometry)
    overlap = ov #add to arg parse (1/overlap)
    #angle
    alpha0 = an #bowtie incline angle (input in degrees) 
    alpha = (math.pi/180)*alpha0 #converting alpha from degrees to rad
    # Spacing in x and y
    spacing_x = sp # add to parse
    sp_y = sp #spacing following the y direction
    spacing_y = sp_y - (brick_width/2)*(math.tan(alpha))  # relationship to keep the same spacing no matter what the angle is
    # Set the number of rows and columns
    cols, rows = round((width/brick_width)) , round((height/(brick_height-2*(brick_width/2)*math.tan(alpha))))


    #defining white image
    image = Image.new('L', (width, height), 'white') 
    # Get a drawing context
    draw = ImageDraw.Draw(image)


    # drawing brick-and-mortar structure
    for row in range(rows):
        for col in range(cols):
            if row % 2 != 0: # to offset every other layer via the overlap
                x1 = round(brick_width/overlap) + col * (brick_width + spacing_x)
                y1 =  row * (brick_height + spacing_y)
                x2, y2 = x1 + brick_width/2 , y1 + (brick_width/2)*(math.tan(alpha))  # first point with an angle
                x3, y3 = x1 + brick_width, y1 
                x4, y4 = x3, y3 + brick_height
                x5, y5 = x2 , y4 - (brick_width/2)*(math.tan(alpha))  # second point with an angle
                x6, y6 = x1, y4
            else:  
                x1, y1 = col * (brick_width + spacing_x), row * (brick_height + spacing_y)
                x2, y2 = x1 + brick_width/2 , y1 + (brick_width/2)*(math.tan(alpha))  # first point with an angle
                x3, y3 = x1 + brick_width, y1 
                x4, y4 = x3, y3 + brick_height
                x5, y5 = x2 , y4 - (brick_width/2)*(math.tan(alpha))  # second point with an angle
                x6, y6 = x1, y4
            # sequence of points to draw 
            points = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5), (x6, y6)]
            # Draw the bowtie
            draw.polygon(points, fill="black")  
    #image.show() #comment before image production 
    nimage = Image.new('L',(image.width,image.height+round(sp_y)),'white') 
    nimage.paste(image,(0,round(sp_y)))
    # flip the image 180 degrees (because we can the soft layer at the bottom to activate crack shilding)
    nimage = nimage.rotate(180) # uncomment when done
    # crop image to delete unwanted artifact in the end (comment to see)
    nimage = nimage.crop((0,sp_y,round(width-brick_width),height+sp_y))
    # Save image
    filename = 'Bnm_sz{}_ar{}_an{}_ov{}_sp{}.png'.format(brick_height,aspect_ratio,alpha0,overlap,spacing_x)
    
    # Create folder for image
    finalname = f"{bh}_{ar}_{ov}_{an}_{sp}"
    if not os.path.exists(finalname):
        os.makedirs(finalname)


    nimage.save(f"{finalname}/{filename}") 


def Img2Particle(input_file,l_x,l_y,l_z,s_,t_,c1,c2):
    '''

    Img2Particle
    
    Img2Particle is a simple, in-house developed python code for conversion from
    gray image to hexagonal packing particle model in LAMMPS Data format.
    
    MIT License
    
    Copyright (c) 6th/June/2020 Yuan Chiang
    
    Modified by Anass Harmal to simulate bending 12/20/2022
    
    3-point bending: top roller is set up on lammps as a fix indent and bottom rollers are 
    replaced by forcing a 5% of molecules on each side to zero (only zeroing y direction for the ends of the specimen)
    
    '''
    import os
    import argparse
    import numpy as np
    import math
    from matplotlib import pyplot as plt
    from PIL import Image
    from datetime import datetime
    
    infile = input_file
    
    lx = l_x
    ly = l_y
    lz = l_z
    s = s_
    t = t_
    crack1 = c1
    crack2 = c2
    
    def import_img(infile):
        img = Image.open(infile).convert('L')
        img.load()
        array = np.asarray(img, dtype="int32")
        return array
    
    img = import_img(infile)
    print("Load image from file: {}".format(infile))
    print("\timage size: {}".format(img.shape))
    
    # ===== Define Model
    
        
    print("Set up model...")
    print("\tlx = {:f}\tly = {:f}\tlz = {:f}".format(lx, ly, lz))
    print("\tlattice constant = {:f}".format(s))
    
    # ===== Define Lattice
    
    # create hexagonal packing unit cell
    unit = s*np.array([[math.sqrt(3.),0,0],[0,1,0],[0,0,0]], dtype=float) # distance unit, changed fro mthe original to rotate the triangle lattice
    
    # Defining the number of lattices
    
    nx = int(round(lx/unit[0,0]))
    ny = int(round(ly/unit[1,1]))
    nz = 1
    
    # Printing number of lattices in the 3 directions
    
    print("Hexagonal lattice:\n{}".format(unit))
    print("\tnumber of lattices in x: %3d " % (nx)) 
    print("\tnumber of lattices in y: %3d " % (ny))
    print("\tnumber of lattices in z: %3d " % (nz))
    
    
    # ===== Create Atoms
    
    dtype = dict(zip(['id', 'molecule', 'type', 'q', 'x', 'y', 'z', 'ix', 'iy', 'iz'], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])) #a dictionary that uses zip: meaning it will match id for examp
    
    natoms = int(2*nx*ny) #total number of atoms (multiplying by 2 because the vector creates atoms in two layers at a time)
    
    #create memory space for atoms (each atom has a vector of 10)
    atoms = np.zeros((natoms,10))
    
    #create id for each molecule
    atoms[:,dtype.get('id')] = np.arange(1,natoms+1) #id needs to be aranges
    
    atoms[:,dtype.get('molecule')] = np.ones((natoms,1)).reshape(-1) #fill molecules memory space with ones and reshape into a 1D array (vector)
    
    
    id_ = 0 # different from id
    colx = dtype.get('x')
    coly = dtype.get('y')
    colz = dtype.get('z')
    for i in range(nx):
        for j in range(ny):
            for k in range(nz): #defining two atoms at a time, unit serves to define the diance between atoms and 0.5 additions on the second atom generate a position at 45 degree ang
                
                atoms[id_,colx:colz+1] = np.dot(np.transpose(unit),[[i],[j],[k]]).reshape(-1) #dot product
                id_ = id_ + 1
                atoms[id_,colx:colz+1] = np.dot(np.transpose(unit),[[i+0.5],[j+0.5],[k]]).reshape(-1)
                id_ = id_ + 1
    
    
    print("Created atoms: {:d}".format(natoms))
    
    # ===== Change Types
    
    ntypes = t # types is in inputs, it is the number of types in the system
    
    rows, cols = img.shape # pixels in row , pixels in columns
    colt = dtype.get('type') #getting the type column in dtype which is where we defined the zip
    
    hist, bin_edges = np.histogram(img.reshape(-1),bins=ntypes) # hist will be the bumber of black dots and the number of whites ; bin edges will be 0 (black) 127.5 (middle) and 255 (
    
    print(hist, bin_edges)
    
    # for loops start from 0 until range - 1 
    
    for i in range(natoms):
        x = int(round(atoms[i,colx]/lx*cols)); # x = (atom position in x / length of system in x direction) * number of pixels in columns #transfering x from atoms ref to pixels ref
        y = rows - 1 - int(round(atoms[i,coly]/ly*rows)); # same as x but starting from the other side of y (not sure why it starts from the other side)
        #do not use this line!# y =  int(round(atoms[i,coly]/ly*rows)); #it can't be like this because it generates a mirror image (inverted) which is not correct
    
        for b in range(ntypes): # b navigates the types (ntypes is the number of types)
            if img[y,x] >= bin_edges[b]: # we read it [y,x] because of the way the matrice is stored by python (row major order)
                # bin_edges of b means it navigates the colors, this if checks if the color of the pixel is larger than the bin edge
                atoms[i,colt] = ntypes - b #this gives every image type 2 (because they're all equal or more than 0), then it gives type 1 to all the pixels that are not black (white)
    
    # Border conditions to freeze left and right ends
    for i in range(natoms):
        if atoms[i,colx] < (0.02*lx): 
            atoms[i,colt]=2 #Top roller
    
    for i in range(natoms):
        if atoms[i,colx] > (0.98*lx): 
            atoms[i,colt]=2 #Top roller
    
    ntypes = int(max(atoms[:,colt])) 
    
    # ===== Create Notch
    
    wn = crack1 # width of notch
    hn = crack2 # height of notch
    
    indices = np.logical_not(np.logical_and(np.logical_and(atoms[:,colx] < lx/2.0 + wn/2.0, atoms[:,colx] > lx/2.0 - wn/2.0), atoms[:,coly] < hn)) #flaging atoms within to be deleted
    #indices is only true for atoms that are outside the notch
    
    atoms = atoms[indices] #only keeps the atoms and discards the notch
    natoms = atoms.shape[0] #shape[0] gives the number of elements in the first direction - this gives the new number of atoms
    atoms[:,dtype.get('id')] = np.arange(1,natoms+1) # Rearranging the atoms since we deleted some
    
    #color is just used for plotting
    color = [str(item/(ntypes+1)) for item in atoms[:,colt]] #item will take values of types of atoms, meaning the item will be either type - then the color will take the string = ite
    
    plt.figure(figsize=[12.8, 9.6], dpi=600, facecolor=None, edgecolor=None, frameon=True)
    
    plt.subplot(121)
    plt.imshow(img,cmap='gist_gray')
    
    plt.subplot(122)
    plt.scatter(atoms[:,colx],atoms[:,coly],s=0.1,marker="o",c=color)
    ax = plt.gca()
    ax.set_aspect(1.0)
    plt.savefig(infile+'_conversion.png')
    
    # this creates output file
    
    outfile = os.path.splitext(infile)[0] + '.data'
    
    print('Write data file: {}'.format(outfile))
    
    # writing in data file
    with open(outfile, "w") as outfile:
        now = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        outfile.write("LAMMPS data created via Img2Particle at "+now+"\n")
        outfile.write("\n%d atoms\n" % natoms)
        outfile.write('# use create_bonds command in lammps to create springs\n')
        outfile.write("\n%d atom types\n" % ntypes)
        
        outfile.write("%d bond types\n" % int(ntypes+1)) 
    
        outfile.write("\n%12.5E %12.5E xlo xhi\n" % (np.min(atoms[:,colx]),np.max(atoms[:,colx])))
        outfile.write("%12.5E %12.5E ylo yhi\n" % (np.min(atoms[:,coly]),np.max(atoms[:,coly])))
        outfile.write("%12.5E %12.5E zlo zhi\n" % ((np.min(atoms[:,colz])-lz/2),(np.max(atoms[:,colz])+lz/2)))
    
        outfile.write("\nMasses\n\n")
        for i in range(ntypes):
            outfile.write("%5d\t%9.3E\n" % (i+1,1))
        
        outfile.write("\nAtoms # bond\n\n")
        for i in range(natoms):
                outfile.write("%5d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%d\n" % (atoms[i,dtype.get('id')],
                                                                              atoms[i,dtype.get('molecule')],
                                                                              atoms[i,dtype.get('type')],
                                                                              atoms[i,dtype.get('x')]*0.0001,
                                                                              atoms[i,dtype.get('y')]*0.0001,
                                                                              atoms[i,dtype.get('z')],
                                                                              atoms[i,dtype.get('ix')],
                                                                              atoms[i,dtype.get('iy')],
                                                                              atoms[i,dtype.get('iz')]))


def SimscriptGen(stiff1,stiff2,stiff3,strength1,strength2,strength3,directory):
    '''
    SimscriptGEN
    
    Code to generate LAMMPS script for bending.
    
    Anass Harmal 
    Worcester Polytechnic Institute
    MIT license 
    
    Copyright (c) 24th of January 2023
    '''
    
    import argparse
    import os
    
    
    # Working in current directory
    
    # Finding .data file
    for file in os.listdir(directory):
        # Check if the file is a .data file
        if file.endswith('.data'):
            # Set the path of the .data file
            datafile = os.path.splitext(file)[0]
            # Exit the loop as soon as the .data file is found
            break
   
    # Stiffness
    stf1 = stiff1
    stf2 = stiff2
    stf3 = stiff3
    # Strength 
    stg1 = strength1
    stg2 = strength2
    stg3 = strength3   
    
    
    #folder_name = f"{stf1}_{stf2}_{stf3}_{stg1}_{stg2}_{stg3}"
    #if not os.path.exists(folder_name):
    #    os.makedirs(folder_name)
    # Name of script
    filename = f"bend_lmps_{stf1}_{stf2}_{stf3}_{stg1}_{stg2}_{stg3}.in" # Create LAMMPS filename
    filepath = os.path.join(directory,filename)
    ##### LAMMPS script writing #####
    # write the part that includes inputs
    LAMMPS_script = """
    ######################################################################################################
    
    ##  LAMMPS code to simulate 3-point bending for bioinspired brick-and-mortar geopolymer composites  ##
    
    #####  VARIABLES #####
    
    # materials variables
    variable        stf1 equal {} #stiffness of the soft phase
    variable        stf2 equal {} #stiffness of the interface
    variable        stf3 equal {} #stiffness of the hard phase
    variable        stg1 equal {} #strength of the soft phase
    variable        stg2 equal {} #strength of the interface
    variable        stg3 equal {} #strength of the hard phase
    variable        lmpsfile equal {} # The name of the lammps executable created
    variable        name index {} # name of the input file
    variable        drctr index {}
    """.format(stf1,stf2,stf3,stg1,stg2,stg3,filename,datafile,directory)
    LAMMPS_script2 = """
    
    # LAMMPS parameters
    variable        dist equal 0.0004    # depends on the distance we set in Img2Particle
    variable        bcheck equal 100 # numbers of checks for broken bonds
    
    ##### LAMMPS CONFIGURATION #####
    
    units           si
    dimension       2
    boundary        s s p
    atom_style      bond
    ##### aquiring data #####
    read_data       ${name}.data extra/bond/per/atom 6 extra/special/per/atom 6
    
    ##### Bonds #####
    
    # groups for hard phase, soft phase, and border conditions
    group           soft type 1
    group           sides type 2
    group           hard type 3
    group           mobile subtract all sides
    group           zeroxy id 1
    pair_style      zero ${dist} nocoeff
    pair_coeff      * *
    neighbor        ${dist} bin
    create_bonds    many hard hard   1 $(v_dist*0.1) $(v_dist*1.1)
    create_bonds    many soft soft   2 $(v_dist*0.1) $(v_dist*1.1)
    create_bonds    many hard soft   3 $(v_dist*0.1) $(v_dist*1.1)
    create_bonds    many hard sides  4 $(v_dist*0.1) $(v_dist*1.1)
    create_bonds    many soft sides  4 $(v_dist*0.1) $(v_dist*1.1)
    create_bonds    many sides sides 4 $(v_dist*0.1) $(v_dist*1.1)
    write_data		${name}_spring.data
    
    ##### Force field #####
    pair_style      zero ${dist} nocoeff
    pair_coeff      * *
    bond_style      harmonic
    bond_coeff      1 ${stf3} ${dist}
    bond_coeff      2 ${stf1} ${dist}
    bond_coeff      3 ${stf2} ${dist}
    bond_coeff      4 ${stf1} ${dist}
    neighbor        ${dist} nsq
    neigh_modify    delay 10
    comm_modify     cutoff $(100*v_dist)
    fix             bondbreak1 all bond/break ${bcheck} 1 ${stg3} #critical stretch
    fix             bondbreak2 all bond/break ${bcheck} 2 ${stg1} #critical stretch
    fix             bondbreak3 all bond/break ${bcheck} 3 ${stg2} #critical stretch
    fix             0 all enforce2d

    ##### Computing #####
    compute         pea all pe/atom bond #potential energy of atoms
    compute         satom all stress/atom NULL # stress of atoms #we can compute stress in each group of atoms as well, might be interesting to compute
    compute         stred all reduce sum c_satom[*]
    compute         bonds all property/local btype batom1 batom2
    compute         bonds2 all bond/local dist engpot
    
    ##### simulation #####
    fix             1 sides setforce NULL 0.0 NULL #border condition on sides of sample to stop y motion 
    fix             0xy zeroxy setforce 0.0 0.0 0.0 # to fix the body in space, which helps avoid divergence of the model
    timestep        1
    fix             nve all nve
    thermo          1000
    thermo_style	custom step time temp etotal pe ke ebond pxx pyy lx ly 
    restart         100000 restart.*
    variable        i loop 120
    variable        stp equal step
    
    #load-min loop      
    label           loop
    variable        cx equal 0.076875
    variable        cy equal $(0.051+0.01-0.00000001*v_stp) 
    variable        force_y equal fy[27725]

    fix				ss all print 10000 "${stp}   ${cx}    ${cy}    ${force_y}" append ss.txt screen no title " " # prints indenter position every step
    
    dump            atom all custom 10000 atom.dump.* id type x y c_pea c_satom[1] c_satom[2] c_satom[4] fx fy
    dump_modify		atom first yes
    dump            bond all local 10000 bond.dump.* index c_bonds[*] c_bonds2[*]
    fix             Indenter all indent 0.01 cylinder z v_cx v_cy 0.01    #try different sizes
    fix             vsc all viscous 0.0001
    run             10000
    undump          atom
    undump          bond
    unfix           Indenter 
    unfix           vsc
    next            i
    
    """
    LAMMPS_script3 = """
    jump            {} loop
    """.format(filename)
    with open(filepath, 'w') as file:
        file.write(LAMMPS_script + "\n")
        file.write(LAMMPS_script2 + "\n")
        file.write(LAMMPS_script3 + "\n")
    