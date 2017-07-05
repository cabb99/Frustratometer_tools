#!/usr/bin/python
import subprocess, shutil, os, glob

#Locate installation folder
path=os.path.realpath(__file__)
path='/'.join(path.split('/')[:-1])
print 'Running frustratometer from', path


#Cleaner='/home/cab22/Programs/mmtsb/perl/convpdb.pl'
PDB2lammps='%s/PdbCoords2Lammps'%path
strideexec='%s/stride/stride'%path
stride2ssweight='%s/stride2ssweight.py'%path
parameters='%s/para/'%path
awsem_para='%s/frust_fix_backbone_coeff.data'%path
lammps_exec='%s/exec/lmp_serial_552_mod_frust_elec'%path

def MakeDir(directory):    
    if not os.path.exists(directory):
        os.makedirs(directory)

def Clean(pdb_file):
    #-nohetero ignores HETATM recodrd
    #-renumber renumbers all the residues from 1
    #-out Specifies the output format, in this case generic
    #-fixcoo fixes C-terminal atoms  
    #commands=[Cleaner,'-nohetero','-renumber','1','-out','generic','-fixcoo',pdb_file]
    #pdb_data=subprocess.check_output(commands)
    #g = open(pdb_file,'w')
    #g.writelines(pdb_data)
    #g.close()
    return pdb_file

def Renumber_resid(pdb_file):
    from Bio.PDB import PDBParser
    from Bio.PDB import PDBIO
    parser=PDBParser()
    structure=parser.get_structure('Renumbered',pdb_file)
    for model in structure:
        i=1
        for chain in model:
            for residue in chain:
                if residue.id == (' ', i, ' '):
                    pass
                else:
                    residue.id = (' ', i, ' ')
                i=i+1
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def Pdb2Lammps(pdb_file,name):
    commands=[PDB2lammps,pdb_file[:-4], name]
    f = open(os.devnull, 'w')
    return subprocess.check_output(commands,stderr=f)

def SSweight(pdb_file):
    commands=[strideexec,pdb_file]
    stride=subprocess.check_output(commands)    
    g = open('ssweight.stride','w')
    g.writelines(stride)
    g.close()
    commands=['python',stride2ssweight]
    ssweight=subprocess.check_output(commands)  
    g = open('ssweight','w')
    g.writelines(ssweight)
    g.close()

def Frustration(pdb,chain,frustration,args):
    awsem_para='%s/frust_fix_backbone_coeff.data'%path
    lammps_exec=args.lammps[0]
       
    
    print "Calculating frustration"
    pdb_loc=pdb    
    pdb_name=pdb_loc.split('/')[-1]
    pdb_cname=pdb_name    
    if pdb[-4:] == ".pdb":
        pdb_cname = pdb[:-4]

    #Make a directory
    print "Creating a directory"    
    directory='%s_%s_%s'%(pdb_cname,chain,frustration)    
    MakeDir(directory)
    
    #Obtain the pdb
    shutil.copy(pdb_loc,'%s/%s'%(directory,pdb_name))
    pdb_loc='%s/%s'%(directory,pdb_name)
    
    return_dir=os.getcwd()    
    os.chdir(directory)
    pdb_name=pdb_loc.split('/')[-1] 
    
    print "Cleaning pdb" 
    if chain<>'All':
        #Select_chain(pdb_name,chain)
        pass
    
    #Clean the pdb    
    Clean(pdb_name)

    Renumber_resid(pdb_name)    
    
    #Convert the pdb to lammps input
    print "Creating lammps input"    
    Pdb2Lammps(pdb_name,pdb_cname)
    with open('%s.in'%pdb_cname) as handle:
        in_data=handle.read()
    with open('%s.in'%pdb_cname,'w+') as handle:
        for line in in_data.split('\n'):
            if line=='run		10000':
                handle.write('run  0\n')
            elif line=='neighbor\t30 bin':
                handle.write('neighbor 15 bin\n')
            else:
                handle.write('%s\n'%line)

    #Write the parameters for the AWSEM simulation
    for filename in glob.glob(os.path.join(parameters, '*')):
        print "Copying",filename        
        shutil.copy(filename, '.')
    os.chdir(return_dir)
    if args.coeff<>None:
        awsem_para=args.coeff
    with open(awsem_para) as handle_in:
        with open('%s/fix_backbone_coeff.data'%directory,'w+') as handle_out:
            for line in handle_in:
                #print line,frustration                
                line=line.replace('configurational',frustration)
                handle_out.write(line)
    
    #Overwrite parameters
    if args.gamma<>None:
        shutil.copy(args.gamma, '%s/gamma.dat'%directory)
    if args.burial_gamma<>None:
        shutil.copy(args.burial_gamma, '%s/burial_gamma.dat'%directory)
    os.chdir(directory)
    #Recalculate the ssweight
    print "Recalculating SSweight"    
    SSweight(pdb_name)

    #Run the frustratometer
    print "Running AWSEM"    
    with open('%s.in'%pdb_cname) as lin:
        subprocess.call([lammps_exec],stdin=lin,stderr=open('lammps.err','w+'),stdout=open('lammps.out','w+'))

    os.chdir(return_dir)
    #return frustration_data
    

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Analyzes the frustration of the residues of a protein')
    parser.add_argument('PDB', type=str, nargs='*', help='PDB files to be analyzed')
    #parser.add_argument('-c','--chain',nargs='*',default=['All'] ,action='store', help='Selects the chains for each PDB')
    parser.add_argument('-f','--frustration', choices=['singleresidue','configurational','mutational','All'],default='mutational', help='Sets the frustration type used in the analysis')
    parser.add_argument('--gamma',help='Specify a gamma file')
    parser.add_argument('--burial-gamma',help='Specify a burial gamma file') 
    parser.add_argument('--coeff',help='Specify a fix backbone coeficient file') 
    parser.add_argument('--lammps',nargs=1,help='Uses a custom lammps executable tu run the frustratometer' ,default=[lammps_exec])
    #parser.add_argument('args', nargs=argparse.REMAINDER)
    
    args = parser.parse_args()
    print args
    args.chain='All'
    
    #Length of chains can be the same as the length of the pdb, 1 or 0.
    if len(args.chain)==1:
        args.chain=[args.chain[0] for i in range(len(args.PDB))]
    elif len(args.PDB)==1:
        args.PDB=[args.PDB[0] for i in range(len(args.chain))]
    elif len(args.chain)==len(args.PDB):
        pass
    else:
        parser.error('Number of chain selections do not correspond to number of PDBs')

    #All files should exist if download is not indicated
    Non_existent_files=[]        
    for p in set(args.PDB):
        try:
            temp=open(p)
            temp.close()
        except IOError:
            Non_existent_files+=[p]
    if len(Non_existent_files)>0:
        parser.error('The following PDB files do not exist:\n%s\n'%(','.join(Non_existent_files)))
    
    if args.frustration=='All':
        args.frustration=['singleresidue','configurational','mutational']
    elif type(args.frustration)==str:
        args.frustration=[args.frustration]
    
    for frustration in args.frustration:
        for pdb,chain in zip(args.PDB,args.chain):    
            print "Processing %s %s %s"%(pdb,chain,frustration)
            Frustration(pdb,chain,frustration,args)

if __name__=='__main__':
    main()

    
