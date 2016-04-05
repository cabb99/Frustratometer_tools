#!/usr/bin/python
import subprocess, shutil, os, glob
Temp_directory='/home/cab22/Programs/CFrustratometer/Temp'
Selector='pdb_subset.py'
Cleaner='/home/cab22/Programs/mmtsb/perl/convpdb.pl'
PDB2lammps='/home/cab22/Programs/Frustratometer/PdbCoords2Lammps'
strideexec='/home/cab22/Programs/Frustratometer/stride/stride'
stride2ssweight='/home/cab22/Programs/Frustratometer/stride2ssweight.py'
parameters='/home/cab22/Programs/Frustratometer/para/'
awsem_para='/home/cab22/Programs/CFrustratometer/frust_fix_backbone_coeff.data'
lammps_exec='/home/cab22/Programs/CFrustratometer/exec/lmp_serial_552_mod_frust_elec'

def MakeDir(directory):    
    if not os.path.exists(directory):
        os.makedirs(directory)

def Download(pdb,directory):
    from pdb_download import pdbDownload
    pdb=[pdb.lower()]  
    pdbDownload(pdb)
    pdb=pdb[0]
    pdb_location='%s/%s.pdb'%(directory,pdb)
    shutil.move('%s.pdb'%pdb,'%s/%s.pdb'%(directory,pdb))
    return pdb_location
     
def Select_chain(pdb_file,chains):
    from pdb_subset import pdbSubset
    # Read file
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    # Pop out subset
    pdb, chain, residues = pdbSubset(pdb,chains,[0,0])
  
    # Write to file 
    g = open(pdb_file,'w')
    g.writelines(pdb)
    g.close()
    return pdb_file
        

def Clean(pdb_file):
    #-nohetero ignores HETATM recodrd
    #-renumber renumbers all the residues from 1
    #-out Specifies the output format, in this case generic
    #-fixcoo fixes C-terminal atoms  
    commands=[Cleaner,'-nohetero','-renumber','1','-out','generic','-fixcoo',pdb_file]
    pdb_data=subprocess.check_output(commands)
    g = open(pdb_file,'w')
    g.writelines(pdb_data)
    g.close()
    return pdb_file

def Renumber_resid(pdb_file):
    from Bio.PDB import *
    parser=PDBParser()
    structure=parser.get_structure('Renumbered',pdb_file)
    for model in structure:
        i=1
        for chain in model:
            for residue in chain:
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

def Frustration(pdb,chain='All',frustration='mutational',download=True,custom_aa_freq=False,clean=True,
lammps_exec=lammps_exec):
    print "Calculating frustration"
    pdb_loc=pdb    
    pdb_name=pdb_loc.split('/')[-1]
    pdb_cname=pdb_name    
    if pdb[-4:] == ".pdb":
        pdb_cname = pdb[:-4]

    #Make a directory
    print "Creating a directory"    
    directory='%s_%s_%s'%(pdb_cname,chain,frustration)    
    if clean:
        directory='%s/%s_%s_%s'%(Temp_directory,pdb_cname,chain,frustration)
    MakeDir(directory)
    
    #Obtain the pdb
    if download:
        "Downloading %s"%pdb_cname        
        pdb_loc=Download(pdb_cname,directory)
    else:
        shutil.copy(pdb_loc,'%s/%s'%(directory,pdb_name))
        pdb_loc='%s/%s'%(directory,pdb_name)
    
    return_dir=os.getcwd()    
    os.chdir(directory)
    pdb_name=pdb_loc.split('/')[-1] 
    
    print "Cleaning pdb" 
    if chain<>'All':
        Select_chain(pdb_name,chain)
    
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
            else:
                handle.write('%s\n'%line)

    #Write the parameters for the AWSEM simulation
    for filename in glob.glob(os.path.join(parameters, '*')):
        print "Copying",filename        
        shutil.copy(filename, '.')
    with open(awsem_para) as handle_in:
        with open('fix_backbone_coeff.data','w+') as handle_out:
            for line in handle_in:
                #print line,frustration                
                line=line.replace('configurational',frustration)
                handle_out.write(line)   


    #Recalculate the ssweight
    print "Recalculating SSweight"    
    SSweight(pdb_name)

    #Run the frustratometer
    print "Running AWSEM"    
    with open('%s.in'%pdb_cname) as lin:
        subprocess.call([lammps_exec],stdin=lin)

    os.chdir(return_dir)
    #return frustration_data
    

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Analyzes the frustration of the residues of a protein')
    parser.add_argument('PDB', type=str, nargs='*', help='PDB files to be analyzed')
    parser.add_argument('-c','--chain',nargs='*',default=['All'] ,action='store', help='Selects the chains for each PDB')
    parser.add_argument('-f','--frustration', choices=['singleresidue','configurational','mutational','All'],default='mutational', help='Sets the frustration type used in the analysis')
    parser.add_argument('--download',help='Downloads the PDB file if the file does not exists',action='store_true')
    parser.add_argument('--freq',nargs=1,help='Uses a custom aminoacid frequency from a csv table',action='store',default=None)
    #parser.add_argument('--awsem',nargs=1,help='Uses custom awsem parameters as input' ,action='store',default=None)
    parser.add_argument('--clean',help='Only returns the output if the frustration',action='store_true')
    parser.add_argument('--lammps',nargs=1,help='Uses a custom lammps executable tu run the frustratometer' ,default=[lammps_exec])
    #parser.add_argument('args', nargs=argparse.REMAINDER)
    
    args = parser.parse_args()
    print args
    
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
    if not args.download:
        Non_existent_files=[]        
        for p in set(args.PDB):
            try:
                temp=open(p)
                temp.close()
            except IOError:
                Non_existent_files+=[p]
        if len(Non_existent_files)>0:
            parser.error('The following PDB files do not exist:\n%s\nIf you want to download them use the --download option'%(','.join(Non_existent_files)))
    
    #If the frequency is specified, try to open it
    if args.freq<>None:
         try:
            temp=open(args.freq)
            temp.close(args.freq)
         except IOError:
            parser.error('Could not open %s'%args.freq)

    if args.frustration=='All':
        args.frustration=['singleresidue','configurational','mutational']
    elif type(args.frustration)==str:
        args.frustration=[args.frustration]
    
    for frustration in args.frustration:
        for pdb,chain in zip(args.PDB,args.chain):    
            print "Processing %s %s %s"%(pdb,chain,frustration)
            Frustration(pdb,chain,frustration,args.download,args.freq,args.clean,lammps_exec=args.lammps[0])

if __name__=='__main__':
    main()

    
