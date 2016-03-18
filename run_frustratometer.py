import subprocess, shutil, os
Temp_directory='/home/carlos/Programas/CFrustratometer/Temp'
Selector='pdb_subset.py'
Cleaner='/home/carlos/Programas/mmtsb/perl/convpdb.pl'
PDB2lammps='/home/carlos/Programas/Frustratometer/PdbCoords2Lammps'
strideexec='/home/carlos/Programas/Frustratometer/stride/stride'
stride2ssweight='/home/carlos/Programas/Frustratometer/stride2ssweight.py'
parameters='/home/carlos/Programas/Frustratometer/para'


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
    commands=[stride2ssweight]
    ssweight=subprocess.check_output(commands)  
    g = open('ssweight','w')
    g.writelines(ssweight)
    g.close()

def Frustration(pdb,chain='All',frustration='mutational',download=True,custom_aa_freq=False,clean=True):
    print "Calculating frustration"
    pdb_loc=pdb    
    pdb_name=pdb_loc.split('/')[-1]
    pdb_cname=pdb_name    
    if pdb[-4:] == ".pdb":
        pdb_cname = pdb[:-4]
    chain='AC'
    frustration='mutational'
    download=True
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
    
    #Open the pdb and select only the correct chains
    return_dir=os.getcwd()    
    os.chdir(directory)
    pdb_name=pdb_loc.split('/')[-1] 
    
    if chain<>'All':
        Select_chain(pdb_name,chain)

    #Clean the pdb    
    Clean(pdb_name)

    #Convert the pdb to lammps input
    Pdb2Lammps(pdb_name,pdb_cname)
    
    #Write the lammps input for the frustratometer
    
    
    #Write the parameters for the AWSEM simulation
    

    #Recalculate the ssweight
    SSweight(pdb_name)

    #Run the frustratometer
    

    #return frustration_data
    

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Analyzes the frustration of the residues of a protein')
    parser.add_argument('PDB', type=str, nargs='*', help='PDB files to be analyzed')
    parser.add_argument('-c','--chain',nargs='*',default=['All'] ,action='store', help='Selects the chains for each PDB')
    parser.add_argument('-f','--frustration', choices=['singleresidue','configurational','mutational','All'], help='Sets the frustration type used in the analysis')
    parser.add_argument('--download',help='Downloads the PDB file if the file does not exists',action='store_true')
    parser.add_argument('--freq',nargs=1,help='Uses a custom aminoacid frequency from a csv table',action='store',default=None)
    parser.add_argument('--awsem',nargs=1,help='Uses custom awsem parameters as input' ,action='store',default=None)
    parser.add_argument('--awsem',nargs=1,help='Uses custom awsem parameters as input' ,action='store',default=None)
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

    Frustration(pdb,chain,args.frustration,args.download)

if __name__=='__main__':
    Frustration('1jge','ABC','mutational')

    
