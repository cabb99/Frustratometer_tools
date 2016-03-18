Temp_directory
import subprocess
Working_directory=
Cleaner='/home/carlos/Programas/mmtsb/perl/convpdb.pl'

def MakeTemp(pdb)

def Clean_pdb(pdb):
    #-nohetero ignores HETATM recodrd
    #-renumber renumbers all the residues from 1
    #-out Specifies the output format, in this case generic
    #-fixcoo fixes C-terminal atoms  
    commands=[Cleaner,'-nohetero','-renumber','1','-out','generic','-fixcoo','$1_$2.pdb']
    subprocess
    $convpdbscript -nohetero -renumber 1 -out generic -fixcoo $1_$2.pdb > $1_cleaned.pdb
    subprocess.call    
    return cleaned_pdb





def Frustration(pdb,chain,frustration,download,custom_aa_freq):
    print "Calculating frustration"
    pdb='1jge'
    chain='ABC'
    frustration='mutational'
    download=True
    #Make a directory
    MakeTemp(
    #Try to download the pdb
    
    #Open the pdb and select only the correct chains

    #Clean the pdb    

    #Convert the pdb to lammps input
    
    #Write the lammps input for the frustratometer
    
    #Write the parameters for the AWSEM simulation

    #Recalculate the ssweight

    #Run the frustratometer

    return frustration data



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Analyzes the frustration of the residues of a protein')
    parser.add_argument('PDB', type=str, nargs='*', help='PDB files to be analyzed')
    parser.add_argument('-c','--chain',nargs='*',default=['All'] ,action='store', help='Selects the chains for each PDB')
    parser.add_argument('-f','--frustration', choices=['singleresidue','configurational','mutational','All'], help='Sets the frustration type used in the analysis')
    parser.add_argument('--download',help='Downloads the PDB file if the file does not exists',action='store_true')
    parser.add_argument('--freq',nargs=1,help='Uses a custom aminoacid frequency from a csv table',action='store',default=None)
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
