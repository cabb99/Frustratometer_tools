if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Analyzes the frustration of the residues of a protein')
    parser.add_argument('PDB', type=str, nargs='*', help='PDB files to be analyzed')
    parser.add_argument('-c','--chain', dest='accumulate', action='store_const',
                       const=sum, default=max, help='sum the integers (default: find the max)')
    parser.add_argument('-f','-frustration-type', choices=['singleresidue','configurational','mutational'], help='Sets the frustration type used in the analysis')
    parser.add_argument('--download',help='Downloads the PDB file if the file does not exists')
    parser.add_argument('args', nargs=argparse.REMAINDER)

    args = parser.parse_args()
    print args
