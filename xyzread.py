from argparse import ArgumentParser

def readxyz(filename):
    f=open(filename)
    lines=f.readlines()[2:]
    X=[]
    atoms=[]
    for l in lines:
        atom,x,y,z=l.split()
        X.append([float(x),float(y),float(z)])
        atoms.append(atom.lower())
    f.close()
    return X, atoms

def numberofatoms(filename):
    parser = ArgumentParser()
    parser.add_argument("number", type=int, help = "number of atoms in the first molecule")
    args = parser.parse_args()
    N = args.number
    return N

def unit_readxyz():
    X, atoms = readxyz("furan-if-1.xyz")
    N = numberofatoms("furan-if-1.xyz")
    mol1 = X[:N]
    if len(mol1) == 9:
        print("Readxyz success.")
    else:
        print("Readxyz failed.")

if __name__ == "__main__":
    unit_readxyz()

