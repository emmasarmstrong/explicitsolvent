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

def unit_readxyz():
    X, atoms = readxyz("test.xyz")
    if len(X) == 102:
        print("Readxyz success.")
    else:
        print("Readxyz failed.")

if __name__ == "__main__":
    unit_readxyz()

