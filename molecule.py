TO_ANG_CUBED = 1.6605390404
import numpy as np
from data import atomic_masses

class Molecule:

    def __init__(self, name, coords, atoms, density=-1.0):
        self.name = name 
        self.coords = coords
        self.atoms = atoms
        self.density = density 

        if (self.density < 0):
            v = 4.0 * np.pi * self.radius_from_coords() ** 3 / 3.0
            self.density = self.mass() / v
    
    #total molar mass of molecule in g/mol
    def mass(self):
        total = 0
        for atom in self.atoms:
            total += atomic_masses[atom]
        return total

    #determining radius as longest distance between atoms in molecule in angstrom 
    def radius_from_coords(self):
        def separation(x1,x2):
            return np.linalg.norm(x1 - x2)
            
        n = self.coords.shape[0]
        c = self.coords
        r = np.array([ [separation(c[i], c[j]) for j in range(n)] for i in range(n) ])        

        return np.amax(r)

    #volume of molecule in angstrom cubed   
    def volume(self):
        return ( self.mass() * TO_ANG_CUBED / self.density )

    #radius in angstrom
    def radius(self):
        return ( 0.75 * self.volume() / np.pi ) ** (1.0/3.0)
        
    def pretty_print(self):
        print("Name:", self.name)
        print("Density:", self.density)
        print("Mass:", self.mass() )
        print("Radius:", self.radius() )





if __name__ == '__main__':
    m1 = Molecule('h2o', np.array([[1,1,1], [0, 0, 0], [-1, -1, -1]]),['h','o','h'],)
    m2 = Molecule('hcl', np.array([[0, 0, 0], [1.8, 0, 0]]),['h','cl'], 1)
    molecules = [m1, m2]

    for m in molecules:
        m.pretty_print()
