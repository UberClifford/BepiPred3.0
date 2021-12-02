from Bio.PDB import *

parser = PDBParser()
s = parser.get_structure("1ob1", "1ob1.pdb")
io = PDBIO()
keepAltID = "A"

#class NMROutputSelector2(Select):  # Inherit methods from Select class
#    def accept_atom(self, atom):
#        if (not atom.is_disordered()) or atom.get_altloc() == keepAltID:
#            atom.set_altloc(" ")  # Eliminate alt location ID before output.
#            return True
#        else:  # Alt location was not one to be output.
#            return False

class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"

io = PDBIO()
io.set_structure(s)
#io.save("1ob1_ordered.pdb", select=NMROutputSelector2())
io.save("1ob1_ordered.pdb", select=NotDisordered())
