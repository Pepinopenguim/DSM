from models import Node, Element, Truss

# define nodes
n1 = Node((0,0), support=(True, True)) # True means limited!
n2 = Node((1,0), load=(200, 150)) # loads in Newton! Always SI
n3 = Node((0,1), support=(True, False))
n4 = Node((1,1), load=(0,-100))

# define elements of Truss
# mesh
el1 = Element(n1, n2)
el2 = Element(n1, n3)
el3 = Element(n1, n4)
el4 = Element(n2, n4)
el5 = Element(n3, n4)

# create truss
truss = Truss(
    elements=[el1, el2, el3, el4, el5], # nodes parsed automatically
    A_global=1,
    E_global=1 # global only overwrites elements not created with custom values
)

truss.calc_global_stiffness()

truss.solve()

print(truss.resultants)
