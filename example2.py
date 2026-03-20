from models import Node, Element, Truss

# define nodes
n1 = Node((0,0), support=(True, True))
n2 = Node((1,0), load=(0,-20))
n3 = Node((2,0), load=(0,-20))
n4 = Node((3,0), load=(0,-20))
n5 = Node((4,0), load=(0,-20))
n6 = Node((5,0), load=(0,-20))
n7 = Node((6,0), load=(0,-20))
n8 = Node((7,0), load=(0,-20))
n9 = Node((8,0), support=(True, True))
n10 = Node.as_polar((1,0), radius=.5, angle=90)
n11 = Node.as_polar((3,0), radius=.5, angle=90)
n12 = Node.as_polar((5,0), radius=.5, angle=90)
n13 = Node.as_polar((7,0), radius=.5, angle=90)
n14 = Node.as_polar((2,0), radius=1, angle=90)
n15 = Node.as_polar((6,0), radius=1, angle=90)
n16 = Node.as_polar((4,0), radius=2, angle=90)



# define elements of Truss
# mesh


# create truss
truss = Truss(
    elements=[
        Element(n1 , n2),
        Element(n2 , n3),
        Element(n3 , n4),
        Element(n4 , n5),
        Element(n5 , n6),
        Element(n6 , n7),
        Element(n7 , n8),
        Element(n8 , n9),
        Element(n1 , n10),
        Element(n2 , n10),
        Element(n3 , n10),
        Element(n3 , n14),
        Element(n3 , n11),
        Element(n10, n14),
        Element(n4, n11),
        Element(n11, n14),
        Element(n14, n16),
        Element(n5 , n11),
        Element(n5 , n16),
        Element(n5 , n12),
        Element(n6, n12),
        Element(n7, n12),
        Element(n7, n15),
        Element(n7, n13),
        Element(n8, n13),
        Element(n9, n13),
        Element(n13, n15),
        Element(n15, n16),
    ], 
)

truss.calc_global_stiffness()

truss.solve()

print(truss.resultants)

truss.draw()
