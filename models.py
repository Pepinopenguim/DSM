
from typing import List, Tuple
import numpy as np

class Node(object):
    def __init__(self, coord:Tuple[float, float], support:Tuple[bool, bool]=(False, False), load:Tuple[float, float]=(0.0,0.0)):
        self.coord = coord
        self.x = coord[0]
        self.y = coord[1]
        self.support = support 
        self.load = load

class Element(object):
    _id_counter = 0
    def __init__(self, node1:Node, node2:Node, E:float|None=None, A:float|None=None):

        self.node1 = node1
        self.node2 = node2

        self.E, self._overwrite_E = (1.0, True) if E is None else (E, False)
        self.A, self._overwrite_A = (1.0, True) if A is None else (A, False)

        # define distance vector
        self.v = np.array(self.node2.coord) - np.array(self.node1.coord)

        # length
        self.L = np.linalg.norm(self.v)

        # unitary vector
        self.n = self.v / self.L

        # local stiffness
        # does not consider material or section (E, A)
        c, s = self.n[0], self.n[1]
        self.k_mod = [
            [ c**2,    s*c ,   -c**2,   -s*c   ],
            [ s*c ,    s**2,   -s*c ,   -s**2  ],
            [-c**2,   -s*c ,    c**2,    s*c   ],
            [-s*c ,   -s**2,    s*c ,    s**2  ],
        ] / self.L
    

class Truss(object):
    def __init__(self, elements: List[Element], E_global:float=1.0, A_global:float=1.0):
        self.elements = elements
        self.nodes = []
        for id, elt in enumerate(elements):
            elt.id = id

            for nd in (elt.node1, elt.node2):
                if nd not in self.nodes:
                    self.nodes.append(nd)

        
        self.num_nodes = len(self.nodes)
        
        # extremely important, will work as matrix indexes
        for id, nd in enumerate(self.nodes):
            nd.id = id

        # since its a 2d solver, each node creates 2 lines in the matrix
        self.K = np.zeros(shape=(self.num_nodes*2, self.num_nodes*2))

        # define supports for truss (as booleans)
        self.supports_array = np.array([
            b
            for nd in self.nodes
            for b in nd.support
        ])

        # define load values
        self.loads_array = np.array([
            l
            for nd in self.nodes
            for l in nd.load
        ])

        # overwrite values if with globals if applicable
        for elt in self.elements:
            elt.E = E_global if elt._overwrite_E else elt.E
            elt.A = A_global if elt._overwrite_A else elt.A

        self.calc_global_stiffness() # updates K

    def calc_global_stiffness(self):
        # Para cada elemento...
        for elt in self.elements:
            # multiplicar contribuições de seção e material
            # considerando E,A_global

            elt.k = elt.k_mod * elt.A * elt.E

            n1, n2 = elt.node1, elt.node2

            global_mapper = [
                2 * n1.id, 2 * n1.id + 1, 2 * n2.id, 2 * n2.id + 1
            ]

            for i in range(4):
                for j in range(4):
                    global_index = (global_mapper[i], global_mapper[j])
                    self.K[*global_index] += elt.k[i, j]

    def solve(self):
        self.K_sys = self.K.copy() # to not overwrite original
        self.F_sys = self.loads_array.copy()

        for i in range(self.num_nodes*2):
            # if is support
            if self.supports_array[i]:
                # we are subbing the lines with reactions,
                # with boundary conditions
                # this makes the linear system solvable
                self.K_sys[i, :] = np.zeros(self.num_nodes*2)
                self.K_sys[:, i] = np.zeros(self.num_nodes*2)
                self.K_sys[i, i] = 1
                self.F_sys[i] = 0

        self.u = np.linalg.solve(self.K_sys, self.F_sys)
        self.resultants = np.linalg.multi_dot([self.K, self.u])
        # with these values, solve for reactions
        





        





    