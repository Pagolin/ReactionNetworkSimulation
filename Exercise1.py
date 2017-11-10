
import matplotlib.pyplot as plt
import networkx as nx

from ReactionNetwork import SimpleReactionNetwork as sRN

species = ["Ec","Em","Rt","M","Rd","RDIm","I","RDIc"]
species2 = ["A","B","C","D"]
inputs = [[1,0,0,0],[0,1,1,0],[0,0,0,1]]
outputs = [[0,2,0,0],[0,0,0,1],[1,0,0,0]]
'''
newNetwork = sRN(species2, inputs, outputs)

listIn= [('C',2)]
listOut=[('A',1)]
newNetwork.add_reaction((listIn, listOut))
newNetwork.get_degree()
newNetwork.calculate_Cons_Law()
'''
yeastNetwork = sRN(species)
yeastNetwork.add_reaction(([("Ec",1),("Rd",1)], [("Em",1)]))
yeastNetwork.add_reaction(([("Em",1)], [("Ec",1)]))
yeastNetwork.add_reaction(([("Em",1), ("Rt", 1)], [("M",1)]))
yeastNetwork.add_reaction(([("M",1)],[("Em",1), ("Rt", 1)]))
yeastNetwork.add_reaction(([("Em",1), ("Rd", 1)], [("Em",1), ("Rt", 1)]))
yeastNetwork.add_reaction(([("M",1),("Rd",1)], [("Rt", 1), ("M",1)]))
yeastNetwork.add_reaction(([("Rt",1)], [("Rd",1)]))
yeastNetwork.add_reaction(([("Ec",1), ("Rt", 1)], [("M",1)]))
yeastNetwork.add_reaction(([("Rd",1), ("I", 1)], [("RDIm",1)]))
yeastNetwork.add_reaction(([("RDIm",1)], [("Rd",1), ("I", 1)]))
yeastNetwork.add_reaction(([("RDIm",1)], [("RDIc",1)]))
yeastNetwork.add_reaction(([("RDIc",1)], [("RDIm",1)]))

yeastNetwork.plot_dependencies("yeastNetworkGraphs")