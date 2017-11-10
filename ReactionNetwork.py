import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from sympy import Matrix


class SimpleReactionNetwork:

    def __init__(self, species, inputs = None, outputs = None):
        self.species = species
        self.input_array = np.array(inputs).transpose() if inputs != None else np.zeros((len(species),1))
        print("InputArray", self.input_array)
        self.output_array = np.array(outputs).transpose() if outputs != None else np.zeros((len(species),1))
        self.total = np.subtract(self.output_array, self.input_array)
        # A reaction r1 in row i influences a reaction r2 in column j
        # if any species from r1 occurs as in the input matrix row for r2
        # => matrix product at (i,j != 0)
        self.reactionCoupling = np.dot(np.add(self.input_array, self.output_array).T, self.input_array) > 0

    def add_reaction(self, reaction):
        in_tuples, out_tuples = reaction
        print("In tuples", in_tuples)
        newColumn = np.zeros((len(self.species),1))
        for (spec, count) in in_tuples:
            pos = self.species.index(spec)
            newColumn[pos,0] = count
        self.input_array = np.c_[self.input_array, newColumn]

        newColumn = np.zeros((len(self.species), 1))
        for (spec, count) in out_tuples:
            pos = self.species.index(spec)
            newColumn[pos, 0] = count
        self.output_array = np.c_[self.output_array, newColumn]
        #Update
        self.total = np.subtract(self.output_array, self.input_array)
        self.reactionCoupling = np.dot(np.add(self.input_array, self.output_array).T, self.input_array) > 0
        print("New Input", self.input_array, "New Output", self.output_array)

    def plot_dependencies(self, name):
        plt.figure()
        plt.subplot(221)
        plt.axis('off')
        plt.title("Reactions-Dependency")
        nx.draw_networkx(self.get_reactions_dependency_graph(), node_color="blue", edge_color="darkgrey", nodesize=3)
        plt.subplot(222)
        plt.axis('off')
        plt.title("Species-Dependency")
        nx.draw_networkx(self.get_species_dependency_graph(), node_color="blue", edge_color="darkgrey", nodesize=3)
        plt.show()
        plt.savefig("./graphicOutput/"+str(name))

    def calculate_total(self):
        return np.subtract(self.output_array, self.input_array)


    def get_reaction_coupling(self):
        allInvolvedSpecies = np.add(self.input_array, self.output_array)
        reactionCouplingMatrix = np.dot(allInvolvedSpecies.transpose(),self.input_array) > 0

    def get_degree(self):
        print( np.sum(self.reactionCoupling, axis=0).max())

    def calculate_Cons_Law(self):
        totalMatrix = Matrix(self.total.T)
        nullSpace = np.array(totalMatrix.nullspace()[0]).flatten()
        nullSpace = ((1/nullSpace.min())*nullSpace).astype(int)
        consLaw = [((str(nullSpace[i])+"*"+ self.species[i])) for i in range(nullSpace.shape[0])]
        print(consLaw)


    def get_species_dependency_graph(self):
        speciesCouplingMatrix = np.dot(self.input_array, self.output_array.T)
        edges = []
        rows, columns = speciesCouplingMatrix.shape
        for i in range(rows):
            print(speciesCouplingMatrix[i], i)
            for j in range(columns):
                if speciesCouplingMatrix[i,j] > 0:
                    edges.append((self.species[i], self.species[j]))
        print("Species Edges", edges)
        graph = nx.DiGraph()
        graph.add_nodes_from(self.species)
        graph.add_edges_from(edges)
        return graph

    def get_reactions_dependency_graph(self):
        edges = []
        rows, columns = self.reactionCoupling.shape
        for i in range(rows):
            print("ReactionCouplingMatrix", self.reactionCoupling[i], i)
            for j in range(columns):
                if self.reactionCoupling[i,j]:
                    edges.append(("R"+str(i), "R"+str(j)))
        print("Reactions Edges",edges)
        graph = nx.DiGraph()
        graph.add_edges_from(edges)
        return graph