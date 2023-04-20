import numpy as np
import random
def cliqueGenerator(E, dim):
    class cliqueGenClass:
        def __init__(self, Ematrix):
            self.Ematrix = Ematrix
            self.size = Ematrix.shape[0]
            self.max = 0
            self.moplex = []
            self.labels = [""] * self.size
            self.cliques = []


        def findAdjacent(self, node):
            adjacents = []
            for i in range(self.size):
                if self.Ematrix[node][i] == 1:
                    adjacents.append(i)
            return adjacents
            


        def LEXBFS(self):
            stop = False
            counter = self.size-1
            while(not stop):
                if counter == self.size-1:
                    node = random.choice(range(0, self.size))
                adjacents = self.findAdjacent(node)
                self.labels[node] = '0'
                if counter == self.size-1:
                    self.cliques.append([node])
                else:
                    self.cliques[-1].append(node)
                if self.labels.count('0') == self.size:
                    stop = True
                for j in adjacents:
                    if self.labels[j] != '0':
                        self.labels[j]= self.labels[j] + str(counter+1)
                currentMax = int(max(self.labels))
                if currentMax <= self.max:
                    self.cliques.append([])
                self.max = currentMax
                index_value = [index for index in range(len(self.labels)) if self.labels[index] == str(currentMax)]
                node = random.choice(index_value)
                counter = counter - 1
            
            return self.cliques[:-1]

        def cliquesGen(self):
            order = self.LEXBFS()
            alpha = np.zeros(self.size)
            cliques = [[]]
            dummy = []
            for element in order:
                for item in element:
                    dummy.append(item)
            for element in order:
                adj = self.findAdjacent(element[-1])
                for i in adj:
                    if dummy.index(i) <= dummy.index(element[-1]):
                        cliques[-1].append(i)
                
                cliques.append([])
            #print(cliques[:-1])
            return cliques[:-1]


    Ec = []
    cliqueGenClass = cliqueGenClass(np.array(E))
    maxCliques = cliqueGenClass.cliquesGen()
    # print(maxCliques)
    # for clique in maxCliques:
    #     dummy = np.zeros((np.shape(clique)[0], dim))
    #     print('dummy')
    #     for i in range(np.shape(dummy)[0]):
    #         dummy[i][clique[i]] = 1
    #     Ec.append(dummy)
    return maxCliques


res = cliqueGenerator(E, dim)