import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs

# ---------------------------------------------------------- #
# CSV FILES OF FEM MESH AND RESULTS INTO FEM VISUALISATIONS
# SAMUEL DEMIAN
# ---------------------------------------------------------- #
# inputs :
#name = 'Trade-off - A'
name = 'truss2_60'
#name = 'cvs_test'
factor = 10
# ---------------------------------------------------------- #

def ImportCsv(name):
    nameorig = name + " - original.csv"
    nametruss = name + " - truss.csv"
    namedisp = name + " - displacements.csv"
    nodes = np.loadtxt(nameorig, delimiter=',')
    disp = np.loadtxt(namedisp, delimiter=',')
    truss = np.loadtxt(nametruss, delimiter=',')
    truss -= 1

    return nodes, disp, truss

def CreateDisp(nodes, disp, factor=1):
    deformed_nodes = nodes + (disp * factor)
    return deformed_nodes

def GenTruss(node, trusses):
    truss = trusses
    nodes = np.transpose(node)
    Trussgraph = []
    for i in truss:
        truss_i = [[nodes[0, int(i[0])], nodes[0, int(i[1])]], [nodes[1, int(i[0])], nodes[1, int(i[1])]]]
        Trussgraph.append(truss_i)

    return Trussgraph

def LenTruss(truss):
    truss = np.array(truss)
    x1, x2, y1, y2 = truss[0, 1], truss[0, 1], truss[1, 0], truss[1, 1]
    length = np.sqrt((x2-x1)**2+(y2-y1)**2)
    return length

def Graph(truss_graph, deformed_truss_graph):

    # computing deformations
    deformation = []
    for i in range(len(truss_graph)):

        truss_undeformed = truss_graph[i]
        truss_deformed = deformed_truss_graph[i]
        length_undef = LenTruss(truss_undeformed)
        length_def = LenTruss(truss_deformed)
        strain = length_def / length_undef
        deformation.append(strain)

    max_strain = 4.8e+5#np.amax(deformation)
    min_strain = -3.1e+5 #np.amin(deformation)
    Emod = 1.04e11
    stress = Emod * np.array(deformation)


    # drawing original truss
    for i in truss_graph:

        x = i[0]
        y = i[1]
        plt.plot(x, y, color='black', linewidth=1)

    # drawing deformed truss
    for i in range(len(deformed_truss_graph)):
        deformed_truss_i = deformed_truss_graph[i]

        strain_i = deformation[i]
        stress_i = stress[i]
        # percentage = (strain_i - min_strain)/(max_strain - min_strain)
        percentage = (stress_i - min_strain) / (max_strain - min_strain)
        perc_color = percentage
        print(perc_color)

        if perc_color >= max_strain:
            color_i = 'yellow'
        elif perc_color <= min_strain:
            color_i = 'yellow'
        else:
            color_i = 'green'#(perc_color, 0, 1-perc_color)

        x = deformed_truss_i[0]
        y = deformed_truss_i[1]

        plt.plot(x, y, color=color_i, linewidth=2)

    plt.axis('equal')
    plt.grid()
    plt.xlabel('m')
    plt.ylabel('m')
    plt.title('Truss, '+ name)
    plt.show()

def main(name, factor):
    nodes, disp, truss_pairs = ImportCsv(name)

    nodes_def = CreateDisp(nodes, disp, factor)
    truss_graph = GenTruss(nodes, truss_pairs)
    deformed_truss_graph = GenTruss(nodes_def, truss_pairs)
    Graph(truss_graph, deformed_truss_graph)

main(name, factor)
