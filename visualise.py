import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


G=nx.Graph()
H=nx.Graph()
nodes=np.loadtxt('nodedata.dat',usecols=[0])
nodepos=np.loadtxt('nodedata.dat',usecols=[1,2])
print (np.shape(nodes))
n_vtx=np.shape(nodes)[0]
for i in range(n_vtx) :
    G.add_node(nodes[i],pos=(nodepos[i,0],nodepos[i,1]))
    H.add_node(nodes[i],pos=(nodepos[i,0],nodepos[i,1]))
edges=np.loadtxt('edgedata.dat')
#medges=np.loadtxt('matched_edges.dat')
nedges=np.shape(edges)[0]
for i in range(nedges) :
    G.add_edge(edges[i,0],edges[i,1])
#nmedges=np.shape(medges)[0]
#for i in range(nmedges) :
#    H.add_edge(medges[i,0],medges[i,1])
pos=nx.get_node_attributes(G,'pos')
#nodes1=list(np.loadtxt('./regions/clust0.dat',usecols=[0],dtype=int))
#nodes2=list(np.loadtxt('./regions/clust4.dat',usecols=[0],dtype=int))
#nodes4=list(np.loadtxt('./regions/clust1.dat',usecols=[0],dtype=int))
#nodes3=list(np.loadtxt('./regions/clust3.dat',usecols=[0],dtype=int))
#H=G.subgraph(nodes1)
#nx.draw(H,pos,node_color='r',node_size=10)
#H=G.subgraph(nodes2)
#nx.draw(H,pos,node_color='r',node_size=10)
#H=G.subgraph(nodes3)
#nx.draw(H,pos,node_color='b',node_size=10)
#H=G.subgraph(nodes4)
#nx.draw(H,pos,node_color='g',node_size=10,with_labels=False)
#sG=G.subgraph(nodes1,)
#sH=H.subgraph(nodes1,)
#nx.draw(H,pos,node_color='y',with_labels=True,font_size=3,node_size=0,width=0.2)
nx.draw(G,pos,node_color='y',font_size=3,node_size=0,width=0.2)
#nx.draw(sH,pos,node_color='y',edge_color='b',font_size=3,node_size=0,width=3.2)
plt.savefig("test.pdf",format='pdf')
#nx.draw_networkxs(G,pos,node_color='b',nodelist=list(nodes1),node_size=10)
#nx.draw_networkx_nodes(G,pos,node_color='r',nodelist=list(nodes2),node_size=10)
#nx.draw_networkx_nodes(G,pos,node_color='y',nodelist=list(nodes3),node_size=10)
#nx.draw_networkx_nodes(G,pos,node_color='g',nodelist=list(nodes4),node_size=10)
