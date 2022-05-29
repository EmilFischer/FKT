import networkx as nx

def halin(r, n):
    G = nx.full_rary_tree(r, n)
    leafs = [x for x in G.nodes() if len(list(G.neighbors(x))) == 1]
    for i in range(len(leafs)-1):
        j = i+1
        G.add_edge(leafs[i], leafs[j])
    return G