# Project:  CSE 450 Code Project 2
# Filename: tcrdiversity.py
# Authors:  Kushagra Kshatri
 
"""
tcrdiversity: Analyze MST-based diversity metrics for T cell repertoires.
"""
 
import csv
import json
import os.path as osp
import numpy as np
 
 
# Import BLOSUM-62 matrix as a global variable.
with open('blosum62.json') as f:
    BLOSUM62 = json.load(f)
 
 
def sw_align(seq1, seq2, gap=-1.0):
    """
    Computes the Smith-Waterman local alignment between two amino acid sequences
    using a linear gap penalty and the BLOSUM-62 substitution matrix.
 
    :param seq1: a string amino acid sequence, possibly containing wildcards *
    :param seq2: a string amino acid sequence, possibly containing wildcards *
    :param gap: a float gap penalty per unmmatched amino acid
    :returns: the float Smith-Waterman alignment score for the two sequences
    :returns: a string representation of the resulting local alignment
    """
    
    nrows, ncols = len(seq1) + 1, len(seq2) + 1
    OPT = np.zeros((nrows, ncols))
    max_score, max_i, max_j = 0, 0, 0
 
    for i in range(1, nrows):
        for j in range(1, ncols):
            match = BLOSUM62[seq1[i - 1]][seq2[j - 1]]
            diagonal = OPT[i - 1, j - 1] + match
            up = OPT[i - 1, j] + gap
            left = OPT[i, j - 1] + gap
            OPT[i, j] = max(0, diagonal, up, left)
 
            if OPT[i, j] > max_score:
                max_score = OPT[i, j]
                max_i, max_j = i, j
 
    # Trace back the OPT matrix, starting from the maximum OPT score. If there
    # are multiple maximum OPT scores, choose the upmost, leftmost one. If there
    # are ties during the trace, prefer matching (up-left) over gaps and prefer
    # gaps in seq2 (up) over gaps in seq1 (left).
    
    aligned_seq1, aligned_seq2 = [], []
    i, j = max_i, max_j
    while OPT[i, j] > 0:
        if OPT[i, j] == OPT[i - 1, j - 1] + BLOSUM62[seq1[i - 1]][seq2[j - 1]]:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i, j = i - 1, j - 1
        elif OPT[i, j] == OPT[i - 1, j] + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
 
    aligned_seq1 = ''.join(aligned_seq1[::-1])
    aligned_seq2 = ''.join(aligned_seq2[::-1])
 
    return max_score, f"{aligned_seq1}\n{aligned_seq2}"
 
 
class Graph:
    """
    A TCR repertoire represented as an undirected graph with TCR nodes and edge
    weights corresponding to pairwise Smith-Waterment local alignment scores.
    """
 
    def __init__(self, id, gap=-1.0):
        """
        Initializes a new Graph with nodes, edges, and edge weights based on the
        input TCR sequences. An edge (i, j) exists between TCR sequences i and j
        if and only if their Smith-Waterman local alignment score is strictly
        positive. If so, the weight of this edge is 1 / sw_align(i, j, gap).
 
        :param id: a string identifier for the input sequence data
        :param gap: a float gap penalty per unmmatched amino acid
        :returns: an instantiated Graph object
        """
        input_file = osp.join("inputs", f"{id}.csv")
        with open(input_file, "r") as f:
            tcr_reader = csv.reader(f)
            self.tcr_sequences = [seq[0] for seq in tcr_reader]
 
        
        num_sequences = len(self.tcr_sequences)
        self.adj_matrix = np.zeros((num_sequences, num_sequences))
        self.num_nodes = len(self.tcr_sequences)
        
        for i in range(num_sequences):
            for j in range(i + 1, num_sequences):
                score, _ = sw_align(self.tcr_sequences[i], self.tcr_sequences[j], gap)
                if score > 0:
                    weight = 1 / score
                    self.adj_matrix[i, j] = weight
                    self.adj_matrix[j, i] = weight
 
 
    def edge_query(self, i, j):
        """
        Determines if the specified edge (i, j) exists in the graph and, if so,
        what its weight is.
 
        :param i: an int index of the i-th TCR sequence
        :param j: an int index of the j-th TCR sequence
        :returns: a float weight of edge (i, j) if it exists; None otherwise
        """
        if self.adj_matrix[i, j] != 0:
            return self.adj_matrix[i, j]
        else:
            return None
 
 
 
    def connected(self):
        """
        Determines if the graph is connected.
 
        :returns: True if the graph is connected; False otherwise.
        """
        def dfs(node, visited):
            visited[node] = True
            for neighbor, edge_weight in enumerate(self.adj_matrix[node]):
                if edge_weight != 0 and not visited[neighbor]:
                    dfs(neighbor, visited)
 
        visited = [False] * self.num_nodes
        dfs(0, visited)
        return all(visited)
 
 
def diversity_mst(graph):
    """
    Reports the diversity of a TCR repertoire by computing the total weight of
    its minimum spanning tree. Raises a ValueError if the graph is disconnected.
 
    :returns: the float total weight of the MST.
    """
    if not graph.connected():
        raise ValueError("Graph is disconnected")
    
    edges = []
    for i in range(graph.num_nodes):
        for j in range(i + 1, graph.num_nodes):
            weight = graph.edge_query(i, j)
            if weight is not None:
                edges.append((i, j, weight))
 
    edges.sort(key=lambda x: x[2])
 
    parent = [i for i in range(graph.num_nodes)]
    rank = [0] * graph.num_nodes
 
    mst_weight = 0
    e = 0
    while e < graph.num_nodes - 1:
        u, v, w = edges.pop(0)
        u_root = find(parent, u)
        v_root = find(parent, v)
 
        if u_root != v_root:
            e += 1
            mst_weight += w
            union(parent, rank, u_root, v_root)
 
    return mst_weight
 
def find(parent, i):
    if parent[i] == i:
        return i
    parent[i] = find(parent, parent[i])
    return parent[i]
 
def union(parent, rank, x, y):
    x_root = find(parent, x)
    y_root = find(parent, y)
    
    if rank[x_root] < rank[y_root]:
        parent[x_root] = y_root
    elif rank[x_root] > rank[y_root]:
        parent[y_root] = x_root
    else:
        parent[y_root] = x_root
        rank[x_root] += 1