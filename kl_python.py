from kernighan_lin import *
import networkx as nx
from networkx.algorithms.cuts import cut_size
from scipy.io import mmread
import time

erdos = mmread('./matrixes/Erdos02/Erdos02.mtx')
graph_erdos = nx.Graph(erdos)
begin_erdos = time.time()
a_erdos, b_erdos = kernighan_lin_bisection(graph_erdos)
end_erdos = time.time()

cut_size_erdos = cut_size(graph_erdos, a_erdos, b_erdos)
print("Erdos:", end_erdos - begin_erdos, "s", cut_size_erdos)


com_DBLP = mmread('./matrixes/com-DBLP/com-DBLP.mtx')
graph_com_DBLP = nx.Graph(com_DBLP)
begin_com_DBLP = time.time()
a_com_DBLP , b_com_DBLP  = kernighan_lin_bisection(graph_com_DBLP)
end_com_DBLP = time.time()

cut_size_com_DBLP = cut_size(graph_com_DBLP, a_com_DBLP, b_com_DBLP)
print("com-DBLP:", end_com_DBLP - begin_com_DBLP, "s", cut_size_com_DBLP)


rgg_n_2_20_s0 = mmread('./matrixes/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx')
graph_rgg_n_2_20_s0 = nx.Graph(rgg_n_2_20_s0)
begin_rgg_n_2_20_s0 = time.time()
a_rgg_n_2_20_s0, b_rgg_n_2_20_s0 = kernighan_lin_bisection(graph_rgg_n_2_20_s0)
end_rgg_n_2_20_s0 = time.time()

cut_size_rgg_n_2_20_s0 = cut_size(graph_rgg_n_2_20_s0, a_rgg_n_2_20_s0, b_rgg_n_2_20_s0)
print("rgg_n_2_20_s0:", end_rgg_n_2_20_s0 - begin_rgg_n_2_20_s0, "s", cut_size_rgg_n_2_20_s0)
