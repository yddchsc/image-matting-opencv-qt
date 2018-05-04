#!/usr/bin/env python
#-*- coding:utf8 -*-
import cv2

class Graph(object):

    def __init__(self, *args, **kwargs):
        self.node_neighbors = {}
        self.visited = {}
        self.bgdArea = []
        self.fgdArea = []
        self.middleArea = []

    def add_nodes(self, nodelist):
        for node in nodelist:
            self.add_node(node)

    def add_node(self, node):
        if not node in self.node_neighbors:
            self.node_neighbors[node] = []

    def add_edge(self, edge):
        u, v = edge
        if (v not in self.node_neighbors[u]) and (u not in self.node_neighbors[v]):
            self.node_neighbors[u].append(v)
            if (u != v):
                self.node_neighbors[v].append(u)

    def nodes(self):
        return self.node_neighbors.keys()

    def depth_first_search(self, root=None):
        order = []
        def dfs(node):
            self.visited[node] = True
            order.append(node)
            for n in self.node_neighbors[node]:
                if not n in  self.visited:
                    dfs(n)
        if root:
            dfs(root)
        for node in self.nodes():
            if not node in self.visited:
                dfs(node)
        #print order
        return order

    def breadth_first_search(self, root=None, hists=[]):
        queue = []
        order = []
        def bfs():
            while len(queue) > 0:
                node = queue.pop(0)
                self.visited[node] = True
                for n in self.node_neighbors[node]:
                    if (not n in self.visited) and (not n in queue):
                        for neighbour in self.node_neighbors[n]:
                            #print(cv2.compareHist(hists[n],hists[neighbour],0))
                            #print('neighbour' + str(neighbour))
                            if cv2.compareHist(hists[n],hists[neighbour],0)>0.5 and (neighbour in self.middleArea):
                                if n in self.fgdArea:    
                                    self.middleArea.remove(neighbour)
                                    self.fgdArea.append(neighbour)
                                elif n in self.bgdArea:
                                    self.middleArea.remove(neighbour)
                                    self.bgdArea.append(neighbour)
                        queue.append(n)
                        order.append(n)
        if root:
            queue.append(root)
            order.append(root)
            bfs()
        for node in self.nodes():
            if not node in self.visited:
                queue.append(node)
                order.append(node)
                bfs()
        #print order
        return order


# if __name__ == '__main__':
#     g = Graph()
#     g.add_nodes([i+1 for i in range(8)])
#     g.add_edge((1, 2))
#     g.add_edge((1, 3))
#     g.add_edge((2, 4))
#     g.add_edge((2, 5))
#     g.add_edge((4, 8))
#     g.add_edge((5, 8))
#     g.add_edge((3, 6))
#     g.add_edge((3, 7))
#     g.add_edge((6, 7))
#     print "nodes:", g.nodes()
#     order = g.depth_first_search(1)
#     order = g.breadth_first_search(1)