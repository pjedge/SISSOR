# -*- coding: utf-8 -*-

"""
THIS CODE TAKEN FROM http://code.activestate.com/recipes/577225-union-find/,
which in turn is a direct port of the wikipedia pseudocode for disjoint-set-forest
All credit to the original author, Ahmed El Deeb


MakeSet(x) initializes disjoint set for object x
Find(x) returns representative object of the set containing x
Union(x,y) makes two sets containing x and y respectively into one set

Some Applications:
- Kruskal's algorithm for finding minimal spanning trees
- Finding connected components in graphs
- Finding connected components in images (binary)
"""

#def MakeSet(x):
#     x.parent = x
#     x.rank   = 0

def Union(x, y):
     xRoot = Find(x)
     yRoot = Find(y)
     if xRoot.rank > yRoot.rank:
         yRoot.parent = xRoot
     elif xRoot.rank < yRoot.rank:
         xRoot.parent = yRoot
     elif xRoot != yRoot: # Unless x and y are already in same set, merge them
         yRoot.parent = xRoot
         xRoot.rank = xRoot.rank + 1

def Find(x):
     if x.parent == x:
        return x
     else:
        x.parent = Find(x.parent)
        return x.parent
        
class Node:
    def __init__ (self, label):
        self.label = label
        self.parent = self
        self.rank   = 0
    def __str__(self):
        return self.label