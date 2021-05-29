# Python program for implementation of Ford Fulkerson algorithm 
# This code is contributed by Neelam Yadav 
from collections import defaultdict 

#This class represents a directed graph using adjacency matrix representation 
class Graph: 

  def __init__(self,graph): 
    self.graph = graph # residual graph 
    self.ROW = len(graph) 

  '''Returns true if there is a path from source 's' to sink 't' in 
  residual graph. Also fills parent[] to store the path '''
  def BFS(self,s, t, parent): 

    # Mark all the vertices as not visited 
    visited =[False]*(self.ROW) 
    
    # Create a queue for BFS 
    queue=[] 
    
    # Mark the source node as visited and enqueue it 
    queue.append(s) 
    visited[s] = True
    
    # Standard BFS Loop 
    while queue: 

      #Dequeue a vertex from queue and print it 
      u = queue.pop(0) 
    
      # Get all adjacent vertices of the dequeued vertex u 
      # If a adjacent has not been visited, then mark it 
      # visited and enqueue it 
      for ind, val in enumerate(self.graph[u]): 
        if visited[ind] == False and val > 0 : 
          queue.append(ind) 
          visited[ind] = True
          parent[ind] = u 

    # If we reached sink in BFS starting from source, then return 
    # true, else false 
    return True if visited[t] else False
      
  
  # Returns tne maximum flow from s to t in the given graph 
  def FordFulkerson(self, source, sink): 

    # This array is filled by BFS and to store path 
    parent = [-1]*(self.ROW) 

    max_flow = 0 # There is no flow initially 

    # Augment the flow while there is path from source to sink 
    while self.BFS(source, sink, parent) : 

      # Find minimum residual capacity of the edges along the 
      # path filled by BFS. Or we can say find the maximum flow 
      # through the path found. 
      path_flow = float("Inf") 
      s = sink 
      while(s != source): 
        path_flow = min (path_flow, self.graph[parent[s]][s]) 
        s = parent[s] 

      # Add path flow to overall flow 
      max_flow += path_flow 

      # update residual capacities of the edges and reverse edges 
      # along the path 
      v = sink 
      while(v != source): 
        u = parent[v] 
        self.graph[u][v] -= path_flow 
        self.graph[v][u] += path_flow 
        v = parent[v] 

    return max_flow 


def redistribution(a, b, c, x, y, z):
  if (a + b + c) != (x + y + z):
    raise Exception('MAX FLOW - Invalid parameters.')
  T = a + b + c
  graph = [[0, a, b, c, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 0, 0, T, T, 0, 0, 0, 0, 0],
       [0, 0, 0, 0, T, 0, T, 0, 0, 0, 0],
       [0, 0, 0, 0, 0, T, T, 0, 0, 0, 0],
       [0, 0, 0, 0, 0, 0, 0, T, T, 0, 0],
       [0, 0, 0, 0, 0, 0, 0, T, 0, T, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, T, T, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, y],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, z],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      ]
  g = Graph(graph) 
  source = 0
  sink = 10
  g.FordFulkerson(source, sink)
  return [(g.graph[4][1], g.graph[6][2], g.graph[5][3]), (g.graph[7][4], g.graph[9][5], g.graph[8][6])]
