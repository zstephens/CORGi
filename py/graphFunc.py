

def longestPath_DAG(g,s,w):
	reverse_graph = {v[0]:[] for v in g}
	for v in g:
		for v2 in v[1]:
			reverse_graph[v2].append(v[0])
	reverse_graph = {k:list(set(reverse_graph[k])) for k in reverse_graph.keys()}
	start = False
	dist  = {v[0]:0 for v in g}
	dist[s[0]] = 999999999
	w_mat = {v[0]:[] for v in g}
	for v in g[::-1]:
		if start:
			if len(reverse_graph[v[0]]):
				w_mat[v[0]] = sorted([(dist[u]+w[u][v[0]],u) for u in reverse_graph[v[0]]],reverse=True)
				dist[v[0]]  = w_mat[v[0]][0][0]
		if v[0] == s[0]:
			start = True
	maxV = max(dist.values())
	# traceback
	outPath = []
	if maxV > 0:
		maxVal = [v for v in dist.keys() if dist[v] == maxV][0]
		vvv = maxVal
		outPath.append(vvv)
		while True:
			if len(w_mat[vvv]):
				#print vvv,'-->',w_mat[vvv][0][1]
				vvv = w_mat[vvv][0][1]

				outPath.append(vvv)
			else:
				break
	return outPath[::-1]

def exhaustive_DAG(g,w,p):
	topPaths = []
	for start in g:
		topPaths.append([p[start[0]],[start[0]]])
		dag_path = longestPath_DAG(g,start,w)
		#print 'RAWR', start[0], dag_path
		if len(dag_path) >= 2:
			for i in xrange(len(dag_path)):
				for j in xrange(i+1,len(dag_path)+1):
					path    = dag_path[i:j]
					myPrize = p[path[0]]
					for k in xrange(1,len(path)):
						myPrize +=  w[path[k-1]][path[k]]
					topPaths.append([myPrize,[n for n in path]])
	#for n in sorted(topPaths,reverse=True):
	#	print 'asdf:', n[1], n[0]
	return sorted(topPaths,reverse=True)[0]

def dfs_scoreAllPaths(weighted_graph,node_prize):
	visited   = {}
	all_paths = []
	for node in node_prize.keys():
		# dfs all paths
		stack = [(node,[node])]
		paths = []
		while stack:
			(vertex,path) = stack.pop()
			toExplore     = [n for n in weighted_graph[vertex] if (n[0] not in path)]
			if not len(toExplore):
				paths.append([n for n in path])
			else:
				for node2 in toExplore:
					stack.append((node2,path+[node2]))
		all_paths.extend([n for n in paths])
	#print weighted_graph
	#print node_prize
	sorted_paths = []
	for path in all_paths:
		myScore = -len(path)
		for i in xrange(len(path)):
			myScore += node_prize[path[i]]
			if i > 0:
				if path[i] in weighted_graph[path[i-1]]:
					myScore += weighted_graph[path[i-1]][path[i]]
		sorted_paths.append([myScore,[n for n in path]])
	sorted_paths = sorted(sorted_paths,reverse=True)
	#for n in sorted_paths:
	#	print n
	return sorted_paths

# return all visited nodes given an adjaceny matrix and a starting index
def dfs_adj(A,start):
	stack   = [start]
	visited = set()
	while stack:
		vertex = stack.pop()
		visited.add(vertex)
		for i in xrange(len(A[vertex])):
			if i not in visited and A[vertex][i]:
				stack.append(i)
	return sorted(list(visited))

# return a list of all connected subgraphs
def get_connected_subgraphs(A):
	clusters_out = []
	not_visited = {i:True for i in xrange(len(A))}
	for k in not_visited.keys():
		if not_visited[k]:
			newCluster = dfs_adj(A,k)
			clusters_out.append([n for n in newCluster])
			for n in newCluster:
				not_visited[n] = False
	return clusters_out


if __name__ == '__main__':

	test_A = [[0,1,0,0,0],
	          [1,0,0,0,0],
	          [0,0,0,1,1],
	          [0,0,1,0,1],
	          [0,0,0,1,0]]

	print get_connected_subgraphs(test_A)

