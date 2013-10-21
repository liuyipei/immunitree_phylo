import sys, math
infile_name = sys.argv[1]
print "opening %s" % infile_name
infile = open(infile_name)
nodes_per_seq_dict = dict()
for x in infile:
  [node, parent, nodesize, subtreesize, description, seq] = x.split(',')
  node = node.strip()
  seq = seq.strip()
  if not seq in nodes_per_seq_dict:
    nodes_per_seq_dict[seq] = set()
  nodes_per_seq_dict[seq].add(node)

for seq in nodes_per_seq_dict:
  nodes = nodes_per_seq_dict[seq]
  if len(nodes) > 1:
    print(seq)
    print nodes
print 'Searched for multiple nodes having the same concensus: done!'

