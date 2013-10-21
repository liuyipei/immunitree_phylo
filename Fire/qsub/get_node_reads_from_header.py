import sys, math
if len(sys.argv) < 3:
  print 'Usage: \n python get_node_reads_from_header.py <infile> <node name>'
  exit(1)

infile_name = sys.argv[1]
target_node = sys.argv[2]
infile = open(infile_name)

for x in infile:
  [header, node, seq] = x.split(',')
  header = header.strip()
  node = node.strip()
  seq = seq .strip()
  if node == target_node:
    print '>' + header
    print seq
