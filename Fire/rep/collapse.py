import sys

# first argument is the filename of the input fa
input_fa = open(sys.argv[1])
for x in input_fa:
  x = x.strip()
  if len(x) > 0 and x[0] == '>':
    x_splitted = x.split('|')
    name = ''
    for y in x_splitted:
      if y.find('IG') == 0:
        name = y
        break
    print '>' + name
  else:
    print x.replace('.','').upper()
  
