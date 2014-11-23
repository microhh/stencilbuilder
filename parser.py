from StencilBuilder import *

f = file("test.cxx")
lines = f.readlines()
f.close()

record = False
blocks = []

lineindex = -1
block = []
names = []
indent = 0

for n in lines:
  lineindex += 1
  line = n.find("//$")

  # Process the keyword
  if(line != -1):
    indent = line
    kw = n[line + 3::].strip().split(' ')
    if(kw[0] == "SBStart"):
      if(len(kw) == 1):
        raise(RuntimeError)
      else:
        names[0:0] = kw[1::]
      startline = lineindex
      record = True
      continue

    elif(kw[0] == "SBEnd"):
      endline = lineindex
      record = False
      blocks.append((block, startline, endline, names, indent))
      block = []
      names = []
      continue

    else:
      raise(RuntimeError)

  if(record):
    block.append(n.strip()+'\n')

# Execute the blocks in reverse order
blocks.reverse()

for n in blocks:
  # Execute the block and save the output
  exec(''.join(n[0]))

  # Loop over the variables to be printed
  test = []
  for v in n[3]:
    #indent = ''.rjust(n[4])
    indent = len(v) + 10
    evalstring = v + ".getString(0,0,0,{0})".format(indent)
    tmplist = (v + "[i,j,k] = " + eval(evalstring)).split('\n')
    # Indent each line to make it fit properly aligned into the code
    for t in range(len(tmplist)):
      tmplist[t] = ''.rjust(n[4]) + tmplist[t]

    print(tmplist[0])
    test[-1:-1] = tmplist[:]


  # Delete the StencilBuilder lines
  del(lines[n[1]:n[2]+1])
  # Replace it with the new code
  lines[n[1]:n[1]] = test[:]

for n in lines:
  print(n.rstrip())
