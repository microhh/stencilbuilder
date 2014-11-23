from StencilBuilder import *

f = file("test.cxx")
lines = f.readlines()
f.close()

record = False
blocks = []

lineindex = -1
block = []
for n in lines:
  lineindex += 1
  line = n.find("//$")

  # Process the keyword
  if(line != -1):
    kw = n[line + 3::].strip()
    #print("KEYWORD: {0}".format(kw))
    if(kw == "SBStart"):
      startline = lineindex
      record = True
      continue

    elif(kw == "SBEnd"):
      endline = lineindex
      record = False
      blocks.append((block, startline, endline))
      block = []
      continue

    else:
      raise(RuntimeError)

  if(record):
    block.append(n)

# Execute the blocks in reverse order
blocks.reverse()

for n in blocks:
  # Delete the StencilBuilder lines
  del(lines[n[1]:n[2]+1])
  # Replace it with the new code
  lines[n[1]:n[1]] = n[0]

for n in lines:
  print(n.rstrip())
