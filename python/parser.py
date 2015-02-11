#!/usr/bin/python
import sys
import StringIO

from StencilBuilder import *

f = file("test.cxx", "r")
lines = f.readlines()
f.close()

# Remove the line breaks
for n in range(len(lines)):
  lines[n] = lines[n].rstrip()

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
      startline = lineindex
      record = True
      continue

    elif(kw[0] == "SBEnd"):
      endline = lineindex
      record = False
      blocks.append((block, startline, endline, indent))
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
  # Get the code block.
  code = ''.join(n[0])

  # Execute the code block and capture the output.
  buffer = StringIO.StringIO()
  # Temporarily redirect stdout to buffer, restore after call!
  sys.stdout = buffer
  exec(code)
  sys.stdout = sys.__stdout__
  output = buffer.getvalue().splitlines()
  del(code)
  del(buffer)

  # Delete the StencilBuilder lines.
  del(lines[n[1]:n[2]+1])

  # Replace it with the new code.
  lines[n[1]:n[1]] = output

f = file("test_new.cxx", "w")
for n in lines:
  f.write("{0}\n".format(n))
f.close()
