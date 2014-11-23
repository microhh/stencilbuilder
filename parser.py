from cStringIO import StringIO
from StencilBuilder import *
import sys

f = file("test.cxx")
lines = f.readlines()
f.close()

record = False
block  = []
linenumber = 0
for n in lines:
  linenumber += 1
  n = n.strip()
  line = n.find("//$")
  if(line != -1):
    kw = n[line + 3::].strip()
    print("KEYWORD: {0}".format(kw))
    if(kw == "SBStart"):
      startline = linenumber
      record = True
      continue

    elif(kw == "SBEnd"):
      endline = linenumber
      record = False

      # Execute the block
      blockstr = ""
      for n in block:
        print("DSL: {0}".format(n))
        blockstr += n + "\n"
      old_stdout = sys.stdout
      tmp = sys.stdout = StringIO()
      exec(blockstr)
      sys.stdout = old_stdout
      blockout = tmp.getvalue()
      continue

    else:
      raise(RuntimeError)

  if(record):
    block.append(n)

print(blockout)
print("DONE")
