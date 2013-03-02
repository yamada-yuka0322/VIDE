import shutil
import tempfile
import sys
import re

rex = "@FILENAME@"

filename = sys.argv[1]


fh = file("header.txt")
header = fh.read()
header_translated = re.sub(r'@FILENAME@', filename, header)
fh.close()

f = file(filename)
lines = f.read()
f.close()

lines = re.sub(r'(?s)/\*\+.*\+\*/','',lines)
lines = header_translated + lines

with tempfile.NamedTemporaryFile(delete=False) as tmp_sources:
  tmp_sources.write(lines)

shutil.move(tmp_sources.name, filename)
