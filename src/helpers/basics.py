import re

match = re.compile(r"^TCGA-(\w\w)-(\w\w\w\w)$")

print(match.match("TCGA-D8-A1XY"))
print(match.match("TCGA-D8-A1sdfsXY"))
