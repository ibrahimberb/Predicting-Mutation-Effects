# import urllib.parse
# import urllib.request
#
# url = 'https://www.uniprot.org/uploadlists/'
#
# params = {
#     'from': 'PDB_ID',
#     'to': 'ID',  # SWISSPROT, GENENAME
#     'format': 'tab',
#     'query': '1GRN 3BN9 3NPS 1A4Y 1Z7X 1A22 1BP3'
# }
#
# data = urllib.parse.urlencode(params)
# data = data.encode('utf-8')
# req = urllib.request.Request(url, data)
# with urllib.request.urlopen(req) as f:
#     response = f.read()
# print(response.decode('utf-8'))

import re

pattern = re.compile(r"^([A-Z])([A-Z])(\d+)([A-Z])$")

s = "LI45G"

print(pattern.match(s).groups())
original_aa, chain_identifier, residue_number, mutant_aa = pattern.match(s).groups()



