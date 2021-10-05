import re

# pattern = re.compile(r"GN=(\w)+ ")
import requests
from mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


# pattern = re.compile(r"GN=(\S+)(\s)")
#
# S = ">sp|Q9UPU5|UBP24_HUMAN Ubiquitin carboxyl-terminal hydrolase 24 OS=Homo sapiens OX=9606 GN=USP24 PE=1 SV=3"
#
# print(re.search(pattern, S))
# print(re.search(pattern, S).group(0))
# print(re.search(pattern, S).group(1))
# print(re.search(pattern, S).group(2))


def get_gene_from_fasta(fasta_text):
    info_line = fasta_text.split('\n')[0]

    # pattern = re.compile(r"GN=(.+)")
    pattern = re.compile(r"GN=(\S+)(\s)")

    # Does not exists in UniProt server.
    if re.search(pattern, info_line) is None:
        return "N/A"

    gene = re.search(pattern, info_line).group(1)

    return gene


def get_gene_id_from_uniprot(uniprot_id):
    log.debug("Retrieving sequence {} ...".format(uniprot_id))
    address = "http://www.uniprot.org/uniprot/{}.fasta".format(uniprot_id)
    n_attempt = 3
    attempt = 0
    while attempt < n_attempt:
        r = requests.get(address)
        if r.status_code == 200:
            gene = get_gene_from_fasta(r.text)
            return gene

        attempt += 1
        log.warning(attempt)

    log.critical(attempt)
    raise IOError("Cannot reach the results.")


while True:
    # print(get_gene_id_from_uniprot("P62136"))
    print(get_gene_id_from_uniprot("P62158"))
