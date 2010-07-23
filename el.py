"""
CREATE TABLE el_rsids (
  rsid VARCHAR(255) NOT NULL,
  position INT,
  alleles VARCHAR(3),
  cytogenetic_band VARCHAR(255),
  genes VARCHAR(255),
  role VARCHAR(255),
  aa_change VARCHAR(255),
  maf_ceu REAL,
  maf_controls REAL,
  log10_bf REAL,
  n_cases REAL,
  n_controls REAL,
  f_cases_AA REAL,
  f_cases_AB REAL,
  f_cases_BB REAL,
  f_controls_AA REAL,
  f_controls_AB REAL,
  f_controls_BB REAL,
  PRIMARY KEY (rsid)
)
"""

import sys
import os
import genotype_tools

import snp
import mysql

base_map = {
  "A": "T",
  "T": "A",
  "C": "G",
  "G": "C"
}

def get_snp_info(el_snp):
  attribute_map = {
    "aa": [12, 15],
    "ab": [13, 16],
    "bb": [14, 17],
  }

  snp_info = {}

  query = "SELECT * FROM diseases.el_snps WHERE rsid = %s;"# % el_snp.rsid
  db.execute(query, (el_snp.rsid,))
  res = db.fetchone()

  snp_info["alleles"] = res[2].split("/")
  for k, vs in attribute_map.items():
    snp_info[k] = [res[v] for v in vs]

  return snp_info

# Return (f_cases_*, f_controls_*) for this SNP.
def get_probabilities(el_snp, flipped=False):
  print(el_snp.genotype, flipped)
  # From their results.
  snp_info = get_snp_info(el_snp)

  if el_snp.genotype              == 2 * snp_info["alleles"][0]:
    return snp_info["aa"]
  elif sorted(el_snp.genotype)    == sorted(snp_info["alleles"]):
    return snp_info["ab"]
  elif el_snp.genotype            == 2 * snp_info["alleles"][1]:
    return snp_info["bb"]
  # Awww shit...
  else:
    if not flipped:
      el_snp.genotype = ''.join([base_map[a] for a in el_snp.genotype])
      return get_probabilities(el_snp, True)
    print('Hmmm')
    return (1, 1)

database = mysql.connector.Connect(host='marlowe', user='gene210-user', 
  passwd='genomics', buffered=True)
db = database.cursor()

#user_snps_path = sys.argv[1].strip()
user_snps_path = 'data/genome_stuart_kim.txt'
#el_rsids_path   = sys.argv[2].strip()
el_rsids_path   = 'data/EL_SNPs.csv'

population = "CEU"

if not os.path.exists(user_snps_path) or not os.path.exists(el_rsids_path):
  print("Missing a file - check that the filenames you gave exist.")
  sys.exit(2)

user_snps = genotype_tools.FileUtils.read_genotype_file(user_snps_path)
el_rsids   = [line.strip() for line in open(el_rsids_path).readlines()]

# Now, we want to compute 
#   a = prior(EL) * product(i = 1 to 150, probability(SNP-i | EL))
#   b = prior(AL) * product(i = 1 to 150, probability(SNP-i | AL))
#   a / (a + b)
el_probability = 1

# For each extreme longevity SNP, we get the value of that SNP from the input
# and adjust the running EL and AL scores according to those values from the
# paper.
for el_rsid in el_rsids:
  # If SNP in SNP_hash, then use that value.
  el_snp = user_snps.get(el_rsid, None)
  
  # Impute if we have to.
  if el_snp is None:
    el_snp = genotype_tools.impute_rsid_simple(user_snps, el_rsid, population)
  # Could not impute, or some other error occurred.
  if el_snp is None:
    continue
  
  # Update probabilities - multiply each by the probability of having this
  # genotype given each condition (EL or AL).
  # Avoiding underflow at each step...
  probabilities = get_probabilities(el_snp)
  el_probability *= (probabilities[0] / probabilities[1])
  #el_probability *= (probabilities[0] / (probabilities[0] + probabilities[1]))

print("Chance of living to 100: %f" % (el_probability / (1 - el_probability)))
