"""
CREATE TABLE el_rsids (
  db_snp VARCHAR(255) NOT NULL,
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
  f_cases_aa REAL,
  f_cases_ab REAL,
  f_cases_bb REAL,
  f_controls_aa REAL,
  f_controls_ab REAL,
  f_controls_bb REAL,
  PRIMARY KEY (rsid)
)
"""

import sys
import os
import genotype_tools

import snp
import mysql

# Complementarity map.
base_map = {
  "A": "T",
  "T": "A",
  "C": "G",
  "G": "C"
}

def get_snp_info(el_snp):
  # Index of each frequency column in rows returned from the database, e.g.
  # f_cases_aa is at 12 and f_controls_aa is at 15.
  attribute_index_map = {
    "aa": [12, 15],
    "ab": [13, 16],
    "bb": [14, 17],
  }

  snp_info = {}

  # Query params are %s with this library, not ?.
  query = "SELECT * FROM diseases.el_snps WHERE db_snp = %s;"
  db.execute(query, (el_snp.rsid,))
  res = db.fetchone()

  snp_info["alleles"] = res[2].split("/")
  for k, vs in attribute_index_map.items():
    snp_info[k] = [res[v] for v in vs]

  return snp_info

# Return (f_cases_*, f_controls_*) for this SNP.
def get_probabilities(el_snp, flipped = False):
  # From their results.
  snp_info = get_snp_info(el_snp)
  print("Calculating for %s, having genotype %s, vs. %s." % (el_snp.rsid, 
      el_snp.genotype, "".join(snp_info["alleles"])))

  if el_snp.genotype            == 2 * snp_info["alleles"][0]:
    return snp_info["aa"]
  elif sorted(el_snp.genotype)  == sorted(snp_info["alleles"]):
    return snp_info["ab"]
  elif el_snp.genotype          == 2 * snp_info["alleles"][1]:
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

if len(sys.argv) != 2:
  print("Usage: python el.py path/to/user/genome/file.txt")
  sys.exit(2)

user_snps_path  = sys.argv[1].strip()

if not os.path.exists(user_snps_path):
  print("No file found at %s." % (user_snps_path))
  sys.exit(2)

populations = ("CEU", "HPT", "YRI",)
population_map = dict((p[0], p) for p in populations)
population_raw = raw_input("Select a population (C)EU, (Y)RI, (H)PT: ")
population = population_map[population_raw.upper()]


user_snps = genotype_tools.FileUtils.read_genotype_file(user_snps_path)
db.execute("SELECT db_snp FROM diseases.el_snps")
el_rsids  = [row[0] for row in db.fetchall()]

# Now, we want to compute 
#   a = prior(EL) * product(i = 1 to 150, probability(SNP-i | EL))
#   b = prior(AL) * product(i = 1 to 150, probability(SNP-i | AL))
#   a / (a + b)
el_probability = 1

# For each extreme longevity SNP, we get the value of that SNP from the input
# and adjust the running EL and AL scores according to those values from the
# paper.
#import pdb
#pdb.set_trace()
for el_rsid in el_rsids:
  # If we already have a value for this SNP - use it.
  el_snp = user_snps.get(el_rsid, None)
  
  # Impute if we have to...
  if el_snp is None:
    try:
      el_snp = genotype_tools.impute_rsid_simple(user_snps, el_rsid, population)
    # TODO: workaround error with rs2042831 and CEU from Mikolaj Habryn.
    except ValueError, e:
      print("Error occurred imputing for %s: %s" % (el_rsid, e))
      continue

  # Imputation returned None (this should never happen).
  if el_snp is None:
    print "Unable to impute value for el_rsid: imputation returned None."
    continue
  
  # Update probabilities - multiply each by the probability of having this
  # genotype given each condition (EL or AL).
  # Avoiding underflow at each step...
  probabilities = get_probabilities(el_snp)
  print(probabilities)
  el_probability *= (probabilities[0] / probabilities[1])
  #el_probability *= (probabilities[0] / (probabilities[0] + probabilities[1]))

print("Chance of living to 100: %f" % (el_probability / (1 - el_probability)))
