import cgi
import sys
import os
import genotype_tools
import subprocess
import urllib

import snp
import mysql

# Complementarity map.
BASE_MAP = {
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
  db.execute(query, (el_snp.rsid[2:],))
  res = db.fetchone()

  snp_info["alleles"] = res[2].split("/")
  for k, vs in attribute_index_map.items():
    snp_info[k] = [res[v] for v in vs]

  return snp_info

# Return (f_cases_*, f_controls_*) for this SNP.
def get_probabilities(el_snp, flipped=False):

  # Get information about this SNP from the paper's supplementals.
  snp_info = get_snp_info(el_snp)

  #print("Calculating for %s, having genotype %s, vs. %s." % (el_snp.rsid, 
  #    el_snp.genotype, "".join(snp_info["alleles"])))

  # We have two of the "AA" alleles.
  if el_snp.genotype            == 2 * snp_info["alleles"][0]:
    return snp_info["aa"]
  # We have one of each type of allele.
  elif sorted(el_snp.genotype)  == sorted(snp_info["alleles"]):
    return snp_info["ab"]
  # We have two of the "BB" alleles.
  elif el_snp.genotype          == 2 * snp_info["alleles"][1]:
    return snp_info["bb"]

  else:
    # If we haven't flipped the bases of the user's EL SNP, and we haven't 
    # found a match yet, flip now and call this function again.
    if not flipped:
      el_snp.genotype = ''.join([BASE_MAP[a] for a in el_snp.genotype])
      return get_probabilities(el_snp, True)

    print("Error occurred matching user SNP to EL SNP %s!" % el_snp.rsid)
    return (1, 1)

GOOGLE_CHARTS_BASE_URL = "http://chart.apis.google.com/chart"

GOOGLE_CHART_PARAMS = {
  "cht": "lc",
  "chs":  "630x240",
  "chtt": "Probability of Extreme Longevity",
  "chma": "15,5,5,5",
  "chxt": "x,x,y,y",
  "chxr": "0,1,150|2,0,1",
  "chxl": "1:|Number of SNPs|3:|Proability of EL|",
  "chxp": "1,50|3,50",
}

def build_chart_url(el_running_odds):
  odds_series = ",".join(str(round(i, 4)) for i in el_running_odds)

  GOOGLE_CHART_PARAMS["chd"] = "t:%s" % (odds_series,)
  query_string = "&".join(
    e[0] + "=" + e[1] for e in GOOGLE_CHART_PARAMS.items()
  )

  return GOOGLE_CHARTS_BASE_URL + "?" + query_string

  

database = mysql.connector.Connect(host='marlowe', user='gene210-user', 
                                   passwd='genomics', buffered=True)
db = database.cursor()

population = None
populations = ("CEU", "YRI")
population_choices = "(" + ", ".join(populations) + ")"

if len(sys.argv) < 2:
  print("Usage: python el.py path/to/user/genome/file.txt " + 
        "[optional population identifier (you will be prompted if not given)]")
  sys.exit(2)

if len(sys.argv) == 3:
  population = sys.argv[2].strip().upper()

if population is None:
  population = raw_input("Select a population %s: " % population_choices)
  population = population.strip().upper()

if population not in populations:
  print("Population identifier %s is not valid." % population)
  sys.exit(1)

user_genome_path  = sys.argv[1].strip()
if not os.path.exists(user_genome_path):
  print("No file found at %s." % (user_genome_path))
  sys.exit(2)

print("Running (this may take some time)... ")
user_snps = genotype_tools.FileUtils.read_genotype_file(user_genome_path)
db.execute("SELECT db_snp FROM diseases.el_snps")
el_rsids  = [row[0] for row in db.fetchall()]

# Now, we want to compute 
#   a = prior(EL) * product(i = 1 to 150, Pr(SNP-i | EL))
#   b = prior(AL) * product(i = 1 to 150, Pr(SNP-i | AL))
#   a / (a + b)
# To avoid underflow, we don't do precisely this - instead, we start with an
# odds ratio of 1 (= 0.5 / 0.5, from the paper's priors on Pr(EL)
# and Pr(AL)), and then adjust that each time by 
# Pr(SNP-i | EL) / Pr(SNP-I | AL).
el_odds = 1

# Store the intermediate results for plotting, starting with 50%
el_running_odds = [50]

# For each EL-SNP, we get the value of that SNP from the provided genome 
# and adjust the running EL and AL scores according to those values from the
# paper.
for el_rsid in el_rsids:
  # If we already have a value for this SNP - use it.
  user_el_snp = user_snps.get("rs" + el_rsid, None)
  
  # Impute if we have to...
  if user_el_snp is None:
    try:
      user_el_snp = genotype_tools.impute_rsid_simple(user_snps, "rs" + el_rsid,
                                                      population)
    # TODO: workaround error with rs2042831 and CEU from Mikolaj Habryn.
    except ValueError, e:
      print("Error occurred imputing for %s: %s." % (el_rsid, e))
      continue

  # Imputation returned None (this should never happen).
  if user_el_snp is None:
    print "Unable to impute value for %s: imputation returned None." % el_rsid
    continue
  
  # Update probabilities - multiply by the ratio of the probability of EL
  # given this genotype / the probability of AL given this genotype.
  probabilities = get_probabilities(user_el_snp)
  
  el_odds *= (probabilities[0] / probabilities[1])
  # Store the running total for display later
  el_running_odds.append((100 * (el_odds / (1 + el_odds))))

print("Percentage chance of living to 100: %f" % \
      (100 * (el_odds / (1 + el_odds))))

out_file_name = "extreme_longevity_test.html"
out_file = open(out_file_name, "w")

out_file.write("<html><body><h3 style='text-align:center;'>Estimated Percentage Chance of Living to be > 100</h3>")
out_file.write("<img src='" + chart_url + "'</img>")
out_file.write("</body></html>")

out_file.close()
subprocess.Popen(("open", out_file_name)).wait()
