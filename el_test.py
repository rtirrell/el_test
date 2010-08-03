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

# Return (f_cases_*, f_controls_*) for this SNP.
def get_snp_info(el_snp, flipped=False):
  # Index of each frequency column in rows returned from the database, e.g.
  # f_cases_aa is at 12 and f_controls_aa is at 15.
  attribute_index_map = {
    "aa": [12, 15],
    "ab": [13, 16],
    "bb": [14, 17],
    "bayes_factor": [9]
  }

  snp_info = {}

  # Query params are %s with this library, not ?.
  query = "SELECT * FROM diseases.el_snps WHERE db_snp = %s;"
  db.execute(query, (el_snp.rsid[2:],))
  res = db.fetchone()

  snp_info["el_alleles"] = res[2].split("/")
  for k, vs in attribute_index_map.items():
    snp_info[k] = [res[v] for v in vs]
    if len(snp_info[k]) == 1:
      snp_info[k] = snp_info[k][0]

  # Get information about this SNP from the paper"s supplementals.
  snp_info["user_alleles"] = el_snp.genotype

  # We have two of the "AA" alleles.
  if el_snp.genotype            == 2 * snp_info["el_alleles"][0]:
    snp_info["probabilities"] = snp_info["aa"]
  # We have one of each type of allele.
  elif sorted(el_snp.genotype)  == sorted(snp_info["el_alleles"]):
    snp_info["probabilities"] = snp_info["ab"]
  # We have two of the "BB" alleles.
  elif el_snp.genotype          == 2 * snp_info["el_alleles"][1]:
    snp_info["probabilities"] = snp_info["bb"]

  else:
    # If we haven"t flipped the bases of the user"s EL SNP, and we haven"t 
    # found a match yet, flip now and call this function again.
    if not flipped:
      el_snp.genotype = "".join([BASE_MAP[a] for a in el_snp.genotype])
      snp_info = get_snp_info(el_snp, True)

    else:
      print "Error occurred matching user SNP to EL SNP %s!" % (el_snp.rsid,)
      snp_info["probabilities"] = (1, 1)

  el_top_frequency = 2 * snp_info["aa"][0] + snp_info["ab"][0]
  al_top_frequency = 2 * snp_info["aa"][1] + snp_info["ab"][1]

  if el_top_frequency >= al_top_frequency:
    snp_info["centenarian_allele"] = snp_info["el_alleles"][0]
  else:
    snp_info["centenarian_allele"] = snp_info["el_alleles"][1]

  return snp_info

GOOGLE_CHARTS_BASE_URL = "http://chart.apis.google.com/chart"

GOOGLE_CHART_PARAMS = {
  "cht": "lc",
  "chs":  "630x240",
  "chtt": "Probability of Exceptional Longevity",
  "chma": "15,5,5,5",
  "chxt": "x,x,y,y",
  "chxr": "0,1,150|2,0,1",
  "chxl": "1:|Number of SNPs|3:|Probability of EL|",
  "chxp": "1,50|3,50",
}

def build_chart_url(el_running_odds):
  odds_series = ",".join(str(round(i, 4)) for i in el_running_odds)

  GOOGLE_CHART_PARAMS["chd"] = "t:%s" % (odds_series,)
  query_string = "&".join(
    e[0] + "=" + e[1] for e in GOOGLE_CHART_PARAMS.items()
  )

  return GOOGLE_CHARTS_BASE_URL + "?" + query_string

database = mysql.connector.Connect(host="marlowe", user="gene210-user", 
                                   passwd="genomics", buffered=True)
db = database.cursor()

population = None
populations = ("CEU", "YRI")
population_choices = "(" + ", ".join(populations) + ")"

codename = None

if len(sys.argv) < 2:
  print "Usage: python el_test.py path/to/user/genome/file.txt",
  print "[population identifier] [codename]"
  print "If the population identifier or codename are not given, you will be",
  print "prompted for their values."
  sys.exit(2)

if len(sys.argv) == 4:
  population = sys.argv[2].strip().upper()
  codename   = sys.argv[3].strip()

if population is None:
  population = raw_input("Select a population %s: " % population_choices)
  population = population.strip().upper()

if population not in populations:
  print "Population identifier %s is not valid." % (population,)
  sys.exit(1)

if codename is None:
  codename = raw_input("Select a codename: ").strip()

user_genome_path  = sys.argv[1].strip()
if not os.path.exists(user_genome_path):
  print "No file found at %s." % (user_genome_path,)
  sys.exit(2)

print "Loading genotype file %s... " % (user_genome_path,),
sys.stdout.flush()
user_snps = genotype_tools.FileUtils.read_genotype_file(user_genome_path)
print("done!")

db.execute("SELECT db_snp FROM diseases.el_snps")
el_rsids  = [row[0] for row in db.fetchall()]

# Now, we want to compute 
#   a = prior(EL) * product(i = 1 to 150, Pr(SNP-i | EL))
#   b = prior(AL) * product(i = 1 to 150, Pr(SNP-i | AL))
#   a / (a + b)
# To avoid underflow, we don"t do precisely this - instead, we start with an
# odds ratio of 1 (= 0.5 / 0.5, from the paper"s priors on Pr(EL)
# and Pr(AL)), and then adjust that each time by 
# Pr(SNP-i | EL) / Pr(SNP-I | AL).
el_odds = 1

# Store the intermediate results for plotting, starting with 50%
el_running_odds = [50]

# For each EL-SNP, we get the value of that SNP from the provided genome 
# and adjust the running EL and AL scores according to those values from the
# paper.

# Track all SNPs for output.
results = {}

for el_rsid in el_rsids:
  # If we already have a value for this SNP - use it.
  user_el_snp = user_snps.get("rs" + el_rsid, None)
  
  # Impute if we have to...
  if user_el_snp is None:
    try:
      user_el_snp = genotype_tools.impute_rsid_simple(
        user_snps, "rs" + el_rsid, population
      )

    # TODO: workaround error with rs2042831 and CEU from Mikolaj Habryn.
    except ValueError, e:
      print "Error occurred imputing for rs%s: %s." % (el_rsid, e)
      continue

    # This SNP was imputed.
    print "Imputed %s." % (user_el_snp.nearest_SNP,)
    
  # Imputation returned None (this should never happen).
  if user_el_snp is None:
    print "Unable to impute value for rs%s: imputation returned None." % el_rsid
    continue
  
  # Update probabilities - multiply by the ratio of the probability of EL
  # given this genotype / the probability of AL given this genotype.
  snp_info = get_snp_info(user_el_snp)
  if user_el_snp.nearest_SNP is not None:
    snp_info["imputed_from"] = user_el_snp.nearest_SNP
  
  el_odds *= (snp_info["probabilities"][0] / snp_info["probabilities"][1])
  # Store the running total for display later
  el_running_odds.append((100 * (el_odds / (1 + el_odds))))

  results[user_el_snp.rsid] = snp_info

el_probability = 100 * (el_odds / (1 + el_odds))
print "Percentage chance of living to 100: %f" % (el_probability,)
      
db.execute("""
  INSERT INTO diseases.class_el_probabilities 
  (codename, population, probability) 
  VALUES ("%s", "%s", %f);
""" % (codename, population, el_probability))


# Name the output HTML file based on the input genome file"s name.
filename = os.path.splitext(os.path.split(user_genome_path)[1])[0]
out_file_name = "%s.html" % (filename,)
out_file = open(out_file_name, "w")

out_file.write("""
<html>
  <head>
    <title>
      Exceptional Longevity Exercise
    </title>
  </head>
  <body>
  <h3 style="text-align:left;">
    Estimated Probability of Exceptional Longevity
  </h3>
  Genome file: %s<br>
  Population given was: %s<br>
  Estimated probability of living > 100 years: %.4f%%<br><br>
""" % (filename, population, el_probability))

chart_url = build_chart_url(el_running_odds)

out_file.write("<img src='%s'></img><br><br>" % (chart_url,))
out_file.write("""
  <table>
    <tr>
      <th>rsID</th>
      <th>Centenarian Alleles</th>
      <th>Favorable Allele</th>
      <th>Bayes Factor</th>
      <th>User Alleles</th>
      <th>Imputed From</th>
    </tr>
""")

for (rsid, snp) in results.items():
  out_file.write("""
    <tr>
      <td>rs%s</td>
      <td>%s</td>
      <td>%s</td>
      <td>%1.2f</td>
      <td>%s</td>
  """ % (rsid, "".join(snp["el_alleles"]), snp["centenarian_allele"], 
         snp["bayes_factor"], snp["user_alleles"]))

  if "imputed_from" in snp:
    out_file.write("<td>%s</td>" % (snp["imputed_from"]))

  out_file.write("</td>")

out_file.write("</table></body></html>")

out_file.close()
subprocess.Popen(("open", out_file_name)).wait()
