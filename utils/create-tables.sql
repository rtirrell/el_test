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
);

CREATE TABLE class_el_probabilities ( 
  codename VARCHAR(255) NOT NULL, 
  population VARCHAR(255) NOT NULL,
 probability REAL NOT NULL 
); 

