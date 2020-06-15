# Install:
# python3 -m pip install gensim goatools wget
#
# Run:
# python3 go-analysis.py

import os
import wget
from gensim.models import Word2Vec
from goatools import obo_parser
from goatools.semantic import deepest_common_ancestor
from opencog.atomspace import AtomSpace, types
from opencog.type_constructors import *
from opencog.scheme_wrapper import scheme_eval
from scipy.spatial import distance

go_obo_url = "http://geneontology.org/ontology/go-basic.obo"
go_obo_path = os.getcwd() + "/kbs/go-basic.obo"
w2v_model_path = os.getcwd() + "/model51.bin"
pln_result_path = os.getcwd() + "/results/pln-intensional-differences.csv"
go_cardinality_path = os.getcwd() + "/results/go-cardinality.csv"
go_pair_path = os.getcwd() + "/go-pairs.txt"
output_csv_path = os.getcwd() + "/results/go.csv"
kbs = [
  "kbs/GO_2020-04-01.scm",
  "kbs/GO_annotation_gene-level_2020-04-01.scm",
  "kbs/Go-Plus-GO_2020-05-04.scm"
]

if not os.path.isfile(go_obo_path):
  print("Downloading go-basic.obo...")
  wget.download(go_obo_url, go_obo_path)

go_dag = obo_parser.GODag(go_obo_path)
w2v_model = Word2Vec.load(w2v_model_path)

pln_results = {}
with open(pln_result_path) as f:
  for line in f:
    content = line.split(",")
    go1 = content[0].strip()
    go2 = content[1].strip()
    tv_strength = content[2].strip()
    tv_confidence = content[3].strip()
    pln_results[go1 + go2] = tv_strength + "," + tv_confidence

go_cardinality = {}
with open(go_cardinality_path) as f:
  for line in f:
    content = line.split(",")
    go = content[0].strip()
    cardinality = content[1].strip()
    go_cardinality[go] = cardinality

def dag_distance(go1, go2):
  common_go = deepest_common_ancestor([go1, go2], go_dag)
  go1_depth = go_dag[go1].depth
  go2_depth = go_dag[go2].depth
  common_go_depth = go_dag[common_go].depth
  return go1_depth + go2_depth - 2 * common_go_depth

def vec_distance(v1, v2):
  vec1 = w2v_model[v1]
  vec2 = w2v_model[v2]
  # vec1 = w2v_model.wv.word_vec(v1, use_norm=True)
  # vec2 = w2v_model.wv.word_vec(v2, use_norm=True)
  return distance.euclidean(vec1, vec2)

'''
atomspace = AtomSpace()

scheme_eval(atomspace, "(add-to-load-path \"/usr/share/guile/site/2.2/opencog\")")
scheme_eval(atomspace, "(add-to-load-path \".\")")

scheme_eval(atomspace, "(use-modules (opencog) (opencog bioscience) (opencog pln) (opencog ure))")
scheme_eval(atomspace, "(load \"bio-as-utils.scm\")")
for kb in kbs:
  scheme_eval(atomspace, "(primitive-load \"{}\")".format(kb))
'''

output_csv_fp = open(output_csv_path, "w")
first_row = ",".join([
  "GO pair",
  "GO_1 cardinality",
  "GO_2 cardinality",
  "GO_1 depth",
  "GO_2 depth",
  "Distance in DAG",
  "IntDiff strength",
  "IntDiff confidence",
  "Vector distance"
])
output_csv_fp.write(first_row + "\n")
with open(go_pair_path) as f:
  for line in f:
    content = line.split(",")
    go1 = content[0].strip()
    go2 = content[1].strip()
    go_pair = go1 + " " + go2
    go1_cardinality = go_cardinality[go1]
    go2_cardinality = go_cardinality[go2]
    if go_dag.get(go1) and go_dag.get(go2):
      go1_depth = go_dag[go1].depth
      go2_depth = go_dag[go2].depth
      tv = pln_results[go1 + go2].split(",")
      tv_strength = tv[0]
      tv_confidence = tv[1]
      dag_dist = dag_distance(go1, go2)
      vec_dist = vec_distance(go1, go2)
      row = ",".join([
        go_pair,
        go1_cardinality,
        go2_cardinality,
        str(go1_depth),
        str(go2_depth),
        str(dag_dist),
        str(tv_strength),
        str(tv_confidence),
        str(vec_dist)])
      # print(row)
      output_csv_fp.write(row + "\n")
