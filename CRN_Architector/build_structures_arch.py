import argparse, json, time, math
import architector.io_obabel as io_obabel
from architector import (build_complex,smiles2Atoms,convert_io_molecule)
from sella.optimize import Sella
import itertools
from random import sample
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--batchnum', type=int)
parser.add_argument('-b', '--batchsize', type=int)
args = parser.parse_args()

batchnum = args.batchnum
batchsize = args.batchsize

opt = Sella

no3_ion_bi = {'smiles':'[N+](=O)([O-])[O-]','coordList':[1,2],'ligType':'bi_cis'}
no3_ion_mono = {'smiles':'[N+](=O)([O-])[O-]','coordList':[1]}
hno3_bi = {'smiles': '[N+](=O)(O)[O-]','coordList':[1,3],'ligType':'bi_cis'}
hno3_mono = {'smiles': '[N+](=O)(O)[O-]','coordList':[1]}
moe_ion_bi = {'smiles':'COCC[O-]','coordList':[1,4],'ligType':'bi_cis'}
moe_ion_mono = {'smiles':'COCC[O-]','coordList':[4]}
moe_0_bi = {'smiles':'COCCO','coordList':[1,4],'ligType':'bi_cis'}
moe_0_mono = {'smiles':'COCCO','coordList':[4]}
moe_0_c_mono = {'smiles':'COCCO','coordList':[1]}

ligand_library = {
"no3_ion_bi": no3_ion_bi,
"no3_ion_mono": no3_ion_mono,
"hno3_bi": hno3_bi,
"hno3_mono": hno3_mono,
"moe_ion_bi": moe_ion_bi,
"moe_ion_mono": moe_ion_mono,
"moe_0_bi": moe_0_bi,
"moe_0_mono": moe_0_mono,
"moe_0_c_mono": moe_0_c_mono
}

directory = "/global/scratch/users/vbaibakova/Architector/"
with open(directory+"filtered_log.json","r") as f:
    recs = json.load(f)

outfile = directory + f"stat{batchnum}.txt"

os.environ["MKL_NUM_THREADS"]="1"
os.environ["NUMEXPR_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"

start = batchnum*batchsize
end = (batchnum+1)*batchsize
if end>len(recs):
    end = len(recs)
for rec in recs[batchnum*batchsize:(batchnum+1)*batchsize]:
    ligs = []
    for key in rec['ligands'].keys():
        ligs = ligs + [ligand_library[key]]*rec['ligands'][key]
    complex_dict = {
    'core':{'metal':'Bi','coreCN':rec['CN']},
    'ligands': ligs,
    'parameters':{
    'ff_preopt':True,
    'xtb_solvent':'octanol',
    'n_symmetries': 20,
    'n_conformers': 10,
    'max_steps': 25,
    'ase_opt_method':opt,
    'ase_opt_kwargs':{'sella_internal_trics':True, 'order':0}}}
    try:
        t0 = time.time()
        out = build_complex(complex_dict)
        t1t0 = (time.time()-t0)//3600
    except:
        print(f"Failed for graph {rec['inter_idx']}!")
        continue
    with open(outfile, 'w') as fout:
        fout.write("Graph | Hours | Structures \n -------------------------------\n")
        fout.write(f" {rec['inter_idx']}  |  {t1t0}  |  {len(out)} ")
    for key in out.keys():
        print(key)
        out[key]['ase_atoms'].write(directory+f"str/graph{rec['inter_idx']}_cn{rec['CN']}_{key}.xyz")
    del out
