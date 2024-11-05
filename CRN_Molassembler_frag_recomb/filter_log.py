import json

with open("log.json","r") as f:
    recs = json.load(f)

no3_ion_bi = {'smiles':'[N+](=O)([O-])[O-]', 
          'coordList':[1,2],
          'ligType':'bi_cis'}
no3_ion_mono = {'smiles':'[N+](=O)([O-])[O-]', 
          'coordList':[1]}
hno3_bi = {
    'smiles': '[N+](=O)(O)[O-]',
    'coordList':[1,3],
    'ligType':'bi_cis'
}
hno3_mono = {
    'smiles': '[N+](=O)(O)[O-]',
    'coordList':[1]
}
moe_ion_bi = {'smiles':'COCC[O-]', 
          'coordList':[1,4], 
          'ligType':'bi_cis'}
moe_ion_mono = {'smiles':'COCC[O-]', 
          'coordList':[4]}
moe_0_bi = {'smiles':'COCCO', 
          'coordList':[1,4], 
          'ligType':'bi_cis'}
moe_0_mono = {'smiles':'COCCO', 
          'coordList':[4]}
moe_0_c_mono = {'smiles':'COCCO', 
          'coordList':[1]}

ligand_library = {"no3_ion_bi": no3_ion_bi,
           "no3_ion_mono": no3_ion_mono,
           "hno3_bi": hno3_bi,
           "hno3_mono": hno3_mono,
           "moe_ion_bi": moe_ion_bi,
           "moe_ion_mono": moe_ion_mono,
           "moe_0_bi": moe_0_bi,
           "moe_0_mono": moe_0_mono,
           "moe_0_c_mono": moe_0_c_mono}

rec_compl = []
rec_idx = []
for rec in recs:
    if 'moe_ion_c_mono' in rec['ligands'].keys():
        continue
    ligs = []
    srt_ligands = dict(sorted(rec['ligands'].items()))
    for key in srt_ligands.keys():
        ligs = ligs + [ligand_library[key]]*srt_ligands[key]
    complex_dict = {'core':{'metal':'Bi',
                    'coreCN':rec['CN']},
                   'ligands': ligs
                   }
    if complex_dict not in rec_compl:
        rec_compl.append(complex_dict)
        rec_idx.append(rec['inter_idx'])
    #print(complex_dict)
print(len(rec_compl)==len(rec_idx))
print(len(rec_compl))

filtered_rec = [recs[i] for i in rec_idx]

with open("filtered_log.json","w") as f:
    json.dump(filtered_rec, f)