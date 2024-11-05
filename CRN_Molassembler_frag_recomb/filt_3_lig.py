import json

with open("filtered_log.json","r") as f:
    recs = json.load(f)

frecs = []
for rec in recs:
    keys_moe = [key for key in rec['ligands'].keys() if "moe" in key]
    keys_no3 = [key for key in rec['ligands'].keys() if "no3" in key]
    n_moe = sum([rec['ligands'][key] for key in keys_moe])
    n_no3 = sum([rec['ligands'][key] for key in keys_no3])
    if (n_moe > 3) or (n_no3 > 3):
        print(rec)
        continue
    else:
        frecs.append(rec)
print(len(frecs))
with open("filtered_3_log.json", "w") as outf:
  json.dump(frecs, outf)