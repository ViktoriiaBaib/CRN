import pickle
import time
import math
from copy import deepcopy
from scine_molassembler import *
import scine_utilities as sc_utils
from collections import defaultdict, Counter
import json

"""# Functions

#### Fragmentation functions
"""

def remove_ligand(scine_molecule_input, bi_position, bonded_position):
  """
  scine_molecule_input: scine_molassembler.Molecule
  bi_position: int
  bonded_position: int
  This function breaks input molecule in two along the provided Bi-O bond
  output is tuple with 0 element being bigger molecule and 1 element being ligand
  """
  scine_molecules = editing.cleave(scine_molecule_input, BondIndex(bi_position, bonded_position))
  fragment = scine_molecules[0]
  if 'Bi' not in {str(a) for a in fragment.graph.elements()}:
    fragment = scine_molecules[1]
  _ = fragment.canonicalize()
  return fragment

def break_bond(scine_molecule_input, bi_position, bonded_position):
  """
  scine_molecule_input: scine_molassembler.Molecule
  bi_position: int
  bonded_position: int
  Breaks bond, returning connected graph or bigger fragment
  """
  if not scine_molecule_input.graph.can_remove(BondIndex(bi_position, bonded_position)):
    #print("remove ligand, break bond between: ", bi_position, bonded_position)
    return remove_ligand(scine_molecule_input, bi_position, bonded_position)
  else:
    #print("remove bidentality")
    fragment = deepcopy(scine_molecule_input)
    fragment.remove_bond(BondIndex(bi_position, bonded_position))
    _ = fragment.canonicalize()
    return fragment

def generate_fragments_for_single_recombinant(scine_molecule_input):
  """
  scine_molecule_input: scine_molassembler.Molecule
  all_fragments: set of scine_molassembler.Molecule objects
  This function checks if Bi connected to only O, then locates all such bonds and generates fragments
  """
  #locate Bi
  bi_poss = scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.Bi)
  if len(bi_poss) == 1:
    bi_pos = bi_poss[0]
  else:
    raise ValueError('Number of Bi atoms is different then 1. Check your graph!')
  bonded_positions = [a for a in scine_molecule_input.graph.adjacents(bi_pos)]
  #check that Bi is bonded to O only:
  bonded_atoms = set([scine_molecule_input.graph.element_type(a) for a in scine_molecule_input.graph.adjacents(bi_pos)])
  if not bonded_atoms == {sc_utils.ElementType.O}:
    raise ValueError('Bi is bonded to something besides O. Check your graph!')
  #print("Bi position: ", bi_pos)
  #print("bonded posistions: ", bonded_positions)
  new_fragments = set([break_bond(scine_molecule_input, bi_pos, bonded_position) for bonded_position in bonded_positions])
  #print("One recombinant is OK!")
  return new_fragments

def generate_fragments(new_recombinants):
  new_fragments = set([])
  for recombinant in new_recombinants:
    new_fragments = new_fragments | generate_fragments_for_single_recombinant(recombinant)
  return new_fragments

"""#### Recombination functions"""

def add_inter_molecular_bonds(scine_molecule_input, bi_position, connectible_position):
  scine_molecule = deepcopy(scine_molecule_input)
  scine_molecule.add_bond(bi_position, connectible_position)
  _ = scine_molecule.canonicalize()
  return scine_molecule

def add_ligand(scine_molecule_input, bi_pos, ligand):
  """
  scine_molecule_input: scine_molassembler.Molecule
  bi_pos: int
  ligand: (scine_molassembler.Molecule, [int, int, ...])
  This function adds ligand (MOE[-] O[H]CH2CH2OCH3) to the Bi in scine molecule by connecting to all pos.
  """
  recombinants = []
  for pos in ligand[1]:
    scine_molecule = editing.connect(scine_molecule_input, ligand[0], bi_pos, pos, BondType.Single)
    _ = scine_molecule.canonicalize()
    recombinants.append(scine_molecule)
  return set(recombinants)

def generate_recombinants_for_single_fragment(scine_molecule_input, ligands):
  """
  scine_molecule_input: scine_molassembler.Molecule
  ligands: [(scine_molassembler.Molecule, [int, int, ...]), (scine_molassembler.Molecule, [int, int, ...])]
  all_recombinants: set of scine_molassembler.Molecule objects
  """
  #locate connectible atoms
  #locate Bi
  bi_poss = scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.Bi)
  if len(bi_poss) == 1:
    bi_pos = bi_poss[0]
  else:
    raise ValueError('Number of Bi atoms is different then 1. Check your graph!')
  #intermolecular bonds
  bonded_positions = {a for a in scine_molecule_input.graph.adjacents(bi_pos)}
  o_positions = set(scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.O))
  connectible_positions = o_positions - bonded_positions
  new_recombinants = set([add_inter_molecular_bonds(scine_molecule_input, bi_pos, con_pos) for con_pos in connectible_positions])
  #ligands
  for ligand in ligands:
    new_recombinants = new_recombinants | add_ligand(scine_molecule_input, bi_pos, ligand)
  return new_recombinants

def connectible_atoms_ligand(scine_molecule_input):
  o_pos_all = scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.O)
  o_pos = o_pos_all
  if len(o_pos_all)>2:
    o_pos_mask = ['H' not in [str(scine_molecule_input.graph.element_type(a)) for a in scine_molecule_input.graph.adjacents(ij)] for ij in o_pos_all]
    o_pos = [o_pos_all[i] for i in range(len(o_pos_all)) if o_pos_mask[i]]
  return (scine_molecule_input, o_pos)

def connectible_atoms_ligand_ion_condition(scine_molecule_input):
  o_pos_all = scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.O)
  o_pos = o_pos_all
  #o_pos_mask = [len([a for a in scine_molecule_input.graph.adjacents(ij)])<2 for ij in o_pos_all]
  #o_pos = [o_pos_all[i] for i in range(len(o_pos_all)) if o_pos_mask[i]]
  if len(o_pos)>2:
    o_pos = o_pos[:2]
  return (scine_molecule_input, o_pos)

def generate_recombinants(filtered_new_fragments, ligands):
  """
  All inputs are sets of scine_molassembler.Molecule objects
  """
  new_recombinants = set([])
  for fragment in filtered_new_fragments:
    new_recombinants = new_recombinants | generate_recombinants_for_single_fragment(fragment, ligands)
  return new_recombinants

def recombinant_cn_ok(scine_molecule_input):
  bi_poss = scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.Bi)
  if len(bi_poss) == 1:
    bi_pos = bi_poss[0]
  else:
    raise ValueError('Number of Bi atoms is different then 1. Check your graph!')
  bonded_positions = [a for a in scine_molecule_input.graph.adjacents(bi_pos)]
  return (len(bonded_positions) <= 8) and (len(bonded_positions) >= 5)

def recombinant_is_neutral(scine_molecule_input, el_charges):
  """
  scine_molecule_input: scine_molassembler.Molecule
  el_charges: dictionary el_charges = {'Bi': 3.0, 'O': -2.0, 'N': 5.0, 'H': 1.0, 'C': -1.333333}
  """
  mol_elements = [str(a) for a in scine_molecule_input.graph.elements()]
  charge = sum([el_charges[me] for me in mol_elements])
  return math.isclose(charge, 0.0, abs_tol=1e-4)

def no_tri_dentate_bonded_atoms(scine_molecule_input):
  """
  scine_molecule_input: scine_molassembler.Molecule
  This function checks if there are atoms which have more than 2 common B-O bonds.
  """
  bi_poss = scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.Bi)
  if len(bi_poss) != 1:
    raise ValueError('Number of Bi atoms is different then 1. Check your graph!')
  bi_pos = bi_poss[0]
  bonded_positions = {a for a in scine_molecule_input.graph.adjacents(bi_pos)}
  second_bonded_positions = []
  _=[second_bonded_positions.extend(list(scine_molecule_input.graph.adjacents(b))) for b in bonded_positions]
  second_bonded_positions_nobi = [pos for pos in second_bonded_positions if pos != bi_pos]
  counts = sum([second_bonded_positions_nobi.count(pos)>2 for pos in set(second_bonded_positions_nobi)])
  return not bool(counts)

# def recombinant_ok_for_3d(scine_molecule_input, n_structures, seed):
#   """
#   scine_molecule_input: scine_molassembler.Molecule
#   n_structures: int
#   seed: int
#   """
#   results = dg.generate_ensemble(scine_molecule_input, n_structures, seed)
#   n_fails = sum([1 if isinstance(r, dg.Error) else 0 for r in results])
#   return n_fails < n_structures

def filter_out_recombinants(new_recombinants, el_charges, n_structures = 20, seed = 1010):
  """
  new_recombinants: set of scine_molassembler.Molecule objects
  This function filters out recombinants with CN>8, charge != 0 and failing to get 3d
  """
  #check Bi CN
  CN_recombinants = [recombinant for recombinant in new_recombinants if recombinant_cn_ok(recombinant)]
  #print("After CN checking, number of recombinants is: ", len(CN_recombinants))
  #check charge
  filtered_recombinants = [recombinant for recombinant in CN_recombinants if recombinant_is_neutral(recombinant, el_charges)]
  #print("After charge checking, number of recombinants is: ", len(filtered_recombinants))
  #check that 3d is generatable
  filtered_new_recombinants = [recombinant for recombinant in filtered_recombinants if no_tri_dentate_bonded_atoms(recombinant)]
  #print("After 3d structures checking, number of recombinants is: ", len(filtered_new_recombinants))
  #print("Produced number of recombinants is: ", len(filtered_new_recombinants))
  return set(filtered_new_recombinants)


# Functions for ligand parsing

def def_ligand(scine_molecule_input, bi_position, bonded_position, lig_atoms):
  """
  scine_molecule_input: scine_molassembler.Molecule
  bi_position: int
  bonded_position: int
  This function breaks input molecule in two along the provided Bi-O bond
  then by the size of the ligand it determines what it is
  """
  suffix = ''
  scine_molecules = editing.cleave(scine_molecule_input, BondIndex(bi_position, bonded_position))
  adj_at = ''.join([str(scine_molecule_input.graph.element_type(a)) for a in scine_molecule_input.graph.adjacents(bonded_position) if a != bi_position])
  if adj_at == 'CC':
    suffix = '_c'
  fragment = scine_molecules[1]
  if 'Bi' in {str(a) for a in fragment.graph.elements()}:
    fragment = scine_molecules[0]
  return lig_atoms[fragment.graph.V]+suffix

def build_adjacency_list(edges):
  adj_list = defaultdict(list)
  for edge in edges:
    adj_list[edge[0]].append(edge[1])
    adj_list[edge[1]].append(edge[0])  # because it's an undirected graph
  return adj_list

def find_all_cycles(start, current, visited, stack, adj_list, cycles, bi_pos):
  visited[current] = True
  stack.append(current)
  for neighbor in adj_list[current]:
    if visited[neighbor] == False:
      find_all_cycles(start, neighbor, visited, stack, adj_list, cycles, bi_pos)
    elif neighbor == start and len(stack) >= 3:  # check for length >= 3 because start will be added, making it a 4-cycle
      new_cycle = stack + [neighbor]
      new_cycle.remove(bi_pos)
      new_cycle = list(set(new_cycle))
      new_cycle.sort()
      if new_cycle not in cycles:
        cycles.append(new_cycle)
  stack.pop()
  visited[current] = False

def get_all_cycles(edges, bi_pos):
  adj_list = build_adjacency_list(edges)
  visited = {node: False for node in adj_list}
  cycles = []
  for vertex in adj_list:
    find_all_cycles(vertex, vertex, visited, [], adj_list, cycles, bi_pos)
    visited[vertex] = True
  return cycles

def parse_ligand_list(scine_molecule_input, lig_atoms):
  ligands = []
  bi_position = scine_molecule_input.graph.atoms_of_element(sc_utils.ElementType.Bi)[0]
  bi_CN = scine_molecule_input.graph.degree(bi_position)
  #print(bi_position, bi_CN)
  shell1 = [a for a in scine_molecule_input.graph.adjacents(bi_position)]
  #print(shell1)
  mono1 = [not scine_molecule_input.graph.can_remove(BondIndex(bi_position, bonded_position)) for bonded_position in shell1]
  #print(mono1)
  for i in range(bi_CN):
    if mono1[i]:
      bonded_position = shell1[i]
      ligands.append(def_ligand(scine_molecule_input, bi_position, bonded_position,lig_atoms)+"_mono")
  #print(ligands)
  edges = [b for b in scine_molecule_input.graph.bonds()]
  adj_list = build_adjacency_list(edges)
  cycles = get_all_cycles(edges, bi_position)
  for cycle in cycles:
    lig_at = []
    for at in cycle:
      if at != bi_position:
        lig_at = lig_at + [at2 for at1 in adj_list[at] if at1 != bi_position for at2 in adj_list[at1]]
    lig_at = list(set(lig_at))
    lig_at.remove(bi_position)
    ligands.append(lig_atoms[len(lig_at)]+"_bi")
  return {'CN':bi_CN, 'ligands':dict(Counter(ligands))}

"""# Main

#### Prepare structures

"""
Options.chiral_state_preservation = ChiralStatePreservation.Unique
#calculations/molassembler/bi3nitrcoord/bi_2w_3nitr.xyz
start_time = time.time()
# make 2d graph with molassembler
bino3 = io.read("bino3.xyz")
_ = bino3.canonicalize()
#bino3

moe_ion = io.read("moe_ion.xyz")
_ = moe_ion.canonicalize()
#moe_ion

moe = io.read("moe.xyz")
_ = moe.canonicalize()

hno3 = io.read("hno3.xyz")
_ = hno3.canonicalize()

no3_ion = io.read("no3.xyz")
_ = no3_ion.canonicalize()

#print("Preparation of sructures took time --- %s seconds ---" % (time.time() - start_time))

"""## Loop"""

start_time = time.time()
directory = "/global/home/groups/lr_mp/vbaibakova/molassember/reaction2/"
outfile = directory + 'stat.txt'
outfile_lig = directory + 'log.json'
outfile_rec = directory + 'recombinats.json'
el_charges = {'Bi': 3.0, 'O': -2.0, 'N': 5.0, 'H': 1.0, 'C': -1.333333}
#with open(directory + "7_all_recombinants.pkl", "rb") as fr:
all_recombinants = set([bino3])
filtered_new_recombinants = set([bino3])
all_fragments = set([])
ligands = [connectible_atoms_ligand_ion_condition(ligand_ion) for ligand_ion in set([moe_ion,no3_ion])] + [connectible_atoms_ligand(ligand) for ligand in set([moe, hno3])]
iteration = 0
lfnr = len(filtered_new_recombinants)
with open(outfile, 'w') as fout:
  fout.write("Iteration | New frag | Filt new fr | All frag | New recomb | Filt new rec | All recomb | Exec Time | Fragm Time | Filt frag Time | Recomb Time | Filt recomb Time\n------------------------------------------------------------------------------------------------\n")
while lfnr > 0:
  start_iter_time = time.time()
  new_fragments = generate_fragments(filtered_new_recombinants)
  lnf = len(new_fragments)
  frag_time = time.time() - start_iter_time
  start_oper_time = time.time()
  filtered_new_fragments = new_fragments - all_fragments
  del new_fragments
  lfnf = len(filtered_new_fragments)
  filt_frag_time = time.time() - start_oper_time
  all_fragments = all_fragments | filtered_new_fragments
  start_oper_time = time.time()
  new_recombinants = filtered_new_fragments | generate_recombinants(filtered_new_fragments, ligands)
  lnr = len(new_recombinants)
  del filtered_new_fragments
  gen_recomb_time = time.time() - start_oper_time
  start_oper_time = time.time()
  filt_new_recombinants = new_recombinants - all_recombinants
  del new_recombinants
  filtered_new_recombinants = filter_out_recombinants(filt_new_recombinants, el_charges)
  del filt_new_recombinants
  lfnr = len(filtered_new_recombinants)
  filt_recomb_time = time.time() - start_oper_time
  all_recombinants = all_recombinants | filtered_new_recombinants
  #del filtered_new_recombinants
  exec_time = time.time() - start_iter_time
  with open(outfile, 'a') as fout:
    fout.write("    %s    |    %s    |    %s    |    %s    |    %s    |    %s    |    %s    |    %s    |    %s    |    %s    |    %s    |    %s    \n" % (iteration, lnf, lfnf, len(all_fragments), lnr, lfnr, len(all_recombinants), int(exec_time), int(frag_time), int(filt_frag_time), int(gen_recomb_time), int(filt_recomb_time)))
  iteration = iteration + 1
with open(outfile, 'a') as fout:
  fout.write("\n\n -- END --\nTotal execution time --- %s seconds ---" % (time.time() - start_time))

intermediates = [rec for rec in all_recombinants]
records = []
#save
for idx,intermediate in enumerate(intermediates):
  records.append({"inter_idx":idx, "graph": intermediate.dump_graphviz()})
with open(outfile_rec, "w") as outf:
  json.dump(records, outf)

## Parse ligands

lig_atoms = {12: 'moe_ion', 4: 'no3_ion', 13: 'moe_0', 5: 'hno3', 3: 'water'}

parsed_ligs = []

print("Parsing ligans")
for idx,intermediate in enumerate(intermediates):
  print(idx)
  pll = parse_ligand_list(intermediate, lig_atoms)
  pll['inter_idx'] = idx
  print(pll)
  parsed_ligs.append(pll)
with open(outfile_lig, "w") as outf:
  json.dump(parsed_ligs, outf)

#with open(directory+"all_recombinants.pkl", "wb") as fr:
#  pickle.dump(all_recombinants, fr)
