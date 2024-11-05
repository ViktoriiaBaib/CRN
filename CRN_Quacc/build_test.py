from ase.io import read
import jobflow as jf
from jobflow.managers.fireworks import flow_to_workflow
from fireworks import LaunchPad
from quacc.recipes.qchem.core import static_job, relax_job
from quacc.recipes.qchem.core import freq_job as qc_freq_job
from quacc.recipes.tblite.core import freq_job
import os

RUN = True

#directory = "/global/scratch/users/vbaibakova/Architector/str"
directory = "/global/home/groups/lr_mp/vbaibakova/bfo_newquacc/xyz"
classtag = "bfo_test"

print(directory)

for filename in os.listdir(directory):
    name=filename[:-4]
    f = os.path.join(directory, filename)
    atoms = read(f, format='xyz')
    charges = [0]*len(atoms.arrays['numbers'])
    atoms.set_initial_charges(charges)
    #freq TBLite
    f_job = freq_job(atoms, method="GFN2-xTB", pressure=0.986923, calc_swaps={"solvent":"octanol"})
    f_job.update_metadata({"class":classtag})
    f_job.name = name+"_TBLite-Freq_no-opt"
    #freq QChem
    f_qc = qc_freq_job(atoms, charge=0, spin_multiplicity=1, pcm_dielectric="16.93", overwrite_inputs={"rem": {"ecp":"def2-ecp", "resp_charges":"false"}})
    f_qc.update_metadata({"class":classtag})
    f_qc.name =  name+"_QChem-Freq_no-opt"
    #sp QChem
    s_job = static_job(atoms, charge=0,spin_multiplicity=1, pcm_dielectric="16.93", overwrite_inputs={"rem": {"ecp":"def2-ecp", "resp_charges":"false"}})
    s_job.update_metadata({"class":classtag})
    s_job.name = name+"_QChem-SP_no-opt"
    #opt QChem
    opt_job = relax_job(atoms, charge=0,spin_multiplicity=1, pcm_dielectric="16.93", overwrite_inputs={"rem": {"ecp":"def2-ecp", "resp_charges":"false"}})
    opt_job.update_metadata({"class":classtag})
    opt_job.name = name+"_QChem-Opt"
    #freq TBLite
    f_opt_job = freq_job(opt_job.output, method="GFN2-xTB", pressure=0.986923, calc_swaps={"solvent":"octanol"})
    f_opt_job.update_metadata({"class":classtag})
    f_opt_job.name = name+"_TBLite-Freq_afteropt"
    #sp QChem
    s_opt_job = static_job(opt_job.output, charge=0,spin_multiplicity=1, pcm_dielectric="16.93", overwrite_inputs={"rem": {"ecp":"def2-ecp", "resp_charges":"false"}})
    s_opt_job.update_metadata({"class":classtag})
    s_opt_job.name = name+"_QChem-SP_afteropt"
    #workflow
    flow1 = jf.Flow([f_job, f_qc, s_job])
    wf1 = flow_to_workflow(flow1)
    flow2 = jf.Flow([opt_job, f_opt_job, s_opt_job])
    wf2 = flow_to_workflow(flow2)
    if RUN:
        lpad = LaunchPad.auto_load()
        lpad.add_wf(wf1)
        lpad.add_wf(wf2)