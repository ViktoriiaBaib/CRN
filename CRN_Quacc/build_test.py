from ase.io import read
import jobflow as jf
from jobflow.managers.fireworks import flow_to_workflow
from fireworks import LaunchPad
from quacc.recipes.qchem.core import relax_job, static_job
from quacc.recipes.tblite.core import freq_job
import os

RUN = True

directory = "/global/home/groups/lr_mp/vbaibakova/bfo_newquacc/xyz"
classtag = "bfo_test"

print(directory)

ch=0

for filename in os.listdir(directory):
    name=filename[:-4]
    f = os.path.join(directory, filename)
    print(f)
    try:
        atoms = read(f, format='xyz')
        print("OK xyz")
    except:
        print("Bad xyz!")
    charges = [0]*len(atoms.arrays['numbers'])
    atoms.set_initial_charges(charges)
    atoms.set_initial_magnetic_moments(charges)
    #opt QChem
    opt_job = relax_job(atoms, charge=0,spin_multiplicity=1,pcm_dielectric="16.93", overwrite_inputs={"rem": {"ecp":"def2-ecp", "resp_charges":"false"}})
    opt_job.update_metadata({"class":classtag})
    opt_job.name = name+"_QChem-Opt"
    #freq TBLite
    f_job = freq_job(opt_job.output["atoms"], method="GFN2-xTB", pressure=0.986923)
    f_job.update_metadata({"class":classtag})
    f_job.name = name+"_TBLite-Freq"
    #sp QChem
    s_job = static_job(opt_job.output["atoms"], charge=0,spin_multiplicity=1, pcm_dielectric="16.93", overwrite_inputs={"rem": {"ecp":"def2-ecp", "resp_charges":"false"}})
    s_job.update_metadata({"class":classtag})
    s_job.name = name+"_QChem-SP"
    #workflow
    flow = jf.Flow([opt_job,f_job, s_job])
    wf = flow_to_workflow(flow)
    if RUN:
        lpad = LaunchPad.auto_load()
        lpad.add_wf(wf)
