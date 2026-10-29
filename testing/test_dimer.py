from pathlib import Path
import numpy as np
from mmfex import actions

HERE = Path(__file__).parent

first_monomer = 301
last_monomer = 302
total_roots = 5   # states computed in the CASSCF calculation
dim = 5           # states per monomer in the Frenkel Hamiltonian


def test_dimer_matches_reference():
    C, E, TDM = actions.collect(first_monomer=first_monomer,
                                last_monomer=last_monomer,
                                computed_roots=total_roots,
                                wanted_roots=dim,
                                source=str(HERE),
                                custom_names=['301_a', '302_a'])
    F = actions.build(E, C, units="eV")
    w, v, f, tdm = actions.compute(F, TDM, units="eV")

    ref = np.load(HERE / "dimer_reference.npz")
    assert np.allclose(C, ref["C"])
    assert np.allclose(E, ref["E"])
    assert np.allclose(TDM, ref["TDM"])
    assert np.allclose(w, ref["w"], atol=1e-8)
    assert np.allclose(f, ref["f"], atol=1e-8)
    assert np.allclose(np.abs(v), np.abs(ref["v"]), atol=1e-6)   # eigenvector signs are arbitrary