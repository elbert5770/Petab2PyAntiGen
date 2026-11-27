"""
Model builder. Outputs go to antimony_models/{MODEL_NAME}/ and generated/{MODEL_NAME}/.
"""
import os

import AntiGen_paths

from pyantigen.generate.pyantigen import PyAntiGen
from antimony_modules.Signaling.insulin_ros_foxo_Smith_BMCSystBiol2013 import (
    InsulinSignaling_Smith_BMCSystBiol2013, OxidativeStress_Smith_BMCSystBiol2013, Foxo1_Smith_BMCSystBiol2013,
    Observables_Smith_BMCSystBiol2013)


def generate_antimony_model(Isotopes=['']):
    MODEL_NAME = AntiGen_paths.MODEL_NAME
    model = PyAntiGen(name=MODEL_NAME, isotopes=Isotopes)
    # SMITH_PIP3_SWITCH_EPS=0.1 generates the smoothed Akt/PKC switch (an AMICI-friendly variant, not PEtab's model); unset: PEtab's.
    eps = float(os.environ.get("SMITH_PIP3_SWITCH_EPS") or 0) or None
    InsulinSignaling_Smith_BMCSystBiol2013(model, pip3_switch_eps=eps)
    OxidativeStress_Smith_BMCSystBiol2013(model)
    Foxo1_Smith_BMCSystBiol2013(model)
    Observables_Smith_BMCSystBiol2013(model)

    print(f"Reactions generated: {model.counter}")
    print(f"Rules generated: {len(model.rules)}")

    model.generate(__file__, model_name=MODEL_NAME)


def update_antimony_model():
    generate_antimony_model([''])


if __name__ == "__main__":
    update_antimony_model()
