"""
Model builder. Outputs go to antimony_models/{MODEL_NAME}/ and generated/{MODEL_NAME}/.
"""
import AntiGen_paths

from pyantigen.generate.pyantigen import PyAntiGen
from antimony_modules.Signaling.pkd_pi4k_cert_Weber_BMC2015 import (
    SignalingCascade_Weber_BMC2015, Observables_Weber_BMC2015)


def generate_antimony_model(Isotopes=['']):
    MODEL_NAME = AntiGen_paths.MODEL_NAME
    model = PyAntiGen(name=MODEL_NAME, isotopes=Isotopes)
    SignalingCascade_Weber_BMC2015(model)
    Observables_Weber_BMC2015(model)

    print(f"Reactions generated: {model.counter}")
    print(f"Rules generated: {len(model.rules)}")

    model.generate(__file__, model_name=MODEL_NAME)


def update_antimony_model():
    generate_antimony_model([''])


if __name__ == "__main__":
    update_antimony_model()
