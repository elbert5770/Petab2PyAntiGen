"""
Model builder. Outputs go to antimony_models/{MODEL_NAME}/ and generated/{MODEL_NAME}/.
"""
import AntiGen_paths

from pyantigen.generate.pyantigen import PyAntiGen
from antimony_modules.Histone.methylation_Zheng_PNAS2012 import (
    MethylationGrid_Zheng_PNAS2012, Observables_Zheng_PNAS2012)


def generate_antimony_model(Isotopes=['']):
    MODEL_NAME = AntiGen_paths.MODEL_NAME
    model = PyAntiGen(name=MODEL_NAME, isotopes=Isotopes)
    MethylationGrid_Zheng_PNAS2012(model)
    Observables_Zheng_PNAS2012(model)

    print(f"Reactions generated: {model.counter}")
    print(f"Rules generated: {len(model.rules)}")

    model.generate(__file__, model_name=MODEL_NAME)


def update_antimony_model():
    generate_antimony_model([''])


if __name__ == "__main__":
    update_antimony_model()
