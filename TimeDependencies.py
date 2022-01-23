from Const import *
from numpy import exp


def none(ACompartment, ATimeStep):
    return None


def R0_mitigating(ACompartment, ATimeStep, r0=10, k=1 / 50, r_bar=0.5):
    if ACompartment[0][1] / (ACompartment[0][0] + ACompartment[0][1] + ACompartment[0][2]) * 100000 > 150:
        R0 = r0 * exp(- k * ATimeStep) + (1 - exp(- k * ATimeStep)) * r_bar
    else:
        R0 = 10.0

    return round(R0 * 0.035, 2)
