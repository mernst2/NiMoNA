from Const import *
from numpy import exp

#Tried to add counter for lockdown and nolockdown
#Tried to add time difference between lockdown and no lockdown
#runtime increased to two years
#gamma value changed to realistic value

isinnolockdown = True


def none(ACompartment, ATimeStep):
    return None

def lockdown():
    lockdown.counter -= 1

lockdown.counter = 440*20

def nolockdown():
    nolockdown.counter += 1

nolockdown.counter=0


def R0_mitigating(ACompartment, ATimeStep, r0=3.35, k=0.25, r_bar=1.8, isinnolockdown=None):
    if (nolockdown.counter > 440*40) and (lockdown.counter > 1):
        #(nolockdown.counter > 20000) and (lockdown.counter > 1):
        #(ACompartment[0][1] / (ACompartment[0][0] + ACompartment[0][1] + ACompartment[0][2]) * 100000 > 150) and (nolockdown.counter > 20000) and (lockdown.counter > 1): #and (R0_mitigating.counter > 7000) :
        R0 = r0 * exp(- k * ATimeStep) + (1 - exp(- k * ATimeStep)) * r_bar
        lockdown()
        print(lockdown.counter)
        isinnolockdown = False


    elif isinnolockdown == False: #Code doesnt run elif because of isinnolockdown?
            nolockdown.counter = 0
            lockdown.counter = 440 * 20
            isinnolockdown = True
            nolockdown()
            R0 = 3.35
            print("TestTestTestTestTestTest") #Doesnt work


    else:

            nolockdown()
            #R0_mitigating.counter += 1
            R0 = 3.35



    return round(R0 * 0.129, 2)

#R0_mitigating.counter = 0
