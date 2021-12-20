import numpy as np
import matplotlib.pyplot as plt
import siressentials as ess

# Paths
PATH_ADJACENCY_CSV = 'AdjacencyMuenster.csv'
PATH_POPULATIONS_CSV = 'Populations2.csv'

# Model parameters
INFECTION_RATE = 0.35
REMOVAL_RATE = 0.035

TIME_STEP = 0.1
TEND = 150

TOTAL_STEPS = int(TEND / TIME_STEP)
TIME_STEPS = np.linspace(start=0, stop=TEND, num=TOTAL_STEPS)

ADJACENCY_MUENSTER_CSV = np.loadtxt(PATH_ADJACENCY_CSV, delimiter=',', skiprows=1,
                                    usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))
POPULATION_CSV = np.loadtxt(PATH_POPULATIONS_CSV, delimiter=',', usecols=(1, 2))
POPULATION = POPULATION_CSV[:, 0]
INITIAL_INFECTED = POPULATION_CSV[:, 1]

ADJACENCY_MATRIX = ADJACENCY_MUENSTER_CSV + np.transpose(ADJACENCY_MUENSTER_CSV)
NUMBER_OF_CITIES = len(ADJACENCY_MATRIX)

# SIR subplot parameters
SIR_NETWORK_ROWS = 3
SIR_NETWORK_COLUMNS = 4
SIR_NETWORK_TOTAL = 12

# Parameters for network plot
CONNECTION_COLOR = 'lightgrey'
ANNOTATION_COLOR = 'black'

ANNOTATION_X_OFFSET = -0.01
ANNOTATION_Y_OFFSET = 0.02


def sir(compartment, infection_rate, removal_rate, adjacency_matrix):
    Susceptible, Infected, Removed = compartment
    TotalPopulation = np.zeros(NUMBER_OF_CITIES)
    # are Susceptible, Infected, and Removed arrays? if yes, you can just add them like S + I + R. 
    # if not arrays, make them arrays, numpy addition of arrays is almost always faster than a pure pyhton for loop
    for iterator in range(0, NUMBER_OF_CITIES):
        TotalPopulation[iterator] = SusceptibleN[0, iterator] + InfectedN[0, iterator] + RemovedN[0, iterator]

    dS = -infection_rate * Susceptible * Infected / TotalPopulation
    dI = infection_rate * Susceptible * Infected / TotalPopulation - removal_rate * Infected
    dR = removal_rate * Infected

    for n in range(0, NUMBER_OF_CITIES):
        for m in range(0, NUMBER_OF_CITIES):
            if m != n:
                dS[n] += adjacency_matrix[n, m] / TotalPopulation[m] * Susceptible[m] - adjacency_matrix[m, n] / \
                         TotalPopulation[n] * Susceptible[n]
                dI[n] += adjacency_matrix[n, m] / TotalPopulation[m] * Infected[m] - adjacency_matrix[m, n] / \
                         TotalPopulation[n] * Infected[n]
                dR[n] += adjacency_matrix[n, m] / TotalPopulation[m] * Removed[m] - adjacency_matrix[m, n] / \
                         TotalPopulation[n] * Removed[n]
    return np.array([dS, dI, dR])


SusceptibleN = np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))
SusceptibleN[0, :] = POPULATION
SusceptibleN[0, :] -= INITIAL_INFECTED

InfectedN = np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))
InfectedN[0, :] = INITIAL_INFECTED

RemovedN = np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))

Compartment = [SusceptibleN[0, :], InfectedN[0, :], RemovedN[0, :]]

sumS = np.zeros(np.shape(TIME_STEPS))  # array for total numbers of all cities
sumI = np.zeros(np.shape(TIME_STEPS))
sumR = np.zeros(np.shape(TIME_STEPS))

for i in range(0, TOTAL_STEPS):
    Compartment = ess.rk4_step(sir, Compartment, [INFECTION_RATE, REMOVAL_RATE, ADJACENCY_MATRIX], TIME_STEP)
    SusceptibleN[i], InfectedN[i], RemovedN[i, :] = Compartment  # SusceptibleN[i] works like SusceptibleN[i,:]
    sumS[i] = sum(SusceptibleN[i, :])
    sumI[i] = sum(InfectedN[i, :])
    sumR[i] = sum(RemovedN[i, :])

popCSV = ess.read_population_csv(PATH_POPULATIONS_CSV, NUMBER_OF_CITIES)
DistrictNameList = popCSV[0]  # to read first row of populations2.csv, to use for plot titles
CitiesPositionX = popCSV[1]
CitiesPositionY = popCSV[2]


def create_network_plot(save_figs=False, only_calc_last=True):
    for iterator in range(0, NUMBER_OF_CITIES):
        CitiesPositionX[iterator] = float(CitiesPositionX[iterator])
        CitiesPositionY[iterator] = float(CitiesPositionY[iterator])

    beginningIterations = 0

    if only_calc_last:
        beginningIterations = TOTAL_STEPS - 1

    for j in range(beginningIterations, TOTAL_STEPS):
        fig, ax = plt.subplots()
        ax.axis('off')

        for n in range(0, NUMBER_OF_CITIES):
            for m in range(1, NUMBER_OF_CITIES):
                ax.plot([CitiesPositionX[m], CitiesPositionX[n]], [CitiesPositionY[m], CitiesPositionY[n]],
                        color=CONNECTION_COLOR, zorder=1)

        for iterator in range(0, NUMBER_OF_CITIES):
            ax.scatter(CitiesPositionX[iterator], CitiesPositionY[iterator], marker='o', color=ess.infection_gradient(
                InfectedN[j, iterator], POPULATION[iterator]), zorder=2)

        for k, txt in enumerate(DistrictNameList):
            ax.annotate(txt, (CitiesPositionX[k] + ANNOTATION_X_OFFSET * len(DistrictNameList[k]), CitiesPositionY[k] +
                              ANNOTATION_Y_OFFSET), weight='bold', color=ANNOTATION_COLOR, zorder=3)

        if only_calc_last:
            if j == TOTAL_STEPS-1:
                plt.show()
        else:
            plt.show()

        if save_figs:
            plt.savefig('./Network/Plot{number}.png'.format(number=j))

        print("{:.2f} %".format(j / TOTAL_STEPS * 100.0))
        plt.clf()


def plot_sir_network():
    for c in range(0, NUMBER_OF_CITIES):
        plt.subplot(SIR_NETWORK_ROWS, SIR_NETWORK_COLUMNS, c + 1)
        plt.plot(TIME_STEPS, SusceptibleN[:, c], 'C0-')
        plt.plot(TIME_STEPS, InfectedN[:, c], 'C3-')
        plt.plot(TIME_STEPS, RemovedN[:, c], 'C2-')
        plt.title(DistrictNameList[c], fontsize=9)

    for c in range(0, SIR_NETWORK_TOTAL - SIR_NETWORK_COLUMNS):  # remove x label for first two rows
        plt.subplot(SIR_NETWORK_ROWS, SIR_NETWORK_COLUMNS, c + 1)
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=True,  # ticks along the bottom edge are on
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off

    plt.subplot(SIR_NETWORK_ROWS, SIR_NETWORK_COLUMNS, SIR_NETWORK_TOTAL)
    plt.plot(TIME_STEPS, sumS, 'C0-')
    plt.plot(TIME_STEPS, sumI, 'C3-')
    plt.plot(TIME_STEPS, sumR, 'C2-')
    plt.xlabel('t in days')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.tight_layout()

    plt.savefig('./SIRNetwork.png')
    plt.show()
    plt.clf()


create_network_plot()
plot_sir_network()
ess.create_video_of_images("./Network")
