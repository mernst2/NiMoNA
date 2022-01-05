import numpy
import numpy as np
import matplotlib.pyplot as plt
import siressentials as ess

# Paths
PATH_ADJACENCY_CSV = 'AdjacencyMuenster.csv'
PATH_POPULATIONS_CSV = 'Populations2.csv'
PATH_SIR_ADJACENCY = 'SIR.csv'
PATH_SIR_PLACEHOLDERS = 'Placeholder.csv'

PLACEHOLDER_COL = 0
VALUE_COL = 1
CATEGORY_COL = 2

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


def read_sir_csv(path_to_input_network, columns=None, ADelimiter=','):
    return np.loadtxt(path_to_input_network, dtype=numpy.str, delimiter=ADelimiter, usecols=columns)


def replace_placeholders(adjacency_matrix_str, placeholder_list, value_list, category_list, category='p'):
    adjacency_matrix_final = np.copy(adjacency_matrix_str)

    for z in range(0, len(adjacency_matrix_str)):
        for j in range(0, len(adjacency_matrix_str)):
            current_entry = adjacency_matrix_str[z, j]
            for p in range(0, len(placeholder_list)):
                if np.char.find(current_entry, placeholder_list[p], start=0, end=None) != -1 \
                        and category_list[p] == category:
                    current_entry = np.char.replace(current_entry, placeholder_list[p], str(value_list[p]))

            adjacency_matrix_final[z, j] = current_entry

    return adjacency_matrix_final


def replace_variables(adjacency_matrix_str, placeholder_list, value_list, category_list):
    return replace_placeholders(adjacency_matrix_str, placeholder_list, value_list, category_list, category='v')


def evaluate_string_matrix(matrix):
    result = np.zeros(np.shape(matrix), dtype=np.float)
    for r in range(0, len(result)):
        for j in range(0, len(result)):
            result[r][j] = ess.eval_expr(str(matrix[r][j]))
    return result

# TODO: Make the population stream dependent on the count of compartments
# TODO: Refactor and optimize the code
def sir_as_network(compartment, infection_rate, removal_rate, adjacency_matrix_population):
    workingCompartment = np.transpose(compartment)

    sir_adjacency_matrix = read_sir_csv(PATH_SIR_ADJACENCY)
    sir_placeholders = read_sir_csv(PATH_SIR_PLACEHOLDERS, columns=PLACEHOLDER_COL)
    sir_placeholder_values = read_sir_csv(PATH_SIR_PLACEHOLDERS, columns=VALUE_COL)
    sir_placeholder_categories = read_sir_csv(PATH_SIR_PLACEHOLDERS, columns=CATEGORY_COL)

    sir_adjacency_matrix = replace_placeholders(sir_adjacency_matrix, sir_placeholders, sir_placeholder_values,
                                                sir_placeholder_categories)

    count_p = 0
    count_v = 0

    for cat in sir_placeholder_categories:
        if cat == 'v':
            count_v = count_v + 1
        elif cat == 'p':
            count_p = count_p + 1
        else:
            raise RuntimeError

    total_population = []
    for o in range(0, NUMBER_OF_CITIES):
        total_population.append(sum(workingCompartment[o, :]))

    total_population = np.copy(total_population)
    # Vals for compartments needs to be in the same order as the compartments itself,
    # eg "SIR" -> "Val S, Val I, Val R"
    # Sum of Compartments are at the end
    value_list_variables = np.eye(NUMBER_OF_CITIES, count_p + count_v)

    for h in range(0, NUMBER_OF_CITIES):
        temp = workingCompartment[h]
        temp = np.append(temp, total_population[h])
        for u in range(0, count_p):
            value_list_variables[h][u] = sir_placeholder_values[u]
        for g in range(count_p, count_v + count_p):
            value_list_variables[h][g] = temp[g - count_p]

    value_list_variables = np.array(value_list_variables)
    sir_adjacency_matrix_cities = []

    for k in range(0, NUMBER_OF_CITIES):
        val = value_list_variables[k]
        temp = replace_variables(sir_adjacency_matrix, sir_placeholders, val,
                                 sir_placeholder_categories)
        sir_adjacency_matrix_cities.append(temp)

    sir_adjacency_matrix_evaluated = []
    for w in range(0, NUMBER_OF_CITIES):
        sir_adjacency_matrix_evaluated.append(evaluate_string_matrix(np.array(sir_adjacency_matrix_cities[w])))

    sir_adjacency_matrix_evaluated = np.array(sir_adjacency_matrix_evaluated)

    result = np.zeros([11, 3])

    for t in range(0, NUMBER_OF_CITIES):
        adj = sir_adjacency_matrix_evaluated[t]
        result[t] = np.matmul(np.transpose(adj), workingCompartment[t])

    result_transposed = np.transpose(result)

    dS = result_transposed[0]
    dI = result_transposed[1]
    dR = result_transposed[2]

    for n in range(0, NUMBER_OF_CITIES):
        for m in range(0, NUMBER_OF_CITIES):
            if m != n:
                dS[n] += adjacency_matrix_population[n, m] / total_population[m] * compartment[0][m] - \
                         adjacency_matrix_population[m, n] / \
                         total_population[n] * compartment[0][n]
                dI[n] += adjacency_matrix_population[n, m] / total_population[m] * compartment[1][m] - \
                         adjacency_matrix_population[m, n] / \
                         total_population[n] * compartment[1][n]
                dR[n] += adjacency_matrix_population[n, m] / total_population[m] * compartment[2][m] - \
                         adjacency_matrix_population[m, n] / \
                         total_population[n] * compartment[2][n]

    return np.array([dS, dI, dR])


def sir(compartment, infection_rate, removal_rate, adjacency_matrix):
    Susceptible, Infected, Removed = compartment
    TotalPopulation = SusceptibleN[0, :] + InfectedN[0, :] + RemovedN[0, :]

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


# Format of the compartments U[x, y]
# x: Index for compartment value at the x-th time-step
# y: The y-th city
SusceptibleN = np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))
SusceptibleN[0, :] = POPULATION
SusceptibleN[0, :] -= INITIAL_INFECTED

InfectedN = np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))
InfectedN[0, :] = INITIAL_INFECTED

RemovedN = np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))

Compartment = [SusceptibleN[0, :], InfectedN[0, :], RemovedN[0, :]]

# array for total numbers of all cities
sumS = np.zeros(np.shape(TIME_STEPS))
sumI = np.zeros(np.shape(TIME_STEPS))
sumR = np.zeros(np.shape(TIME_STEPS))

for i in range(0, TOTAL_STEPS):
    Compartment = ess.rk4_step(sir_as_network, Compartment, [INFECTION_RATE, REMOVAL_RATE, ADJACENCY_MATRIX], TIME_STEP)
    SusceptibleN[i], InfectedN[i], RemovedN[i, :] = Compartment  # SusceptibleN[i] works like SusceptibleN[i,:]
    sumS[i] = sum(SusceptibleN[i, :])
    sumI[i] = sum(InfectedN[i, :])
    sumR[i] = sum(RemovedN[i, :])
    print("Iteration {it} finished".format(it=i))

popCSV = ess.read_population_csv(PATH_POPULATIONS_CSV, NUMBER_OF_CITIES)
DistrictNameList = popCSV[0]  # to read first row of populations2.csv, to use for plot titles
CitiesPositionX = popCSV[1]
CitiesPositionY = popCSV[2]


def create_network_plot(directory, save_figs=False, only_calc_last=True, plot_figures=False):
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

        if plot_figures:
            plt.show()

        if save_figs:
            ess.check_directory_exists(directory)
            plt.savefig('{dir}/Plot{number}.png'.format(dir=directory, number=j))

        ess.print_progress(j, TOTAL_STEPS)
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


# create_network_plot("./Network", save_figs=True, only_calc_last=False, plot_figures=False)
plot_sir_network()
# ess.create_video_of_images("./Network")
