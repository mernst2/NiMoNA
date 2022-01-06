import numpy
import numpy as np
import matplotlib.pyplot as plt
import siressentials as ess
import time

# Parameters for variable count of compartments
# ONLY CHANGE HERE
COLS_TO_USE = (1, 2, 3)
TOTAL_COLUMNS = 6

# Constants
NOT_CONTAINS = -1
DEFAULT_DELIMITER = ','

# Paths
PATH_ADJACENCY_CSV = 'AdjacencyMuenster.csv'
PATH_POPULATIONS_CSV = 'Populations2.csv'
PATH_SIR_ADJACENCY = 'SIR.csv'
PATH_SIR_PLACEHOLDERS = 'Placeholder.csv'

PLACEHOLDER_COL = 0
VALUE_COL = 1
CATEGORY_COL = 2

CATEGORY_PLACEHOLDER = 'p'
CATEGORY_VARIABLE = 'v'

INDEX_CITY_NAMES = 0
INDEX_POS_X = 1
INDEX_POS_Y = 2

TIME_STEP = 0.1
TEND = 150

TOTAL_STEPS = int(TEND / TIME_STEP)
TIME_STEPS = np.linspace(start=0, stop=TEND, num=TOTAL_STEPS)

ADJACENCY_MUENSTER_CSV = np.loadtxt(PATH_ADJACENCY_CSV, delimiter=DEFAULT_DELIMITER, skiprows=1,
                                    usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

# TODO: Write function to calculate the columns that should be used in POPULATION_CSV depending on the total number of
#       compartments
# When more than 3 compartments are used, the parameter 'usecols' needs to be changed. City positions always
# need to be in the last two columns in the order x, y
POPULATION_CSV = np.loadtxt(PATH_POPULATIONS_CSV, delimiter=DEFAULT_DELIMITER, usecols=COLS_TO_USE)
POPULATION = POPULATION_CSV[:, 0]
NUMBER_OF_COMPARTMENTS = len(POPULATION_CSV[0])

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

Z_ORDER_NETWORK_LINES = 1
Z_ORDER_CITY_POINTS = 2
Z_ORDER_CITY_NAMES = 3


def read_sir_csv(path_to_input_network, columns=None, ADelimiter=DEFAULT_DELIMITER):
    return np.loadtxt(path_to_input_network, dtype=numpy.str, delimiter=ADelimiter, usecols=columns)


def replace_placeholders(adjacency_matrix_str, placeholder_list, value_list, category_list,
                         category=CATEGORY_PLACEHOLDER):
    adjacency_matrix_final = np.copy(adjacency_matrix_str)

    for z in range(0, len(adjacency_matrix_str)):
        for j in range(0, len(adjacency_matrix_str)):
            current_entry = adjacency_matrix_str[z, j]
            for p in range(0, len(placeholder_list)):
                if np.char.find(current_entry, placeholder_list[p], start=0, end=None) != NOT_CONTAINS \
                        and category_list[p] == category:
                    current_entry = np.char.replace(current_entry, placeholder_list[p], str(value_list[p]))

            adjacency_matrix_final[z, j] = current_entry

    return adjacency_matrix_final


def replace_variables(adjacency_matrix_str, placeholder_list, value_list, category_list):
    return replace_placeholders(adjacency_matrix_str, placeholder_list, value_list, category_list,
                                category=CATEGORY_VARIABLE)


def evaluate_string_matrix(matrix):
    result = np.zeros(np.shape(matrix), dtype=np.float)
    for r in range(0, len(result)):
        for j in range(0, len(result)):
            result[r][j] = ess.eval_expr(str(matrix[r][j]))
    return result


# TODO: Make the population stream dependent on the count of compartments
# TODO: Refactor and optimize the code
def sir_as_network(compartment, adjacency_matrix_population):
    working_compartment = np.transpose(compartment)

    sir_adjacency_matrix = read_sir_csv(PATH_SIR_ADJACENCY)
    sir_placeholders = read_sir_csv(PATH_SIR_PLACEHOLDERS, columns=PLACEHOLDER_COL)
    sir_placeholder_values = read_sir_csv(PATH_SIR_PLACEHOLDERS, columns=VALUE_COL)
    sir_placeholder_categories = read_sir_csv(PATH_SIR_PLACEHOLDERS, columns=CATEGORY_COL)

    sir_adjacency_matrix = replace_placeholders(sir_adjacency_matrix, sir_placeholders, sir_placeholder_values,
                                                sir_placeholder_categories)

    number_of_placeholders = 0
    number_of_variables = 0

    for cat in sir_placeholder_categories:
        if cat == CATEGORY_VARIABLE:
            number_of_variables += 1
        elif cat == CATEGORY_PLACEHOLDER:
            number_of_placeholders += 1
        else:
            raise RuntimeError

    total_population = []
    for o in range(0, NUMBER_OF_CITIES):
        total_population.append(sum(working_compartment[o, :]))

    total_population = np.copy(total_population)
    # Vals for compartments needs to be in the same order as the compartments itself,
    # eg "SIR" -> "Val S, Val I, Val R"
    # Sum of compartments are at the end
    value_list_variables = np.eye(NUMBER_OF_CITIES, number_of_placeholders + number_of_variables)

    for h in range(0, NUMBER_OF_CITIES):
        temp = working_compartment[h]
        temp = np.append(temp, total_population[h])
        for u in range(0, number_of_placeholders):
            value_list_variables[h][u] = sir_placeholder_values[u]
        for g in range(number_of_placeholders, number_of_variables + number_of_placeholders):
            value_list_variables[h][g] = temp[g - number_of_placeholders]

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

    result = np.zeros([NUMBER_OF_CITIES, NUMBER_OF_COMPARTMENTS])

    for t in range(0, NUMBER_OF_CITIES):
        adj = sir_adjacency_matrix_evaluated[t]
        result[t] = np.matmul(np.transpose(adj), working_compartment[t])

    result_transposed = np.transpose(result)

    for n in range(0, NUMBER_OF_CITIES):
        for m in range(0, NUMBER_OF_CITIES):
            if m != n:
                for v in range(0, len(result_transposed)):
                    result_transposed[v][n] += adjacency_matrix_population[n, m] / total_population[m] * compartment[v][
                        m] - adjacency_matrix_population[m, n] / total_population[n] * compartment[v][n]
    return result_transposed


Compartment = np.array(np.transpose(POPULATION_CSV))
sum_of_cities = []
time_steps_of_compartments = []

for q in range(0, NUMBER_OF_COMPARTMENTS):
    sum_of_cities.append(np.zeros(np.shape(TIME_STEPS)))
    time_steps_of_compartments.append((np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))))

sum_of_cities = np.array(sum_of_cities)
time_steps_of_compartments = np.array(time_steps_of_compartments)

time_before = time.time()
for i in range(0, TOTAL_STEPS):
    Compartment = ess.rk4_step(sir_as_network, Compartment, [ADJACENCY_MATRIX], TIME_STEP)
    for c in range(0, NUMBER_OF_COMPARTMENTS):
        time_steps_of_compartments[c, i] = Compartment[c]
        sum_of_cities[c, i] = sum(Compartment[c])

    print("Iteration {it} finished".format(it=i))

time_after = time.time()
print("Duration: {dt}s".format(dt=time_after-time_before))

# TODO: Make read_population_csv read the total number of columns in the csv and remove the argument
pop_csv = ess.read_population_csv(PATH_POPULATIONS_CSV, NUMBER_OF_CITIES, TOTAL_COLUMNS)
district_name_list = pop_csv[INDEX_CITY_NAMES]  # to read first row of populations2.csv, to use for plot titles
cities_position_x = pop_csv[INDEX_POS_X]
cities_position_y = pop_csv[INDEX_POS_Y]


def create_network_plot(directory, save_figs=False, only_calc_last=True, plot_figures=False):
    for iterator in range(0, NUMBER_OF_CITIES):
        cities_position_x[iterator] = float(cities_position_x[iterator])
        cities_position_y[iterator] = float(cities_position_y[iterator])

    beginning_iterations = 0

    if only_calc_last:
        beginning_iterations = TOTAL_STEPS - 1

    for j in range(beginning_iterations, TOTAL_STEPS):
        fig, ax = plt.subplots()
        ax.axis('off')

        for n in range(0, NUMBER_OF_CITIES):
            for m in range(1, NUMBER_OF_CITIES):
                ax.plot([cities_position_x[m], cities_position_x[n]], [cities_position_y[m], cities_position_y[n]],
                        color=CONNECTION_COLOR, zorder=Z_ORDER_NETWORK_LINES)

        for iterator in range(0, NUMBER_OF_CITIES):
            ax.scatter(cities_position_x[iterator], cities_position_y[iterator], marker='o',
                       color=ess.infection_gradient(time_steps_of_compartments[1, j, iterator], POPULATION[iterator]),
                       zorder=Z_ORDER_CITY_POINTS)

        for k, txt in enumerate(district_name_list):
            ax.annotate(txt, (cities_position_x[k] + ANNOTATION_X_OFFSET * len(district_name_list[k]),
                              cities_position_y[k] + ANNOTATION_Y_OFFSET), weight='bold', color=ANNOTATION_COLOR,
                              zorder=Z_ORDER_CITY_NAMES)

        if plot_figures:
            plt.show()

        if save_figs:
            ess.check_directory_exists(directory)
            plt.savefig('{dir}/Plot{number}.png'.format(dir=directory, number=j))

        ess.print_progress(j, TOTAL_STEPS)
        plt.clf()


def plot_sir_network():
    for f in range(0, NUMBER_OF_CITIES):
        plt.subplot(SIR_NETWORK_ROWS, SIR_NETWORK_COLUMNS, f + 1)
        for g in range(0, NUMBER_OF_COMPARTMENTS):
            current_compartment = time_steps_of_compartments[g]
            plt.plot(TIME_STEPS, current_compartment[:, f])

        plt.title(district_name_list[f], fontsize=9)

    for y in range(0, SIR_NETWORK_TOTAL - SIR_NETWORK_COLUMNS):  # remove x label for first two rows
        plt.subplot(SIR_NETWORK_ROWS, SIR_NETWORK_COLUMNS, y + 1)
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=True,  # ticks along the bottom edge are on
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off

    plt.subplot(SIR_NETWORK_ROWS, SIR_NETWORK_COLUMNS, SIR_NETWORK_TOTAL)
    for s in sum_of_cities:
        plt.plot(TIME_STEPS, s)

    plt.xlabel('t in days')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.tight_layout()

    plt.savefig('./SIRNetwork.png')
    plt.show()
    plt.clf()


create_network_plot("./Network", save_figs=False, only_calc_last=True, plot_figures=True)
plot_sir_network()
# ess.create_video_of_images("./Network")
