import time

import matplotlib.pyplot as plt

import TimeDependencies
from Const import *
import ModelEssentials as ess


def replace_placeholders(AModelAdjacencyMatrixStr,
                         APlaceholderList,
                         AValueList,
                         ACategoryList,
                         ACategory):

    adjacencyMatrixFinal = np.copy(AModelAdjacencyMatrixStr)

    for row in range(0, len(AModelAdjacencyMatrixStr)):
        for column in range(0, len(AModelAdjacencyMatrixStr)):
            currentEntry = AModelAdjacencyMatrixStr[row, column]
            for placeholder in range(0, len(APlaceholderList)):
                if ess.string_contains(currentEntry, APlaceholderList[placeholder]) != NOT_CONTAINS \
                        and ACategoryList[placeholder] == ACategory:
                    currentEntry = np.char.replace(currentEntry,
                                                    APlaceholderList[placeholder],
                                                    str(AValueList[placeholder]))

            adjacencyMatrixFinal[row, column] = currentEntry

    return adjacencyMatrixFinal


def replace_variables(AAdjacencyMatrixStr,
                      APlaceholderList,
                      AValueList,
                      ACategoryList):
    return replace_placeholders(AAdjacencyMatrixStr,
                                APlaceholderList,
                                AValueList,
                                ACategoryList,
                                CATEGORY_VARIABLE)


def replace_time_dependencies(AModelAdjacencyMatrixStr,
                              APlaceholderList,
                              AValueList,
                              ACategoryList,
                              ACompartment,
                              ACurrentTimeStep):

    workingValueList = []
    evaluatedValueList = []

    for value in range(0, len(AValueList)):
        if ACategoryList[value] == CATEGORY_TIME_DEPENDANT:
            workingValueList.append(AValueList[value])
        else:
            workingValueList.append(NONE_FUNCTION)

    for value in range(0, len(workingValueList)):
        evaluatedValueList.append(getattr(TimeDependencies,
                                          workingValueList[value])(ACompartment, ACurrentTimeStep))

    return replace_placeholders(AModelAdjacencyMatrixStr,
                                APlaceholderList,
                                np.copy(evaluatedValueList),
                                ACategoryList,
                                CATEGORY_TIME_DEPENDANT)


def compartment_model_core(ACompartment,
                           AModelAdjacencyMatrix,
                           AModelPlaceholders,
                           AModelPlaceholderValues,
                           AModelPlaceholderCategories,
                           AAdjacencyMatrixPopulation,
                           ACurrentTimeStep):

    numberOfPlaceholders = 0
    numberOfVariables = 0
    numberOfTimeDeps = 0

    for category in AModelPlaceholderCategories:
        if category == CATEGORY_VARIABLE:
            numberOfVariables += 1
        elif category == CATEGORY_PLACEHOLDER:
            numberOfPlaceholders += 1
        elif category == CATEGORY_TIME_DEPENDANT:
            numberOfTimeDeps += 1
        else:
            raise RuntimeError

    totalPopulationCities = np.zeros(NUMBER_OF_CITIES)
    valueListVariables = np.eye(NUMBER_OF_CITIES, numberOfPlaceholders + numberOfVariables + numberOfTimeDeps)
    modelAdjacencyMatrixCitiesVals = []
    modelAdjacencyMatrixCities = []
    modelAdjacencyMatrixEvaluated = []
    result = np.zeros([NUMBER_OF_CITIES, NUMBER_OF_COMPARTMENTS])
    workingCompartment = np.transpose(ACompartment)

    for city_ in CITIES:
        totalPopulationCities[city_] = sum(workingCompartment[city_, :])

        temp = workingCompartment[city_]
        temp = np.append(temp, totalPopulationCities[city_])
        for u in range(0, numberOfPlaceholders):
            valueListVariables[city_][u] = AModelPlaceholderValues[u]
        for g in range(numberOfPlaceholders, numberOfVariables + numberOfPlaceholders):
            valueListVariables[city_][g] = temp[g - numberOfPlaceholders]

        modelAdjacencyMatrixCitiesVals.append(
            replace_variables(AModelAdjacencyMatrix,
                                  AModelPlaceholders,
                                  valueListVariables[city_],
                                  AModelPlaceholderCategories))

        modelAdjacencyMatrixCities.append(replace_time_dependencies(modelAdjacencyMatrixCitiesVals[city_],
                                                                        AModelPlaceholders,
                                                                        AModelPlaceholderValues,
                                                                        AModelPlaceholderCategories,
                                                                        workingCompartment,
                                                                        ACurrentTimeStep))

        modelAdjacencyMatrixEvaluated.append(ess.evaluate_string_matrix(np.array(modelAdjacencyMatrixCities[city_])))
        result[city_] = np.matmul(np.transpose(modelAdjacencyMatrixEvaluated[city_]), workingCompartment[city_])

    result = np.transpose(result)

    for row in CITIES:
        for column in CITIES:
            if column != row:
                for comp_ in range(0, len(result)):
                    result[comp_][row] += \
                        AAdjacencyMatrixPopulation[row, column] / totalPopulationCities[column] * \
                        ACompartment[comp_][column] \
                        - AAdjacencyMatrixPopulation[column, row] / totalPopulationCities[row] * \
                        ACompartment[comp_][row]
    return result


Compartment = np.array(np.transpose(POPULATION_CSV))
sumOfCities = []
timeStepsOfCompartments = []

for comp in COMPARTMENTS:
    sumOfCities.append(np.zeros(np.shape(TIME_STEPS)))
    timeStepsOfCompartments.append((np.zeros((len(TIME_STEPS), NUMBER_OF_CITIES))))

sumOfCities = np.array(sumOfCities)
timeStepsOfCompartments = np.array(timeStepsOfCompartments)

modelAdjMatrix = ess.read_sir_csv(PATH_SIR_ADJACENCY)
modelPlaceholders = ess.read_sir_csv(PATH_SIR_PLACEHOLDERS, AColumns=PLACEHOLDER_COL)
modelPlaceholderValues = ess.read_sir_csv(PATH_SIR_PLACEHOLDERS, AColumns=VALUE_COL)
modelPlaceholderCategories = ess.read_sir_csv(PATH_SIR_PLACEHOLDERS, AColumns=CATEGORY_COL)

modelAdjMatrix = replace_placeholders(modelAdjMatrix, modelPlaceholders, modelPlaceholderValues,
                                          modelPlaceholderCategories, CATEGORY_PLACEHOLDER)

timeBefore = time.time()
for iteration in ITERATION_STEPS:
    Compartment = ess.rk4_step(compartment_model_core,
                               Compartment,
                               [modelAdjMatrix,
                                modelPlaceholders,
                                modelPlaceholderValues,
                                modelPlaceholderCategories,
                                ADJACENCY_MATRIX_CITIES, TIME_STEPS[iteration]],
                               TIME_STEP)

    for city in COMPARTMENTS:
        timeStepsOfCompartments[city, iteration] = Compartment[city]
        sumOfCities[city, iteration] = sum(Compartment[city])

timeAfter = time.time()
print("Duration: {dt}s".format(dt=timeAfter - timeBefore))

populationCSV = ess.read_population_csv(PATH_POPULATIONS_CSV,
                                        NUMBER_OF_CITIES,
                                        TOTAL_COLUMNS)

districtNameList = populationCSV[INDEX_CITY_NAMES]  # to read first row of populations2.csv, to use for plot titles
citiesPositionX = populationCSV[INDEX_POS_X]
citiesPositionY = populationCSV[INDEX_POS_Y]


def create_network_plot(ADirectory,
                        ASaveFigs=False,
                        AOnlyCalcLast=True,
                        APlotFigures=False):

    for city_ in CITIES:
        citiesPositionX[city_] = float(citiesPositionX[city_])
        citiesPositionY[city_] = float(citiesPositionY[city_])

    beginningIterations = 0

    if AOnlyCalcLast:
        beginningIterations = TOTAL_STEPS - 1

    for j in range(beginningIterations, TOTAL_STEPS):
        fig, ax = plt.subplots()
        ax.axis('off')

        for n in CITIES:
            for m in range(1, NUMBER_OF_CITIES):
                ax.plot([citiesPositionX[m], citiesPositionX[n]], [citiesPositionY[m], citiesPositionY[n]],
                        color=CONNECTION_COLOR, zorder=Z_ORDER_NETWORK_LINES)

        for city_ in CITIES:
            ax.scatter(citiesPositionX[city_], citiesPositionY[city_], marker='o',
                       color=ess.infection_gradient(timeStepsOfCompartments[1, j, city_], POPULATION[city_]),
                       zorder=Z_ORDER_CITY_POINTS)

        for k, txt in enumerate(districtNameList):
            ax.annotate(txt, (citiesPositionX[k] + ANNOTATION_X_OFFSET * len(districtNameList[k]),
                              citiesPositionY[k] + ANNOTATION_Y_OFFSET), weight='bold', color=ANNOTATION_COLOR,
                        zorder=Z_ORDER_CITY_NAMES)

        if APlotFigures:
            plt.show()

        if ASaveFigs:
            ess.check_directory_exists(ADirectory)
            plt.savefig('{dir}/Plot{number}.png'.format(dir=ADirectory, number=j))

        ess.print_progress(j, TOTAL_STEPS)
        plt.clf()


def plot_calculation_results():
    for city_ in CITIES:
        plt.subplot(MODEL_NETWORK_ROWS, MODEL_NETWORK_COLUMNS, city_ + 1)
        for g in COMPARTMENTS:
            current_compartment = timeStepsOfCompartments[g]
            plt.plot(TIME_STEPS, current_compartment[:, city_])

        plt.title(districtNameList[city_], fontsize=DEFAULT_PLOT_FONTSIZE)

    for y in SUBPLOT_TUPLE:  # remove x label for first two rows
        plt.subplot(MODEL_NETWORK_ROWS, MODEL_NETWORK_COLUMNS, y + 1)
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=True,  # ticks along the bottom edge are on
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off

    plt.subplot(MODEL_NETWORK_ROWS, MODEL_NETWORK_COLUMNS, MODEL_NETWORK_TOTAL)
    for s in sumOfCities:
        plt.plot(TIME_STEPS, s)

    plt.xlabel('t in days')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.tight_layout()

    plt.savefig('./ModelNetwork.png')
    plt.show()
    plt.clf()


# create_network_plot("./Network", save_figs=False, only_calc_last=True, plot_figures=True)
plot_calculation_results()
# ess.create_video_of_images("./Network")
