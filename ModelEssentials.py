import ast
import csv
import operator as op
import os
import re
from pathlib import Path

import imageio
import moviepy.video.fx.all as vfx
from moviepy.editor import VideoFileClip

from Const import *

DEFAULT_DELIMITER = ','
PATH_POPULATIONS_CSV = 'CSVs/Populations2.csv'
COLS_NO_COMPARTMENT = 3


def read_sir_csv(APathToInputNetwork,
                 AColumns=None,
                 ADelimiter=DEFAULT_DELIMITER):
    return np.loadtxt(APathToInputNetwork, dtype=str, delimiter=ADelimiter, usecols=AColumns)


def print_progress(ACurrent,
                   ATotal):
    print("{:.2f} %".format(ACurrent / ATotal * 100.0))


def rk4_step(AFunction,
             AVariable,
             AFunctionParameter,
             AStepSize):
    k1 = AFunction(AVariable, *AFunctionParameter)
    k2 = AFunction(AVariable + k1 * AStepSize / 2., *AFunctionParameter)
    k3 = AFunction(AVariable + k2 * AStepSize / 2., *AFunctionParameter)
    k4 = AFunction(AVariable + k3 * AStepSize, *AFunctionParameter)
    return AVariable + AStepSize / 6. * (k1 + 2 * (k2 + k3) + k4)


def read_population_csv(APathToCsv,
                        ACityCount,
                        ANumberOfColumns,
                        ADelimiter=','):
    with open(APathToCsv) as csvfile:
        nameList = [[]] * ACityCount
        positionsX = [[]] * ACityCount
        positionsY = [[]] * ACityCount

        i = 0
        CSV = csv.reader(csvfile, delimiter=ADelimiter)

        for row in CSV:
            nameList[i] = row[0]
            positionsX[i] = row[ANumberOfColumns - 2]
            positionsY[i] = row[ANumberOfColumns - 1]
            i = i + 1
    return [nameList, positionsX, positionsY]


# supported operators
operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
             ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg}


def evaluate_string_matrix(AMatrix):
    result = np.zeros(np.shape(AMatrix), dtype=float)
    for r in range(0, len(result)):
        for j in range(0, len(result)):
            result[r][j] = eval_expr(str(AMatrix[r][j]))
    return result


def string_contains(ASearchString, ASubString):
    return np.char.find(ASearchString, ASubString, start=0, end=None)


def eval_expr(AExpression):
    return eval_(ast.parse(AExpression, mode='eval').body)


def eval_(ANode):
    if isinstance(ANode, ast.Num):  # <number>
        return ANode.n
    elif isinstance(ANode, ast.BinOp):  # <left> <operator> <right>
        return operators[type(ANode.op)](eval_(ANode.left), eval_(ANode.right))
    elif isinstance(ANode, ast.UnaryOp):  # <operator> <operand> e.g., -1
        return operators[type(ANode.op)](eval_(ANode.operand))
    else:
        raise TypeError(ANode)


def infection_gradient(ACurrentInfections,
                       ATotalPopulation):
    percent_infected = ACurrentInfections / ATotalPopulation

    FULL_COLOR = 255.0
    ZERO_COLOR = 0.0

    GREEN = FULL_COLOR - percent_infected * FULL_COLOR
    RED = ZERO_COLOR + percent_infected * FULL_COLOR
    BLUE = 0.0
    ALPHA = 1.0

    return [RED / FULL_COLOR, GREEN / FULL_COLOR, BLUE / FULL_COLOR, ALPHA]


# https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
def atoi(AText):
    return int(AText) if AText.isdigit() else AText


def natural_keys(AText):
    return [atoi(c) for c in re.split(r'(\d+)', AText)]


def check_directory_exists(ADirectory,
                           ACreateDir=True):
    p = Path(ADirectory)

    if not p.exists():
        if ACreateDir:
            p.mkdir()
        return p
    else:
        return p


def count_compartments():
    return get_total_columns_in_population() - COLS_NO_COMPARTMENT


def get_total_columns_in_population():
    return len(np.transpose(read_sir_csv(PATH_POPULATIONS_CSV)))


def tuple_of_compartments_cols(ANumberOfCompartments):
    cols = []
    for i in range(1, ANumberOfCompartments + 1):
        cols.append(i)

    return tuple(cols)


def create_video_of_images(APathToImages,
                           AImageExtension=".png",
                           AMovieName="movie",
                           AMovieExtension=".mp4", fps=30,
                           ASpeed=5):
    check_directory_exists(APathToImages)

    filenames = []

    for file in os.listdir(APathToImages + "/"):
        filename = os.fsdecode(file)
        if filename.endswith(AImageExtension):
            filenames.append('{path}/{name}'.format(path=APathToImages, name=filename))
            continue
        else:
            continue

    if len(filenames) == 0:
        raise RuntimeError('The folder does not contain any images that match the given parameters.')

    filenames.sort(key=natural_keys)

    images = []
    for j in range(0, len(filenames)):
        images.append(imageio.imread(filenames[j]))

    imageio.mimsave('{path}/{movieName}.gif'.format(path=APathToImages, movieName=AMovieName), images)
    clip = VideoFileClip("{path}/{movieName}.gif".format(path=APathToImages, movieName=AMovieName))

    clip = clip.set_fps(clip.fps * fps)
    final = clip.fx(vfx.speedx, ASpeed)
    final.write_videofile("{path}/{movieName}{movieExtension}".format(path=APathToImages, movieName=AMovieName,
                                                                      movieExtension=AMovieExtension))



