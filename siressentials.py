import os
import re
import imageio
import csv
from pathlib import Path
from moviepy.editor import VideoFileClip
import moviepy.video.fx.all as vfx


def rk4_step(rhs, variable, function_parameter, step_size):
    k1 = rhs(variable, *function_parameter)
    k2 = rhs(variable + k1 * step_size / 2., *function_parameter)
    k3 = rhs(variable + k2 * step_size / 2., *function_parameter)
    k4 = rhs(variable + k3 * step_size, *function_parameter)
    return variable + step_size / 6. * (k1 + 2 * (k2 + k3) + k4)


def read_population_csv(pathToCSV, cityCount, delimiter=','):
    with open(pathToCSV) as csvfile:
        nameList = [[]] * cityCount
        positionsX = [[]] * cityCount
        positionsY = [[]] * cityCount

        i = 0
        CSV = csv.reader(csvfile, delimiter=delimiter)
        for row in CSV:
            nameList[i] = row[0]
            positionsX[i] = row[3]
            positionsY[i] = row[4]
            i = i + 1
    return [nameList, positionsX, positionsY]


def infection_gradient(current_infections, total_population):
    percentInfected = current_infections / total_population

    FULL_COLOR = 255.0
    ZERO_COLOR = 0.0

    GREEN = FULL_COLOR - percentInfected * FULL_COLOR
    RED = ZERO_COLOR + percentInfected * FULL_COLOR
    BLUE = 0.0
    ALPHA = 1.0

    return [RED / FULL_COLOR, GREEN / FULL_COLOR, BLUE / FULL_COLOR, ALPHA]


# https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def create_video_of_images(pathToImages, imageExtension=".png", movieName="movie", movieExtension=".mp4", fps=30,
                           speed=5):
    p = Path(pathToImages)
    if not p.exists():
        raise RuntimeError('The path to the image folder does not exists. Please enter an existing folder.')

    filenames = []

    for file in os.listdir(pathToImages + "/"):
        filename = os.fsdecode(file)
        if filename.endswith(imageExtension):
            filenames.append('{path}/{name}'.format(path=pathToImages, name=filename))
            continue
        else:
            continue

    if len(filenames) == 0:
        raise RuntimeError('The folder does not contain any images that match the given parameters.')

    filenames.sort(key=natural_keys)

    images = []
    for j in range(0, len(filenames)):
        images.append(imageio.imread(filenames[j]))

    imageio.mimsave('{path}/{movieName}.gif'.format(path=pathToImages, movieName=movieName), images)
    clip = VideoFileClip("{path}/{movieName}.gif".format(path=pathToImages, movieName=movieName))

    clip = clip.set_fps(clip.fps * fps)
    final = clip.fx(vfx.speedx, speed)
    final.write_videofile("{path}/{movieName}{movieExtension}".format(path=pathToImages, movieName=movieName,
                                                                      movieExtension=movieExtension))
