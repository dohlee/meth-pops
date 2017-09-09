import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

CIRCLE_INTERVAL = 1
READ_INTERVAL = 1
MARGIN = 0.3
RADIUS = 0.47
MAX_HEIGHT_IN_INCHES = 8
MAX_WIDTH_IN_INCHES = 4  # This may not be the max width when there are a lot of CpG sites to show.


def position_circle_center(i, j, numRead, numCircle):
    x = (j + 0.5) * CIRCLE_INTERVAL
    y = numRead - (i + 0.5) * READ_INTERVAL
    return (x, y)


def adjust_figure_size(figure, numRead, numCircle):
    if numRead > MAX_HEIGHT_IN_INCHES:
        ratio = (numCircle + 2 * MARGIN) / (numRead + 2 * MARGIN)
        figure.set_size_inches(MAX_HEIGHT_IN_INCHES * ratio, MAX_HEIGHT_IN_INCHES)
    else:
        ratio = (numRead + 2 * MARGIN) / (numCircle + 2 * MARGIN)
        figure.set_size_inches(min(MAX_WIDTH_IN_INCHES, numCircle), min(MAX_WIDTH_IN_INCHES, numCircle) * ratio)


def set_x_y_limits(ax, numRead, numCircle):
    ax.set_ylim([-0.3, numRead + 0.3])
    ax.set_xlim([-0.3, numCircle + 0.3])


def put_line(ax, numCircle, i):
    y = (i + 0.5) * READ_INTERVAL
    ax.plot([-MARGIN, numCircle + MARGIN], [y, y], zorder=-1, color='black', lw=2)


def generate_circle(ax, data, numRead, numCircle, i, j):
    circle = patches.Circle(position_circle_center(i, j, numRead, numCircle),
                            RADIUS,
                            facecolor=['black', 'white'][int(data[i, j])],
                            edgecolor='black')

    return circle


def draw(ax, data, numRead, numCircle):
    circles = []
    for i in range(numRead):
        put_line(ax, numCircle, i)

        for j in range(numCircle):
            # If we don't have information of this CpG site in this read, just skip it.
            if data[i, j] == -1:
                continue
            # Generate a circle and store it temporarily
            circles.append(generate_circle(ax, data, numRead, numCircle, i, j))

    for circle in circles:
        ax.add_patch(circle)


def save(data, title=None, out=None):
    if out is None:
        out = 'methylation_lollipop.txt'

    numRead, numCircle = data.shape

    figure = plt.figure()
    adjust_figure_size(figure, numRead, numCircle)

    if title:
        plt.title(title)

    ax = figure.add_subplot(111)

    set_x_y_limits(ax, numRead, numCircle)
    draw(ax, data, numRead, numCircle)

    plt.axis('off')
    # bbox_inches='tight' ensures that the entire title is shown
    plt.savefig('test.png', bbox_inches='tight')


def show(data, title=None):
    numRead, numCircle = data.shape

    figure = plt.figure()
    adjust_figure_size(figure, numRead, numCircle)

    if title:
        plt.title(title)

    ax = figure.add_subplot(111)

    set_x_y_limits(ax, numRead, numCircle)
    draw(ax, data, numRead, numCircle)

    plt.axis('off')
    plt.show()
