import numpy as np


def parse_matrix(file_path):
    matrix = []
    with open(file_path) as f:
        lines = f.readlines()
        data = [str(x).split("\n")[0] for x in lines]

        for line in data:
            array = []
            values = line.strip(' ').split(' ')
            for value in values:
                array.append(float(value))

            matrix.append(np.array(array))

    return np.array(matrix)
