import numpy as np, os

preset_path = './parameters/'
test_cases = './tests/test_cases.txt'

AU = 149597870700 # [m]

def preset_list(path):
    return [x[:-4] for x in os.listdir(path)]

def load_preset(file):
    f = open(file, 'r')
    i = 0
    preset = {'parent': {}, 'params': {}}

    for line in f:
        if line[0] != '#':
            if i == 0:
                a_in_AU = bool(int(line))
            elif i == 1:
                use_wb = bool(int(line))
            elif i == 2:
                use_L = bool(int(line))
            else:
                new_line = np.array(line.split(","))
                name = new_line[0].strip()

                preset['parent'][name] = new_line[1].strip()

                params = new_line[2:].astype(float)
                if a_in_AU:
                    params[1] *= AU
                if use_wb:
                    params[4] = (params[4] - params[3] + 360) % 360 # ω = ϖ - Ω
                if use_L:
                    params[5] = (params[5] - params[4] - params[3] + 360) % 360 # M = L - ω - Ω

                preset['params'][name] = params

            i += 1

    f.close()
    return preset

def load_test_cases(file):
    f = open(file, 'r')
    test_cases = {}

    for line in f:
        if line[0] != '#':
            new_line = np.array(line.split(","))
            params = new_line[1:].astype(float)
            test_cases[new_line[0].strip()] = params

    f.close()
    return test_cases
