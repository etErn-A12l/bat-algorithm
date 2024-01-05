# @title BAT Algorithm
import time, os, random, math, copy
import numpy as np
from IPython.display import clear_output


def Initiate():
    global bats
    for bat in bats:
        cnt = 0
        for j in range(ROW):
            bat.pos[j] = cnt
            bat.vel[j] = random.random()
            cnt += 1
        random.shuffle(bat.pos)
        bat.pulse = random.randrange(0, 1)
        bat.loud = random.randrange(1, 2)
        bat.freq = random.randrange(0, 2)


def CountDimension():
    global ROW
    with open(filePath, "r") as f:
        lines = f.readlines()

    for line in lines:
        row = [int(x) for x in line.split()]
        for each in row:
            nums.append(each)

    ROW = int(len(nums) ** 0.5)


def MatrixRead():
    global matrix
    a = 0
    for i in range(ROW):
        for j in range(ROW):
            matrix[i][j] = int(nums[a])
            a += 1


def CalFitness(bats):
    for bat in bats:
        bat.fit = int(Fitness(bat.pos))
    return bats


def Fitness(pos):
    fit = 0
    for i in range(ROW - 1):
        fit += matrix[pos[i]][pos[i + 1]]
    fit += matrix[pos[ROW - 1]][pos[0]]
    return fit


def bestFitness(bats):
    index = 0
    for i in range(POP_NO):
        index = i if bats[i].fit < bats[index].fit else index
    return index


def lin_kernighan():
    global bats
    improved = 1
    n = ROW
    for bat in bats:
        while improved:
            improved = 0
            for i in range(1, n - 2):
                if not improved:
                    for j in range(i + 1, n):
                        if not improved:
                            delta = (
                                int(matrix[bat.pos[i - 1]][bat.pos[i]])
                                + matrix[bat.pos[j]][bat.pos[(j + 1) % n]]
                                - matrix[bat.pos[i - 1]][bat.pos[j]]
                                - matrix[bat.pos[i]][bat.pos[(j + 1) % n]]
                            )
                            if delta > 0:
                                while i < j:
                                    bat.pos[i], bat.pos[j] = bat.pos[j], bat.pos[i]
                                    i += 1
                                    j -= 1
                                improved = 0
                        else:
                            break
                else:
                    break

def updtPos():
    global c_bats

    for s in range(POP_NO):
        for t in range(ROW):
            c_bats[s].pos[t] = bats[s].pos[t]

    for bat in c_bats:
        bat.pos = simulated_annealing(bat.pos, 100, 0.99)

        bat.fit = int(Fitness(bat.pos))


def simulated_annealing(tour, initial_temp, cooling_rate):
    n = ROW
    new_tour = tour.copy()

    temp = initial_temp
    while temp > 1:
        i = random.randint(0, n - 1)
        j = random.randint(0, n - 1)
        while i == j:
            j = random.randint(0, n - 1)
        if i > j:
            i, j = j, i

        new_tour[i : j + 1] = reversed(new_tour[i : j + 1])

        current_cost = 0
        for k in range(n):
            current_cost += matrix[tour[k]][tour[(k + 1) % n]]

        new_cost = 0
        for k in range(n):
            new_cost += matrix[new_tour[k]][new_tour[(k + 1) % n]]

        if new_cost < current_cost:
            tour = new_tour.copy()
        else:
            p = math.exp((current_cost - new_cost) / temp)
            if random.random() < p:
                tour = new_tour.copy()
            else:
                new_tour = tour.copy()

        temp *= cooling_rate

    return tour


def update_loud_pulse(i):
    global bats
    ALPHA, GAMMA = 0.9, 0.9
    for bat in bats:
        bat.loud *= ALPHA
        r0 = bat.pulse * 0.001 / 100  # 0.001 % of initial pulse rate
        bat.pulse = r0 * (1 - math.exp(-GAMMA * i))


ROW = 0
POP_NO = 80  # 44-43 # number of bats
DIM = 100  # number of dimensions
ITERATIONS = 20000  # number of iterations
optimum = 38673
ROW = DIM
filePath = "/content/ftv35.txt"
nums = []


class BAT:
    def __init__(self):
        self.pulse = 0.0
        self.vel = [0.0] * DIM
        self.loud = 0.0
        self.freq = 0.0
        self.fit = 0
        self.pos = [0] * ROW


random.seed(time.time())


CountDimension()
matrix = np.zeros((ROW, ROW), dtype=int)
MatrixRead()
# print(matrix)
bats = [BAT() for _ in range(POP_NO)]
Initiate()
bats = CalFitness(bats)
A = random.random()
BestGen = 0

bst_IDX = bestFitness(bats)
# bestBAT = BAT()
bestBAT = bats[bst_IDX]
# for bat in bats:
#     print(bat.pos)
Itr = 1
while True:
    clear_output()

    c_bats = [BAT() for _ in range(POP_NO)]

    # print(f"\n\n\t\t======= FITNESS OF {o + 1} ITERATION =======\n\n")
    # for k in range(POP_NO):
    #     print(f"  BAT {k} : \t{bats[k].fit}")

    print(f"  Best Fitness : {bestBAT.fit} From BAT {bst_IDX+1}", end="")

    updtPos()

    c_bats = CalFitness(c_bats)

    for z, bat in enumerate(bats):
        r = random.random()
        if bat.pulse > r:
            lin_kernighan()
        if r < bat.loud and c_bats[z].fit < bat.fit:
            update_loud_pulse(z)

    bats = CalFitness(bats)

    for s in range(POP_NO):
        if c_bats[s].fit < bats[s].fit:
            bats[s].fit = c_bats[s].fit
            bats[s].pos = c_bats[s].pos.copy()

    bst_IDX = bestFitness(bats)

    if bats[bst_IDX].fit < bestBAT.fit:
        BestGen = Itr
        bestBAT.fit = bats[bst_IDX].fit
        bestBAT.pos = bats[bst_IDX].pos.copy()
    Itr += 1

    if bestBAT.fit == optimum:
        break

print(f"\nBEST FIT WAS FOUND AT {BestGen + 1} ITERATION WITH {bestBAT.fit} FITNESS !")
print("\nThe GNOME WAS\n")
for p in bestBAT.pos:
    print(f"{p} ", end="")
