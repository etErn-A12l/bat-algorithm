import random
import math

POP_NO = 20
DIM = 2
ITERATIONS = 200

class Bat:
    def __init__(self):
        self.pos = [0] * DIM
        self.pulse = 0
        self.vel = [0] * DIM
        self.loud = 0
        self.freq = 0
        self.fit = 0

def initiate(bats):
    for bat in bats:
        for j in range(DIM):
            bat.vel[j] = random.random()
            bat.pos[0] = random.uniform(3.0, 12.1)
            bat.pos[1] = random.uniform(4.1, 5.8)

        bat.pulse = random.random()
        bat.loud = random.uniform(1, 2)
        bat.freq = random.random()
 
def calc_fitness(bats):
    for bat in bats:
        bat.fit = 1000 - (21.5 + bat.pos[0] * math.sin(4 * math.pi * bat.pos[0]) + bat.pos[1] * math.sin(20 * math.pi * bat.pos[1]))

def best_fitness(bats):
    index = 0
    for i in range(len(bats)):
        if bats[i].fit < bats[index].fit:
            index = i
    return index

def adjust_freq(bats, new_bats):
    for i in range(len(bats)):
        new_bats[i].freq = random.random()

def update_velocity(bats, new_bats, index):
    for i in range(len(bats)):
        for j in range(DIM):
            new_bats[i].vel[j] = bats[i].vel[j] + (bats[i].pos[j] - bats[index].pos[j]) * bats[i].freq

def update_position(bats, new_bats):
    for i in range(len(bats)):
        for j in range(DIM):
            new_bats[i].pos[j] = bats[i].pos[j] + bats[i].vel[j]

def generate_local_solution(bats, index):
    avg_loud = sum([bat.loud for bat in bats]) / POP_NO
    eta = 0.6294

    for j in range(DIM):
        bats[index].pos[j] += eta * avg_loud
        if bats[index].pos[j] < 0:
            bats[index].pos[j] = 0
        elif bats[index].pos[j] > 2:
            bats[index].pos[j] = 2

def update_loud_pulse(bats, iteration):
    ALPHA = 0.9
    r0 = 0
    GAMMA = 0.9

    for bat in bats:
        bat.loud *= ALPHA
        
        r0 = bat.pulse * 0.001 / 100
        bat.pulse = r0 * (1 - math.exp(-GAMMA * iteration))

def main():
    x = [Bat() for _ in range(POP_NO)]
    random.seed(2)

    itr = 0
    best_gen = 0

    initiate(x)
    calc_fitness(x)

    best_idx = best_fitness(x)
    global_best_fit = x[best_idx].fit

    while itr < ITERATIONS:
        y = [Bat() for _ in range(POP_NO)]

        print(f"\n\n\n\t\t\t\t======= ITERATION: {itr + 1} FITNESS ======\n\n\n")
        for i in range(POP_NO):
            print(f"  BAT {i + 1}: {1000 - x[i].fit:.3f}  |", end=" ")

        print(f"\n\n  Best Fitness: {1000 - x[best_idx].fit:.3f}  From BAT {best_idx + 1}")

        adjust_freq(x, y)
        update_velocity(x, y, best_idx)
        update_position(x, y)

        calc_fitness(y)

        for i in range(POP_NO):
            r = random.random()

            if x[i].pulse > r:
                generate_local_solution(x, i)
            
            if r < x[i].loud and y[i].fit < x[i].fit:
                update_loud_pulse(x, itr)
        
        for i in range(POP_NO):
            if y[i].fit < x[i].fit:
                x[i].fit = y[i].fit
                x[i].freq = y[i].freq
                
                for j in range(DIM):
                    x[i].pos[j] = y[i].pos[j]
                    x[i].vel[j] = y[i].vel[j]
        
        best_idx = best_fitness(x)

        if x[best_idx].fit < global_best_fit:
            global_best_fit = x[best_idx].fit
            best_gen = itr

        itr += 1

    print(f"\n\nBEST FIT WAS FOUND AT {best_gen + 1} ITERATION WITH {global_best_fit:.6f} FITNESS !\n\n")

if __name__ == "__main__":
    main()
