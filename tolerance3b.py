"""Settle objects of a tolerance relations in 2D."""

import random
import math
import copy
from fractions import gcd
import matplotlib.pyplot as plt


LIMIT_ATTRACT = 0.5             # distance of similar objects
LIMIT_REPULSE = 25.0            # distance of dissimilar objects
LIMIT_REPULSE2 = 15.0
SCALE = 0.1                     # factor of remission
R_SIZE = 200                    # length of remission process
SIZE_X = 5                      # sizes of the picture
SIZE_Y = 5


def init_points(n):
    """Random starting position.
    n - number of points"""
    return [[random.gauss(0.0, LIMIT_REPULSE/5),
             random.gauss(0.0, LIMIT_REPULSE/5)] for x in range(n)]


def init_points2(n):
    """Objects positioned on a circle.
    n - number of points"""
    r = 2 * LIMIT_REPULSE
    indices = list(range(n))
    random.shuffle(indices)
    points = []
    for i in range(n):
        fi = indices[i]/n * 2 * math.pi
        points.append([r*math.cos(fi), r*math.sin(fi)])
    return points

# ################################################# Linear forces

def linear_force(d, critical):
    """Calculate the force.
    d - distance
    critical - alternating distance"""
    return (critical-d)


def linear_attract(d):
    """Force for similar objects.
    d - distance"""
    return linear_force(d, LIMIT_ATTRACT)


def linear_repulse(d):
    """Force for dissimilar objects.
    d - distance"""
    return linear_force(d, LIMIT_REPULSE)


def linear_repulse2(d):
    """Force for neutral objects.
    d - distance"""
    if d < LIMIT_REPULSE2:
        return linear_force(d, LIMIT_REPULSE2)
    elif d > LIMIT_REPULSE:
        return linear_force(d, LIMIT_REPULSE)
    else:
        return 0.0


def linear_neutral(d):
    """Force for neutral objects. (a variant)
    d - distance"""
    if d < LIMIT_ATTRACT:
        return linear_force(d, LIMIT_ATTRACT)
    elif d > LIMIT_REPULSE:
        return linear_force(d, LIMIT_REPULSE)
    else:
        return 0.0

# ################################################# Atan forces


def atan_force(d, critical):
    """Calculate the force.
    d - distance
    critical - alternate distance"""
    return math.atan(d-critical)


def atan_attract(d):
    """Force for similar objects.
    d - distance"""
    return atan_force(d, LIMIT_ATTRACT)


def atan_repulse(d):
    """Force for dissimilar objects.
    d - distance"""
    return atan_force(d, LIMIT_REPULSE)

########################################### 
def print_points(ps):
    """Positions in readable format
    ps - list of pairs of coordinates."""
    for p in ps:
        x, y = p
        print("({0:4.2f}, {1:4.2f}) ".format(x, y), end=" ")
    print()


def movement(points, relation, attract, repulse, neutral, n):
    """Calculate the superposition of forces.
    points - list of pairs of coordinates
    relation - the tolerance relation
    attract - function for calculate the forces between similar objects
    repulse - function for calculate the forces between dissimilar objects
    neutral - function for calculate the forces between neutral objects
    n - number of objects"""
    
    move = [[0.0, 0.0] for i in range(n)]       # sum of forces
    for i in range(n):
        for j in range(i+1, n):
            dx = points[j][0]-points[i][0]
            dy = points[j][1]-points[i][1]
            d = math.sqrt(dx * dx + dy * dy)    # eucledian distance of objects
            if (i, j) in relation:
                if relation[(i, j)] == 1:
                    force = attract(d)          
                else:
                    force = repulse(d)          
            else:
                force = neutral(d)               
            # add to get the superposition
            move[i][0] -= dx/d ** 1.5 * force
            move[i][1] -= dy/d ** 1.5 * force
            move[j][0] += dx/d ** 1.5 * force
            move[j][1] += dy/d ** 1.5 * force
    return move


def maximal_move(movement):
    """Length of the biggest movement."""
    max = 0.0
    for m in movement:
        d = math.hypot(m[0], m[1])      # length of one vector
        if d > max:                     # search the maximal
            max = d
    return max


def increment(points, movement, scale):
    """Apply the forces on objects.
    points - the position of the objects
    movement - superposition of forces on objects
    """
    new_points = []
    for i, xy in enumerate(points):
        uv = movement[i]
        new_points.append([xy[0]+uv[0]*scale, xy[1]+uv[1]*scale])
    return new_points


def gcd_rel(n):
    """Tolerance relation based on gcd.
    n - number of objects"""
    relation = {}
    for i in range(n):
        for j in range(i+1, n):
            if gcd(i+1, j+1) != 1:
                relation[(i, j)] = 1
                relation[(j, i)] = 1
            else:
                relation[(i, j)] = -1
                relation[(j, i)] = -1
    # print("gcd: ", relation)
    return relation


def load_rel(filename):
    """Load a tolerance relation form a file."""
    relation = {}
    with open(filename) as f:
        lines = f.read().splitlines()
        for i, line in enumerate(lines):
            for j, char in enumerate(line):
                if char == "1":
                    relation[(i, j)] = 1
                elif char == "2":
                    relation[(i, j)] = -1
                else:
                    pass    # 0-t nem tárolunk
    return (len(lines), relation)


def test_rel(relation):
    """Check the relation."""
    good = True
    for ij in relation:
        i, j = ij
        if (j, i) not in relation:
            print("At ({0},{1}) missing value {2}".format(j, i, relation[ij]))
            good = False
        elif relation[(i, j)] != relation[(j, i)]:
            print("Not symmetric at ({0},{1}) {3} and {4}".format(
                i, j, relation[(i, j)], relation[(j, i)]))
            good = False
    return good


def remission(N, attract, repulse, neutral):
    """Calculate the remission process (list of max. changes).
    N - number of objects (numbers)
    attract, repules, neutral - the suitable functions"""
    p = init_points2(N)                                  # random start
    rel = gcd_rel(N)                                     # use the gdc relation
    m = movement(p, rel, attract, repulse, neutral, N)   # one step forward
    mm = maximal_move(m)
    scale = LIMIT_REPULSE / mm / 4
    ms = []
    for i in range(R_SIZE):
        ms.append(mm)                               # store the maximal length
        p = increment(p, m, scale)                  # execute the step
        m = movement(p, rel, attract, repulse, neutral, N)   # one step again
        mm = maximal_move(m)
    return ms


def average_remission(N, count):
    """Draw the average.
    N - number of objects
    count - number of repetition"""
    min = remission(N, linear_attract, linear_repulse2, linear_neutral)
    max = copy.deepcopy(min)
    sum = copy.deepcopy(min)
    for i in range(1, count):
        res = remission(N, linear_attract, linear_repulse2, linear_neutral)
        for j, rj in enumerate(res):
            sum[j] += rj
            if rj < min[j]:
                min[j] = rj
            if rj > max[j]:
                max[j] = rj
    avg = [si/count for si in sum]
    xs = list(range(R_SIZE))
    plt.ylim(0, 80)
    plt.plot(xs, min, 'k--', xs, max, 'k--', xs, avg, 'k')
    plt.savefig('g2b.pdf')


def arrange_from_file(filename):
    """Arrange a relation from a file."""
    (N, rel) = load_rel(filename)
    if test_rel(rel):
        ps = init_points(N)
        m = movement(ps, rel, linear_attract, linear_repulse, linear_neutral, N)
        mm = maximal_move(m)
        scale = LIMIT_REPULSE / mm / 4
        for i in range(500):
            m = movement(ps, rel, linear_attract, linear_repulse, linear_neutral, N)
            ps = increment(ps, m, scale)
        print_points(ps)
        xs, ys = zip(*ps)
        plt.scatter(xs, ys)
        for i, p1 in enumerate(ps):
            for j, p2 in enumerate(ps):
                if (i, j) in rel and rel[(i, j)] == 1:
                    plt.plot([p1[0], p2[0]], [p1[1], p2[1]])
        plt.show()

HEADER = r'''
\documentclass{article}
\usepackage{tikz}
\tikzstyle{cblue}=[circle, draw, thin,fill=cyan!20, scale=0.8]
\begin{document}
\begin{tikzpicture}[auto, thick]
\foreach \place/\x in {'''


def arrange2tex(filename, output):
    """Arrange a relation from a file, and generate a TikZ picture."""
    (N, rel) = load_rel(filename)
    if test_rel(rel):
        with open(output, "w") as out:
            out.write(HEADER)
            # számolás kezdete
            ps = init_points(N)
            m = movement(ps, rel, linear_attract, linear_repulse,
                         linear_neutral, N)
            mm = maximal_move(m)
            scale = LIMIT_REPULSE / mm / 4
            for i in range(R_SIZE):
                ps = increment(ps, m, scale)
                m = movement(ps, rel, linear_attract, linear_repulse,
                             linear_neutral, N)
            # képernyőre zsugorítás adatai
            xs, ys = zip(*ps)
            max_x = max(xs)
            min_x = min(xs)
            max_y = max(ys)
            min_y = min(ys)
            sc = max(max_x-min_x, max_y-min_y)
            # jöhet a kiírás
            for j, p in enumerate(ps[:-1]):
                # print(j, p)
                x = (p[0]-min_x)*SIZE_X/sc
                y = (p[1]-min_y)*SIZE_X/sc
                out.write("{{({0:4.2f},{1:4.2f})/{2}}},".format(x, y, j+1))
            x = (ps[-1][0]-min_x)*SIZE_X/sc
            y = (ps[-1][1]-min_y)*SIZE_X/sc
            out.write("{{({0:4.2f},{1:4.2f})/{2}}}}}\n".format(x, y, len(ps)))
            out.write("\\node[cblue] (a\\x) at \place {\\x};\n")
            # pozitív élek berajzolása
            # for ij in rel:
            #    i, j = ij
            #    if rel[ij] == 1:
            #        out.write("\\path[thin] (a{0}) edge (a{1});\n".format(
            #                  i+1, j+1))
            out.write("\\end{tikzpicture}\n\\end{document}")


if __name__ == "__main__":
    #arrange2tex("gcd.tbl", "gcd.tex")
    arrange_from_file('gcd.tbl')
