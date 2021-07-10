#!/usr/bin/python3
import sys
import json

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __eq__(self, p):
        if p is None or not isinstance(p, Point): return False
        return self.x == p.x and self.y == p.y
    def __ne__(self, p):
        return not self.__eq__(p)
    def __add__(self, p):
        return Point(self.x + p.x, self.y + p.y)
    def __sub__(self, p):
        return Point(self.x - p.x, self.y - p.y)
    def __mul__(self, s):
        return Point(self.x * s, self.y * s)
    def __floordiv__(self, s):
        return Point(self.x // s, self.y // s)

def distance(a, b):
    dx = a.x - b.x
    dy = a.y - b.y
    return dx * dx + dy * dy

def dot(a, b):
    return a.x * b.x + a.y * b.y

def cross(a, b):
    return a.x * b.y - a.y * b.x

def ccw(a, b, c):
    d = Point(b.x - a.x, b.y - a.y)
    e = Point(c.x - a.x, c.y - a.y)
    if cross(d, e) > 0: return  1
    if cross(d, e) < 0: return -1
    return 0

class Segment:
    def __init__(self, a, b):
        self.a = a
        self.b = b

def intersect(a, b):
    if ccw(a.a, a.b, b.a) * ccw(a.a, a.b, b.b) >= 0: return False
    if ccw(b.a, b.b, a.a) * ccw(b.a, b.b, a.b) >= 0: return False
    return True

def touch(s, p):
    dxa = s.a.x - p.x
    dya = s.a.y - p.y
    dxb = p.x - s.b.x
    dyb = p.y - s.b.y
    if dxa * dyb != dxb * dya: return False
    if p.x < min(s.a.x, s.b.x): return False
    if p.x > max(s.a.x, s.b.x): return False
    if p.y < min(s.a.y, s.b.y): return False
    if p.y > max(s.a.y, s.b.y): return False
    return True

def area2(polygon):
    s = 0
    n = len(polygon)
    for i in range(n):
        s += cross(polygon[i], polygon[(i + 1) % n])
    return s

def contains(polygon, p):
    n = len(polygon)
    result = False
    for i in range(n):
        a = polygon[i] - p
        b = polygon[(i + 1) % n] - p
        if a.y > b.y:
            a, b = b, a
        if a.y <= 0 and 0 < b.y:
            if cross(a, b) < 0:
                result = not result
        if cross(a, b) == 0 and dot(a, b) <= 0:
            return True
    return result


class EvaluationResult:
    def __init__(self):
        self.score    = None
        self.messages = []

    def add_message(self, s):
        self.messages.append(s)


def test_length(problem_edges, problem_vertices, solution_vertices, epsilon):
    for e in problem_edges:
        pa = problem_vertices[e[0]]
        pb = problem_vertices[e[1]]
        sa = solution_vertices[e[0]]
        sb = solution_vertices[e[1]]
        pd = distance(pa, pb)
        sd = distance(sa, sb)
        th = epsilon * pd
        v  = 1000000 * sd - 1000000 * pd
        if v < -th or th < v:
            return False
    return True

def test_in_hole(hole_vertices, problem_edges, solution_vertices):
    pose_edges = []
    for e in problem_edges:
        pose_edges.append(Segment(solution_vertices[e[0]], solution_vertices[e[1]]))
    n = len(hole_vertices)
    for p in solution_vertices:
        if not contains(hole_vertices, p):
            return False
    for e in pose_edges:
        if not contains(hole_vertices, (e.a + e.b) // 2):
            return False
        for i in range(n):
            last_ccw = ccw(e.a, e.b, hole_vertices[i])
            if last_ccw != 0:
                offset = i
                break
        for i in range(n):
            ha = hole_vertices[(i + 0 + offset) % n]
            hb = hole_vertices[(i + 1 + offset) % n]
            hc = hole_vertices[(i + 2 + offset) % n]
            if intersect(e, Segment(ha, hb)):
                return False
            if hb != e.a and hb != e.b and touch(e, hb):
                cur_ccw = ccw(e.a, e.b, hc)
                if cur_ccw * last_ccw < 0:
                    return False
                elif cur_ccw != 0:
                    last_ccw = cur_ccw
    return True

def evaluate(problem, solution):
    result = EvaluationResult()

    if len(solution['vertices']) != len(problem['figure']['vertices']):
        result.add_message('The number of vertices must be equal to vertices in figure.')

    pass_dimension = True
    pass_integer   = True
    solution_vertices = []
    for p in solution['vertices']:
        if len(p) != 2:
            pass_dimension = False
            continue
        if type(p[0]) is not int or type(p[1]) is not int:
            pass_integer = False
            continue
        solution_vertices.append(Point(p[0], p[1]) * 2)
    if not pass_dimension:
        result.add_message('The number of values in each points must be 2.')
    if not pass_integer:
        result.add_message('Coordinate must be an integer.')

    if len(result.messages) > 0:
        return result

    problem_vertices = []
    for p in problem['figure']['vertices']:
        problem_vertices.append(Point(p[0], p[1]) * 2)
    problem_edges = problem['figure']['edges']

    hole_vertices = []
    for p in problem['hole']:
        hole_vertices.append(Point(p[0], p[1]) * 2)
    if area2(hole_vertices) < 0:
        hole_vertices.reverse()

    epsilon = problem['epsilon']
    if not test_length(problem_edges, problem_vertices, solution_vertices, epsilon):
        result.add_message('Some edges are too compressed or stretched.')
    if not test_in_hole(hole_vertices, problem_edges, solution_vertices):
        result.add_message('Pose must fit to the hole.')

    if len(result.messages) > 0:
        return result

    score = 0
    for h in hole_vertices:
        best = distance(h, solution_vertices[0])
        for v in solution_vertices:
            best = min(best, distance(h, v))
        score += best
    result.score = score // 4

    return result


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('python3 evaluate.py problem solution')
        sys.exit()
    with open(sys.argv[1], 'r') as f:
        problem = json.load(f)
    with open(sys.argv[2], 'r') as f:
        solution = json.load(f)
    result = evaluate(problem, solution)
    print(result.score)
    print(result.messages)

