import sys

n_hole = int(sys.stdin.readline())
for i in range(n_hole):
    sys.stdin.readline()

line = sys.stdin.readline()
n_figure, m_figure = map(int, line.split())
for i in range(n_figure):
    print(sys.stdin.readline().strip())
for i in range(m_figure):
    sys.stdin.readline()
sys.stdin.readline()
