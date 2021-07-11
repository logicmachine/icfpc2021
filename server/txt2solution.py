import sys
import json
from batch_solve import text2solution

def main():
    if len(sys.argv) < 2:
        print('Usage: txt2solution.py input')
        sys.exit()
    with open(sys.argv[1], 'r') as f:
        text = f.read()
    print(text2solution(text))

if __name__ == '__main__':
    main()
