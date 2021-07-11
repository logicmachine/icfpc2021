import sys
import json
from batch_solve import problem2text

def main():
    if len(sys.argv) < 2:
        print('Usage: problem2txt.py input [bonuses...]')
        sys.exit()

    with open(sys.argv[1], 'r') as f:
        json_obj = json.load(f)
    json_obj['globalist'] = False
    json_obj['break_a_leg'] = False
    json_obj['wallhack'] = False
    for b in sys.argv[2:]:
        t = b.lower()
        json_obj[t] = True

    print(problem2text(json.dumps(json_obj)))

if __name__ == '__main__':
    main()
