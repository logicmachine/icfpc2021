import sys
import json

def convert(json_str):
    obj = json.loads(json_str)
    lines = []
    hole = obj['hole']
    lines.append('{}'.format(len(hole)))
    for p in hole:
        lines.append('{} {}'.format(p[0], p[1]))
    figure = obj['figure']
    lines.append('{} {}'.format(len(figure['vertices']), len(figure['edges'])))
    for p in figure['vertices']:
        lines.append('{} {}'.format(p[0], p[1]))
    for e in figure['edges']:
        lines.append('{} {}'.format(e[0], e[1]))
    lines.append('{}'.format(obj['epsilon']))
    return '\n'.join(lines)

def main():
    if len(sys.argv) < 2:
        print('Usage: problem2txt.py input')
        sys.exit()
    with open(sys.argv[1], 'r') as f:
        json_str = f.read()
    print(convert(json_str))

if __name__ == '__main__':
    main()
