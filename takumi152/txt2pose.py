import sys
import json

def convert(txt_str):
    str_list = txt_str.split('\n')
    json_obj = {'vertices': []}
    for v in str_list:
        vs = v.split()
        if len(vs) == 2:
            json_obj['vertices'].append([int(x) for x in vs])
    return json.dumps(json_obj)

def main():
    if len(sys.argv) < 2:
        print('Usage: txt2pose.py input')
        sys.exit()
    with open(sys.argv[1], 'r') as f:
        txt_str = f.read()
    print(convert(txt_str))

if __name__ == '__main__':
    main()
