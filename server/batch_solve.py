#!/usr/bin/python3
import re
import json
import argparse
import subprocess
import requests

ENDPOINT = 'https://icfpc.logicmachine.jp'
# ENDPOINT = 'http://localhost:5000'
TOKEN = 'G5wh0MEeI4ASstafP4Ih'

API_HEADERS = {
    'Content-Type': 'application/json',
    'Authorization': 'Bearer ' + TOKEN
}

def problem2text(json_str):
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

def text2solution(text):
    vertices = []
    for line in text.split('\n'):
        if line == '':
            continue
        row = []
        for token in line.split():
            if re.match(r'^\d+$', token):
                row.append(int(token))
            else:
                row.append(token)
        vertices.append(row)
    return { 'vertices': vertices }

def fetch_problem_list():
    r = requests.get(ENDPOINT + '/api/problems', headers=API_HEADERS)
    return r.json()

def fetch_problem(problem_id):
    r = requests.get(ENDPOINT + '/api/problems/' + str(problem_id), headers=API_HEADERS)
    return r.json()

def submit_solution(problem_id, solution):
    r = requests.post(ENDPOINT + '/api/submit/' + str(problem_id), json.dumps(solution), headers=API_HEADERS)
    return r.json()

def run(problem_id, solver_name, args):
    problem = fetch_problem(problem_id)
    problem_body = problem2text(json.dumps(problem['problem']))
    proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    outs, _ = proc.communicate(input=problem_body.encode('utf-8'))
    solution_body = text2solution(outs.decode('utf-8'))
    if len(solution_body['vertices']) == 0:
        # Empty output: skip to submit
        return
    solution = {
        'solver': solver_name,
        'solution': solution_body
    }
    ret = submit_solution(problem_id, solution)
    print(problem_id, ret)


def main():
    parser = argparse.ArgumentParser(description='Solver runner')
    # parser.add_argument('--serial', '-s', action='store_true', help='use serial runner')
    parser.add_argument('--name', '-n', required=True, help='solver name')
    args, remains = parser.parse_known_args()

    problems = fetch_problem_list()
    for problem in problems:
        run(problem['id'], args.name, remains)

if __name__ == '__main__':
    main()
