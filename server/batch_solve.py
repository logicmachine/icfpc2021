#!/usr/bin/python3
import re
import json
import argparse
import subprocess
import requests
from multiprocessing import Pool
from evaluate import evaluate

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
    # Hole
    hole = obj['hole']
    lines.append('{}'.format(len(hole)))
    for p in hole:
        lines.append('{} {}'.format(p[0], p[1]))
    # Figure
    figure = obj['figure']
    lines.append('{} {}'.format(len(figure['vertices']), len(figure['edges'])))
    for p in figure['vertices']:
        lines.append('{} {}'.format(p[0], p[1]))
    for e in figure['edges']:
        lines.append('{} {}'.format(e[0], e[1]))
    # Epsilon
    lines.append('{}'.format(obj['epsilon']))
    # Bonuses
    lines.append('{}'.format(len(obj['bonuses'])))
    for b in obj['bonuses']:
        lines.append('{} {} {} {}'.format(b['bonus'], b['problem'], b['position'][0], b['position'][1]))
    # Flags
    if 'globalist' in obj and obj['globalist']:
        lines.append('1')
        lines.append('GLOBALIST')
    elif 'break_leg' in obj and obj['break_leg']:
        lines.append('1')
        lines.append('BREAK_A_LEG')
    elif 'wallhack' in obj and obj['wallhack']:
        lines.append('1')
        lines.append('WALLHACK')
    elif 'superflex' in obj and obj['superflex']:
        lines.append('1')
        lines.append('SUPERFLEX')
    else:
        lines.append('0')
    return '\n'.join(lines)

def text2solution(text):
    vertices = []
    bonuses = []
    for line in text.split('\n'):
        if line == '':
            continue
        tokens = line.strip().split()
        if tokens[0] == 'GLOBALIST':
            bonuses.append({ 'bonus': tokens[0] })
            continue
        elif tokens[0] == 'BREAK_A_LEG':
            bonuses.append({
                'bonus': tokens[0],
                'edge': [int(tokens[1]), int(tokens[2])]
            })
            continue
        elif tokens[0] == 'WALLHACK':
            bonuses.append({ 'bonus': tokens[0] })
            continue
        row = []
        for token in tokens:
            if re.match(r'^\d+$', token):
                row.append(int(token))
            else:
                row.append(token)
        vertices.append(row)
    return { 'vertices': vertices, 'bonuses': bonuses }

def fetch_problem_list():
    r = requests.get(ENDPOINT + '/api/problems', headers=API_HEADERS)
    return r.json()

def fetch_problem(problem_id):
    r = requests.get(ENDPOINT + '/api/problems/' + str(problem_id), headers=API_HEADERS)
    return r.json()

def submit_solution(problem_id, solution):
    r = requests.post(ENDPOINT + '/api/submit/' + str(problem_id), json.dumps(solution), headers=API_HEADERS)
    return r.json()

def evaluate_solution(problem, solution):
    r = evaluate(problem, solution)
    return {
        'messages': r.messages,
        'score': r.score
    }

def run(problem_id, solver_name, dry_run, args):
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
    if dry_run:
        ret = evaluate_solution(problem['problem'], solution_body)
    else:
        ret = submit_solution(problem_id, solution)
    print(problem_id, ret)

def run_multiprocessing(args):
    run(args['id'], args['name'], args['dry'], args['remains'])


def main():
    parser = argparse.ArgumentParser(description='Solver runner')
    parser.add_argument('--parallel', '-p', action='store_true', help='use parallel runner')
    parser.add_argument('--name', '-n', required=True, help='solver name')
    parser.add_argument('--single', '-s', help='problem id')
    parser.add_argument('--dry', '-d', action='store_true', help='do not submit')
    args, remains = parser.parse_known_args()

    if args.single is not None:
        run(args.single, args.name, args.dry, remains)
    elif args.parallel:
        problems = fetch_problem_list()
        tasks = []
        for p in problems:
            if p['best_dislikes'] != 0:
                tasks.append({
                    'id': p['id'],
                    'name': args.name,
                    'dry': args.dry,
                    'remains': remains
                })
        p = Pool()
        p.map(run_multiprocessing, tasks)
    else:
        problems = fetch_problem_list()
        for problem in problems:
            run(problem['id'], args.name, args.dry, remains)

if __name__ == '__main__':
    main()
