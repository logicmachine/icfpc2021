#!/usr/bin/python3
import os
import sys
import json
import requests

# ENDPOINT = 'https://icfpc.logicmachine.jp'
ENDPOINT = 'http://localhost:5000'
TOKEN = 'G5wh0MEeI4ASstafP4Ih'

def main():
    if len(sys.argv) < 4:
        print('python3 submit.py solver_name problem_id solution_json')
        sys.exit()
    with open(sys.argv[3], 'r') as f:
        solution = f.read()
    data = {
        'solver': sys.argv[1],
        'solution': json.loads(solution)
    }
    r = requests.post(
        ENDPOINT + '/api/submit/' + sys.argv[2],
        json.dumps(data),
        headers = {
            'Content-Type': 'application/json',
            'Authorization': 'Bearer ' + TOKEN
        })
    print(r.json())

if __name__ == '__main__':
    main()
