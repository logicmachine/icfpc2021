#!/usr/bin/env python

import os
import sys
import json
import requests

ENDPOINT = "https://poses.live"
API = "/api/problems/{}"
TOKEN = "575ecf40-1834-46ec-9e38-da4124d1fc81"

ALL = range(1, 106+1)

def get_problem(problem_id):
    r = requests.get(
        ENDPOINT + API.format(problem_id),
        headers={
            "Content-Type": "application/json",
            "Authorization": "Bearer {}".format(TOKEN),
        },
    )
    return json.dumps(r.json())


def main(args):
    if len(args) < 2:
        print("python3 {} problem_id".format(args[0]))
        sys.exit()

    problem_id = args[1]
    problem_ids = []
    if problem_id.find("-") >= 0:
        s, e = map(int, problem_id.split("-"))
        problem_ids = range(s, e + 1)
    problem_id = int(problem_id)
    if problem_id > 0:
        problem_ids.append(problem_id)
    else:
        problem_ids = ALL

    for problem_id in problem_ids:
        print("problem_id: {}".format(problem_id))
        problem = get_problem(problem_id)
        with open("{}.problem".format(problem_id), "w") as f:
            f.write(problem)


if __name__ == "__main__":
    main(sys.argv)
