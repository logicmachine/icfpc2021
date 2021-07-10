#!/use/bin/python3
import json
import sqlite3
import requests

DATABASE = './database.db'
TOKEN = '575ecf40-1834-46ec-9e38-da4124d1fc81'

def main():
    db = sqlite3.connect(DATABASE)
    cur = db.cursor()
    problems = []
    for row in cur.execute('select id, remote_id from problems'):
        if row[1] is not None:
            problems.append(row)
    for problem in problems:
        cur.execute(
            'select score, body from submissions where problem_id=? order by score asc limit 1',
            (problem[0],))
        row = cur.fetchone()
        if row is None:
            continue
        if row[0] >= 9223372036854775807:
            continue
        r = requests.post(
            'https://poses.live/api/problems/' + problem[1] + '/solutions',
            row[1],
            headers={
                'Content-Type': 'application/json',
                'Authorization': 'Bearer ' + TOKEN
            })
        print(r.json())

if __name__ == '__main__':
    main()
