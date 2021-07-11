import os
import glob
import json
import sqlite3

DATABASE = './database.db'
GLOBALIST_OFFSET = 1000
BREAK_LEG_OFFSET = 2000
WALLHACK_OFFSET  = 3000

def main():
    db = sqlite3.connect(DATABASE)
    cur = db.cursor()
    for path in glob.iglob('../problems/json/*.problem'):
        basename = os.path.basename(path)
        name, ext = os.path.splitext(basename)
        key = int(name)
        with open(path, 'r') as f:
            body = json.load(f)
        hole     = len(body['hole'])
        edges    = len(body['figure']['edges'])
        vertices = len(body['figure']['vertices'])
        body['break_leg'] = False
        body['wallhack']  = False
        # plain
        cur.execute(
            '''
            insert or ignore into problems
                (id, remote_id, bonus, min_dislikes, hole_size, num_edges, num_vertices, body)
                values (?, ?, 0, 0, ?, ?, ?, ?)
            ''',
            (key, key, hole, edges, vertices, json.dumps(body)))
        # globalist
        body['epsilon'] *= len(body['figure']['edges'])
        cur.execute(
            '''
            insert or ignore into problems
                (id, remote_id, bonus, min_dislikes, hole_size, num_edges, num_vertices, body)
                values (?, ?, 1, 0, ?, ?, ?, ?)
            ''',
            (key + GLOBALIST_OFFSET, key, hole, edges, vertices, json.dumps(body)))
        body['epsilon'] //= len(body['figure']['edges'])
        # break-a-leg
        body['break_leg'] = True
        cur.execute(
            '''
            insert or ignore into problems
                (id, remote_id, bonus, min_dislikes, hole_size, num_edges, num_vertices, body)
                values (?, ?, 2, 0, ?, ?, ?, ?)
            ''',
            (key + BREAK_LEG_OFFSET, key, hole, edges, vertices, json.dumps(body)))
        body['break_leg'] = False
        # wallhack
        body['wallhack'] = True
        cur.execute(
            '''
            insert or ignore into problems
                (id, remote_id, bonus, min_dislikes, hole_size, num_edges, num_vertices, body)
                values (?, ?, 3, 0, ?, ?, ?, ?)
            ''',
            (key + WALLHACK_OFFSET, key, hole, edges, vertices, json.dumps(body)))
        body['wallhack'] = False


    db.commit()
    db.close()

if __name__ == '__main__':
    main()
