import os
import glob
import sqlite3

DATABASE = './database.db'

def main():
    db = sqlite3.connect(DATABASE)
    cur = db.cursor()
    for path in glob.iglob('../problems/json/*.problem'):
        basename = os.path.basename(path)
        name, ext = os.path.splitext(basename)
        with open(path, 'r') as f:
            body = f.read()
        cur.execute('insert into problems (id, remote_id, body) values (?, ?, ?)', (int(name), name, body))
    db.commit()
    db.close()

if __name__ == '__main__':
    main()
