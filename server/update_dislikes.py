import sys
import sqlite3

DATABASE = './database.db'

def main():
    db = sqlite3.connect(DATABASE)
    cur = db.cursor()

    for line in sys.stdin:
        tokens = line.strip().split('\t')
        if len(tokens) == 3 and tokens[0].isdecimal() and tokens[2].isdecimal():
            cur.execute(
                'update problems set min_dislikes=? where remote_id=?',
                (int(tokens[2]), int(tokens[0])))

    db.commit()
    db.close()

if __name__ == '__main__':
    main()
