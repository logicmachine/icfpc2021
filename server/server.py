#!/usr/bin/python3
import os
import json
import sqlite3
from evaluate import evaluate
from flask import Flask, g, request, abort, jsonify

#============================================================================
# Application settings
#============================================================================

DATABASE   = './database.db'
AUTH_TOKEN = 'G5wh0MEeI4ASstafP4Ih'


#============================================================================
# Flask initialization
#============================================================================

app = Flask(__name__)

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
    return db

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()


#============================================================================
# Authentication
#============================================================================

def authenticate():
    auth = request.headers.get('Authorization', None)
    if not auth:
        abort(403, 'forbidden')
    tokens = auth.split()
    if len(tokens) != 2 or tokens[0].lower() != 'bearer' or tokens[1] != AUTH_TOKEN:
        abort(403, 'forbidden')


#============================================================================
# API definitions
#============================================================================

@app.route('/api/problems')
def problem_list():
    authenticate()
    cur = get_db().cursor()
    best_scores = {}
    for row in cur.execute('select problem_id, min(score) from submissions group by problem_id'):
        best_scores[row[0]] = row[1]
    problems = []
    for row in cur.execute('select id, remote_id from problems order by id asc'):
        data = { 'id': row[0], 'name': row[1], 'best_dislikes': None }
        if row[0] in best_scores:
            data['best_dislikes'] = best_scores[row[0]]
        problems.append(data)
    return jsonify(problems)

@app.route('/api/problems/<problem_id>')
def get_problem(problem_id):
    authenticate()
    cur = get_db().cursor()
    cur.execute('select id, remote_id, body from problems where id=?', (problem_id,))
    row = cur.fetchone()
    if row is None:
        abort(404, 'problem not found')
    return jsonify({
        'id': row[0],
        'name': row[1],
        'problem': json.loads(row[2])
    })

@app.route('/api/submissions/<problem_id>')
def get_submissions(problem_id):
    authenticate()
    cur = get_db().cursor()
    submissions = []
    for row in cur.execute('select id, solver, score from submissions where problem_id=?', (problem_id,)):
        submissions.append({
            'id': row[0],
            'solver': row[1],
            'dislikes': row[2]
        })
    return jsonify(submissions)

@app.route('/api/submissions/<problem_id>/<submission_id>')
def get_submission(problem_id, submission_id):
    authenticate()
    cur = get_db().cursor()
    cur.execute(
        'select solver, score, body from submissions where id=? and problem_id=?',
        (submission_id, problem_id))
    row = cur.fetchone()
    if row is None:
        abort(404, 'submission not found')
    return jsonify({
        'solver': row[0],
        'dislikes': row[1],
        'solution': json.loads(row[2])
    })

@app.route('/api/submit/<problem_id>', methods=['POST'])
def submit_solution(problem_id):
    authenticate()
    db  = get_db()
    cur = db.cursor()
    cur.execute('select id, body from problems where id=?', (problem_id,))
    row = cur.fetchone()
    if row is None:
        abort(404, 'problem not found')
    data = json.loads(request.data.decode('utf-8'))
    solver   = data['solver']
    solution = data['solution']
    result   = evaluate(json.loads(row[1]), solution)
    score    = 9223372036854775807 if result.score is None else result.score
    cur.execute(
        'insert into submissions (problem_id, solver, score, body) values (?, ?, ?, ?)',
        (problem_id, solver, score, json.dumps(solution)))
    db.commit()
    return jsonify({
        'score': result.score,
        'messages': result.messages
    })

