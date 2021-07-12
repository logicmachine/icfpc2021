create table problems (
			id           integer primary key,
			remote_id    integer not null,
			bonus        integer not null,
			min_dislikes integer not null,
			hole_size    integer not null,
			num_vertices integer not null,
			num_edges    integer not null,
			body         text not null
);
create index problem_remote_index on problems(remote_id, bonus);

create table submissions (
			id         integer   primary key,
			problem_id integer   not null,
			created_at timestamp not null default CURRENT_TIMESTAMP,
			solver     text      not null,
			score      real      not null,
			body       text      not null
);
create index submission_problem_score_index on submissions(problem_id, score);
create index submission_problem_solver_index on submissions(problem_id, solver);
create index submission_problem_solver_score_index on submissions(problem_id, solver, score);

create table bonuses (
			id            integer primary key,
			submission_id integer not null,
			bonus         text    not null,
			problem_id    integer not null
);
create index bonus_submission_index on bonuses(submission_id);
