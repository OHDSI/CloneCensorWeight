-- The results data model of CloneCensorWeight, as created by createResultsDataModel(). This file is
-- frozen: it creates the first version of the tables. Later changes go only through migrations in
-- inst/sql/sql_server/migrations/ (Migration_<n>-v<version>.sql), which migrateDataModel() applies and
-- records in ccw_migration.
{DEFAULT @table_prefix = ''}

CREATE TABLE @database_schema.@table_prefixccw_analysis (
	analysis_id INT NOT NULL,
	description VARCHAR,
	definition VARCHAR,
	PRIMARY KEY(analysis_id)
);

CREATE TABLE @database_schema.@table_prefixccw_strategy (
	analysis_id INT NOT NULL,
	strategy_id INT NOT NULL,
	label VARCHAR,
	type VARCHAR,
	treatment VARCHAR,
	grace_start INT,
	grace_end INT,
	prohibited_treatments VARCHAR,
	definition VARCHAR,
	PRIMARY KEY(analysis_id, strategy_id)
);

CREATE TABLE @database_schema.@table_prefixccw_target_treatment (
	target_treatment_id BIGINT NOT NULL,
	target_id BIGINT,
	treatment_id BIGINT,
	comparator_treatment_id BIGINT,
	PRIMARY KEY(target_treatment_id)
);

CREATE TABLE @database_schema.@table_prefixccw_target_treatment_outcome (
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	outcome_of_interest INT,
	true_effect_size FLOAT,
	PRIMARY KEY(target_treatment_id, outcome_id)
);

CREATE TABLE @database_schema.@table_prefixccw_attrition (
	sequence_number INT NOT NULL,
	description VARCHAR,
	subjects INT,
	strategy_id INT NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	database_id VARCHAR NOT NULL,
	PRIMARY KEY(sequence_number, strategy_id, analysis_id, target_treatment_id, outcome_id, database_id)
);

CREATE TABLE @database_schema.@table_prefixccw_follow_up_dist (
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	strategy_id INT NOT NULL,
	min_days FLOAT,
	p_10_days FLOAT,
	p_25_days FLOAT,
	median_days FLOAT,
	p_75_days FLOAT,
	p_90_days FLOAT,
	max_days FLOAT,
	database_id VARCHAR NOT NULL,
	PRIMARY KEY(analysis_id, target_treatment_id, outcome_id, strategy_id, database_id)
);

CREATE TABLE @database_schema.@table_prefixccw_covariate (
	covariate_id BIGINT NOT NULL,
	covariate_name VARCHAR,
	analysis_id INT NOT NULL,
	covariate_analysis_id INT,
	database_id VARCHAR NOT NULL,
	PRIMARY KEY(covariate_id, analysis_id, database_id)
);

CREATE TABLE @database_schema.@table_prefixccw_covariate_analysis (
	covariate_analysis_id INT NOT NULL,
	covariate_analysis_name VARCHAR,
	analysis_id INT NOT NULL,
	database_id VARCHAR NOT NULL,
	PRIMARY KEY(covariate_analysis_id, analysis_id, database_id)
);

CREATE TABLE @database_schema.@table_prefixccw_covariate_balance (
	database_id VARCHAR NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	target_strategy_id INT NOT NULL,
	comparator_strategy_id INT NOT NULL,
	covariate_id BIGINT NOT NULL,
	target_mean_before FLOAT,
	comparator_mean_before FLOAT,
	mean_before FLOAT,
	std_diff_before FLOAT,
	std_diff_var_before FLOAT,
	balanced_before INT,
	target_mean_after FLOAT,
	comparator_mean_after FLOAT,
	mean_after FLOAT,
	std_diff_after FLOAT,
	std_diff_var_after FLOAT,
	balanced_after INT,
	target_std_diff FLOAT,
	comparator_std_diff FLOAT,
	target_comparator_std_diff FLOAT,
	PRIMARY KEY(database_id, analysis_id, target_treatment_id, outcome_id, target_strategy_id, comparator_strategy_id, covariate_id)
);

CREATE TABLE @database_schema.@table_prefixccw_censoring_model (
	database_id VARCHAR NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	strategy_id INT NOT NULL,
	model VARCHAR(20) NOT NULL,
	term VARCHAR NOT NULL,
	covariate_id BIGINT,
	coefficient FLOAT,
	PRIMARY KEY(database_id, analysis_id, target_treatment_id, outcome_id, strategy_id, model, term)
);

CREATE TABLE @database_schema.@table_prefixccw_censoring_model_dist (
	database_id VARCHAR NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	strategy_id INT NOT NULL,
	model VARCHAR(20) NOT NULL,
	horizon INT,
	probability FLOAT NOT NULL,
	density FLOAT,
	PRIMARY KEY(database_id, analysis_id, target_treatment_id, outcome_id, strategy_id, model, probability)
);

CREATE TABLE @database_schema.@table_prefixccw_weight_dist (
	database_id VARCHAR NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	strategy_id INT NOT NULL,
	scale VARCHAR(20) NOT NULL,
	truncation VARCHAR(20) NOT NULL,
	band_start INT NOT NULL,
	day_band VARCHAR,
	clone_days BIGINT,
	days INT,
	mean_weight FLOAT,
	min_weight FLOAT,
	p_25_weight FLOAT,
	median_weight FLOAT,
	p_75_weight FLOAT,
	p_95_weight FLOAT,
	p_99_weight FLOAT,
	max_weight FLOAT,
	effective_sample_size FLOAT,
	effective_sample_share FLOAT,
	share_truncated FLOAT,
	cap FLOAT,
	PRIMARY KEY(database_id, analysis_id, target_treatment_id, outcome_id, strategy_id, scale, truncation, band_start)
);

CREATE TABLE @database_schema.@table_prefixccw_cumulative_incidence (
	database_id VARCHAR NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	strategy_id INT NOT NULL,
	time_day INT NOT NULL,
	risk FLOAT,
	risk_lb FLOAT,
	risk_ub FLOAT,
	at_risk_weighted FLOAT,
	at_risk INT,
	PRIMARY KEY(database_id, analysis_id, target_treatment_id, outcome_id, strategy_id, time_day)
);

CREATE TABLE @database_schema.@table_prefixccw_result (
	database_id VARCHAR NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	target_strategy_id INT NOT NULL,
	comparator_strategy_id INT NOT NULL,
	estimand VARCHAR(50) NOT NULL,
	estimate FLOAT,
	ci_95_lb FLOAT,
	ci_95_ub FLOAT,
	p FLOAT,
	one_sided_p FLOAT,
	log_rr FLOAT,
	se_log_rr FLOAT,
	se FLOAT,
	calibrated_rr FLOAT,
	calibrated_ci_95_lb FLOAT,
	calibrated_ci_95_ub FLOAT,
	calibrated_p FLOAT,
	calibrated_one_sided_p FLOAT,
	calibrated_log_rr FLOAT,
	calibrated_se_log_rr FLOAT,
	target_subjects INT,
	comparator_subjects INT,
	target_clones INT,
	comparator_clones INT,
	target_days BIGINT,
	comparator_days BIGINT,
	target_outcomes INT,
	comparator_outcomes INT,
	target_competing_events INT,
	comparator_competing_events INT,
	target_weighted_days FLOAT,
	comparator_weighted_days FLOAT,
	target_weighted_outcomes FLOAT,
	comparator_weighted_outcomes FLOAT,
	PRIMARY KEY(database_id, analysis_id, target_treatment_id, outcome_id, target_strategy_id, comparator_strategy_id, estimand)
);

CREATE TABLE @database_schema.@table_prefixccw_diagnostics_summary (
	database_id VARCHAR NOT NULL,
	analysis_id INT NOT NULL,
	target_treatment_id BIGINT NOT NULL,
	outcome_id BIGINT NOT NULL,
	target_strategy_id INT NOT NULL,
	comparator_strategy_id INT NOT NULL,
	max_sdm FLOAT,
	balance_diagnostic VARCHAR(20),
	mdrr FLOAT,
	mdrr_diagnostic VARCHAR(20),
	ease FLOAT,
	ease_diagnostic VARCHAR(20),
	ess_fraction FLOAT,
	ess_fraction_diagnostic VARCHAR(20),
	max_weight FLOAT,
	max_weight_before_truncation FLOAT,
	max_weight_diagnostic VARCHAR(20),
	positivity FLOAT,
	positivity_share_extreme FLOAT,
	positivity_diagnostic VARCHAR(20),
	unblind INT,
	unblind_for_evidence_synthesis INT,
	PRIMARY KEY(database_id, analysis_id, target_treatment_id, outcome_id, target_strategy_id, comparator_strategy_id)
);

CREATE TABLE @database_schema.@table_prefixccw_migration (
	migration_file VARCHAR PRIMARY KEY,
	migration_order INT NOT NULL UNIQUE
);

