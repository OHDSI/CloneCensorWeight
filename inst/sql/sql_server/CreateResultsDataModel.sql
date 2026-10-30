-- Create the CloneCensorWeight results data model.
--
-- Parameters (SqlRender):
--   @database_schema   Target schema, e.g. 'my_results' (include the dot separator)
--   @table_prefix      Optional prefix for all table names (default '')
--
-- Usage via SqlRender + DatabaseConnector:
--   sql <- SqlRender::loadRenderTranslateSql(
--     "CreateResultsDataModel.sql",
--     packageName    = "CloneCensorWeight",
--     database_schema = "results.",
--     table_prefix   = ""
--   )
--   DatabaseConnector::executeSql(connection, sql)

IF OBJECT_ID('@database_schema@{table_prefix}ccw_analysis', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_analysis;

CREATE TABLE @database_schema@{table_prefix}ccw_analysis (
  analysis_id          VARCHAR(50)  NOT NULL,
  description          VARCHAR(255) NULL,
  created_date         DATE         NOT NULL,
  package_version      VARCHAR(50)  NULL,
  database_id          VARCHAR(100) NOT NULL,
  database_name        VARCHAR(255) NULL,
  vocabulary_version   VARCHAR(100) NULL,
  min_observation_date DATE         NULL,
  max_observation_date DATE         NULL,
  PRIMARY KEY (analysis_id)
);

IF OBJECT_ID('@database_schema@{table_prefix}ccw_cohort_definition', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_cohort_definition;

CREATE TABLE @database_schema@{table_prefix}ccw_cohort_definition (
  analysis_id          VARCHAR(50)  NOT NULL,
  cohort_definition_id INT          NOT NULL,
  cohort_name          VARCHAR(255) NOT NULL,
  PRIMARY KEY (analysis_id, cohort_definition_id)
);

IF OBJECT_ID('@database_schema@{table_prefix}ccw_cohort_count', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_cohort_count;

CREATE TABLE @database_schema@{table_prefix}ccw_cohort_count (
  analysis_id          VARCHAR(50)  NOT NULL,
  cohort_definition_id INT          NOT NULL,
  arm                  VARCHAR(100) NOT NULL,
  n_persons            INT          NOT NULL,
  n_events             INT          NULL,
  PRIMARY KEY (analysis_id, cohort_definition_id, arm)
);

IF OBJECT_ID('@database_schema@{table_prefix}ccw_outcome_model', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_outcome_model;

CREATE TABLE @database_schema@{table_prefix}ccw_outcome_model (
  analysis_id           VARCHAR(50)  NOT NULL,
  model_type            VARCHAR(50)  NOT NULL,
  estimand              VARCHAR(100) NOT NULL,
  variance_estimator    VARCHAR(50)  NOT NULL,
  arm1                  VARCHAR(100) NOT NULL,
  arm2                  VARCHAR(100) NOT NULL,
  estimate              FLOAT        NOT NULL,
  lower_95              FLOAT        NOT NULL,
  upper_95              FLOAT        NOT NULL,
  p_value               FLOAT        NULL,
  log_likelihood        FLOAT        NULL,
  calibrated_p          FLOAT        NULL,
  calibrated_lower_95   FLOAT        NULL,
  calibrated_upper_95   FLOAT        NULL,
  PRIMARY KEY (analysis_id)
);

IF OBJECT_ID('@database_schema@{table_prefix}ccw_standardized_difference', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_standardized_difference;

CREATE TABLE @database_schema@{table_prefix}ccw_standardized_difference (
  analysis_id  VARCHAR(50)  NOT NULL,
  covariate_id INT          NOT NULL,
  weighted     INT          NOT NULL,
  covariate    VARCHAR(255) NOT NULL,
  mean_arm1    FLOAT        NOT NULL,
  sd_arm1      FLOAT        NULL,
  mean_arm2    FLOAT        NOT NULL,
  sd_arm2      FLOAT        NULL,
  smd          FLOAT        NOT NULL,
  abs_smd      FLOAT        NOT NULL,
  PRIMARY KEY (analysis_id, covariate_id, weighted)
);

IF OBJECT_ID('@database_schema@{table_prefix}ccw_kaplan_meier', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_kaplan_meier;

CREATE TABLE @database_schema@{table_prefix}ccw_kaplan_meier (
  analysis_id VARCHAR(50)  NOT NULL,
  arm         VARCHAR(100) NOT NULL,
  time        FLOAT        NOT NULL,
  survival    FLOAT        NOT NULL,
  lower       FLOAT        NULL,
  upper       FLOAT        NULL,
  n_risk      FLOAT        NULL,
  n_event     FLOAT        NOT NULL,
  PRIMARY KEY (analysis_id, arm, time)
);

IF OBJECT_ID('@database_schema@{table_prefix}ccw_power_mdrr', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_power_mdrr;

CREATE TABLE @database_schema@{table_prefix}ccw_power_mdrr (
  analysis_id       VARCHAR(50)  NOT NULL,
  arm1              VARCHAR(100) NOT NULL,
  arm2              VARCHAR(100) NOT NULL,
  alpha             FLOAT        NOT NULL,
  power             FLOAT        NOT NULL,
  n_events_arm1     INT          NOT NULL,
  n_events_arm2     INT          NOT NULL,
  n_events_eff_arm1 FLOAT        NOT NULL,
  n_events_eff_arm2 FLOAT        NOT NULL,
  mdrr              FLOAT        NOT NULL,
  PRIMARY KEY (analysis_id)
);

IF OBJECT_ID('@database_schema@{table_prefix}ccw_negative_control', 'U') IS NOT NULL
  DROP TABLE @database_schema@{table_prefix}ccw_negative_control;

CREATE TABLE @database_schema@{table_prefix}ccw_negative_control (
  analysis_id VARCHAR(50)  NOT NULL,
  control_id  VARCHAR(255) NOT NULL,
  log_rr      FLOAT        NOT NULL,
  se_log_rr   FLOAT        NOT NULL,
  PRIMARY KEY (analysis_id, control_id)
);
