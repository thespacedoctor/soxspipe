-- DY-1513: replace the cross-joining `cal_<type>` views with reads of a real `calibration_matches` table,
-- which `soxspipe.commonutils.calibration_matching.refresh_calibration_matches` fills, and add the indexes the
-- SOF build needs. Applied once to the template `soxspipe/resources/soxspipe.db`; workspaces made from an older
-- template are rebuilt with `soxspipe prep --refresh`.

BEGIN;

CREATE TABLE IF NOT EXISTS "calibration_matches" (
  "cal_type" text NOT NULL,
  "sof" text NOT NULL,
  "file" text,
  "tag" text,
  "filepath" text,
  "complete" int,
  "upstream_status" text,
  "upstream_recipe" text,
  "upstream_recipe_order" integer,
  "downstream_recipe" text,
  "downstream_recipe_order" integer,
  "mjd_diff" real,
  "rank" integer
);
CREATE INDEX IF NOT EXISTS idx_calibration_matches_type_sof ON "calibration_matches" ("cal_type", "sof");

CREATE INDEX IF NOT EXISTS idx_sof_map_base_sof ON "sof_map_base" ("sof");
CREATE INDEX IF NOT EXISTS idx_sof_map_base_filepath ON "sof_map_base" ("filepath");
CREATE INDEX IF NOT EXISTS idx_product_frames_filepath ON "product_frames" ("filepath");

DROP VIEW IF EXISTS "cal_mbias";
CREATE VIEW "cal_mbias" AS
SELECT
	sof,
	file,
	tag,
	filepath,
	complete,
	file AS upstream_file,
	tag AS upstream_tag,
	filepath AS upstream_filepath,
	upstream_recipe,
	upstream_status,
	1 AS upstream_complete,
	upstream_recipe_order,
	sof AS downstream_sof,
	complete AS downstream_complete,
	downstream_recipe,
	downstream_recipe_order,
	mjd_diff,
	rank AS all_downstream_ranked,
	rank AS matched_upstream_ranked
FROM
	calibration_matches
WHERE
	cal_type = 'mbias'
ORDER BY
	upstream_recipe_order,
	sof,
	rank;

DROP VIEW IF EXISTS "cal_mdark";
CREATE VIEW "cal_mdark" AS
SELECT
	sof,
	file,
	tag,
	filepath,
	complete,
	file AS upstream_file,
	tag AS upstream_tag,
	filepath AS upstream_filepath,
	upstream_recipe,
	upstream_status,
	1 AS upstream_complete,
	upstream_recipe_order,
	sof AS downstream_sof,
	complete AS downstream_complete,
	downstream_recipe,
	downstream_recipe_order,
	mjd_diff,
	rank AS all_downstream_ranked,
	rank AS matched_upstream_ranked
FROM
	calibration_matches
WHERE
	cal_type = 'mdark'
ORDER BY
	upstream_recipe_order,
	sof,
	rank;

DROP VIEW IF EXISTS "cal_disp_solution";
CREATE VIEW "cal_disp_solution" AS
SELECT
	sof,
	file,
	tag,
	filepath,
	complete,
	file AS upstream_file,
	tag AS upstream_tag,
	filepath AS upstream_filepath,
	upstream_recipe,
	upstream_status,
	1 AS upstream_complete,
	upstream_recipe_order,
	sof AS downstream_sof,
	complete AS downstream_complete,
	downstream_recipe,
	downstream_recipe_order,
	mjd_diff,
	rank AS all_downstream_ranked,
	rank AS matched_upstream_ranked
FROM
	calibration_matches
WHERE
	cal_type = 'disp_solution'
ORDER BY
	upstream_recipe_order,
	sof,
	rank;

DROP VIEW IF EXISTS "cal_order_centres";
CREATE VIEW "cal_order_centres" AS
SELECT
	sof,
	file,
	tag,
	filepath,
	complete,
	file AS upstream_file,
	tag AS upstream_tag,
	filepath AS upstream_filepath,
	upstream_recipe,
	upstream_status,
	1 AS upstream_complete,
	upstream_recipe_order,
	sof AS downstream_sof,
	complete AS downstream_complete,
	downstream_recipe,
	downstream_recipe_order,
	mjd_diff,
	rank AS all_downstream_ranked,
	rank AS matched_upstream_ranked
FROM
	calibration_matches
WHERE
	cal_type = 'order_centres'
ORDER BY
	upstream_recipe_order,
	sof,
	rank;

DROP VIEW IF EXISTS "cal_mflat";
CREATE VIEW "cal_mflat" AS
SELECT
	sof,
	file,
	tag,
	filepath,
	complete,
	file AS upstream_file,
	tag AS upstream_tag,
	filepath AS upstream_filepath,
	upstream_recipe,
	upstream_status,
	1 AS upstream_complete,
	upstream_recipe_order,
	sof AS downstream_sof,
	complete AS downstream_complete,
	downstream_recipe,
	downstream_recipe_order,
	mjd_diff,
	rank AS all_downstream_ranked,
	rank AS matched_upstream_ranked
FROM
	calibration_matches
WHERE
	cal_type = 'mflat'
ORDER BY
	upstream_recipe_order,
	sof,
	rank;

DROP VIEW IF EXISTS "cal_spat_solution";
CREATE VIEW "cal_spat_solution" AS
SELECT
	sof,
	file,
	tag,
	filepath,
	complete,
	file AS upstream_file,
	tag AS upstream_tag,
	filepath AS upstream_filepath,
	upstream_recipe,
	upstream_status,
	1 AS upstream_complete,
	upstream_recipe_order,
	sof AS downstream_sof,
	complete AS downstream_complete,
	downstream_recipe,
	downstream_recipe_order,
	mjd_diff,
	rank AS all_downstream_ranked,
	rank AS matched_upstream_ranked
FROM
	calibration_matches
WHERE
	cal_type = 'spat_solution'
ORDER BY
	upstream_recipe_order,
	sof,
	rank;

DROP VIEW IF EXISTS "cal_std_flux";
CREATE VIEW "cal_std_flux" AS
SELECT
	sof,
	file,
	tag,
	filepath,
	complete,
	file AS upstream_file,
	tag AS upstream_tag,
	filepath AS upstream_filepath,
	upstream_recipe,
	upstream_status,
	1 AS upstream_complete,
	upstream_recipe_order,
	sof AS downstream_sof,
	complete AS downstream_complete,
	downstream_recipe,
	downstream_recipe_order,
	mjd_diff,
	rank AS all_downstream_ranked,
	rank AS matched_upstream_ranked
FROM
	calibration_matches
WHERE
	cal_type = 'std_flux'
ORDER BY
	upstream_recipe_order,
	sof,
	rank;

COMMIT;
