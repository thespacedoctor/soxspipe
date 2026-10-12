-- The `cal_<type>` view SQL from the template `soxspipe.db` before DY-1513, renamed `legacy_cal_<type>`.
-- Kept only so tests can check the nearest-k matcher reproduces it. Do not use it in the pipeline.

CREATE VIEW "legacy_cal_mbias" AS
WITH closest_upstream_match AS (
	-- THIS PARTITION GETS THE CLOSEST MATCHING upstream CALIBRATION FOR EACH downstream PRODUCT FRAME
	SELECT
		*,
		ROW_NUMBER() OVER (PARTITION BY downstream_sof,
			upstream_recipe ORDER BY all_downstream_ranked) AS matched_upstream_ranked
	FROM (
	SELECT
		upstream.file AS upstream_file,
		upstream. 'eso pro catg' AS upstream_tag,
		upstream.filepath AS upstream_filepath,
		upstream.recipe AS upstream_recipe,
		upstream.status as upstream_status,
		upstream.complete as upstream_complete,
		upstream.recipe_order AS upstream_recipe_order,
		downstream.sof AS downstream_sof,
		downstream.complete as downstream_complete,
		downstream.recipe AS downstream_recipe,
		downstream.recipe_order AS downstream_recipe_order,
		ABS(downstream. `mjd-obs` - upstream. `mjd-obs`) AS mjd_diff,
		ROW_NUMBER() OVER (PARTITION BY downstream.binning,
			downstream.slit,
			downstream.recipe,
			downstream.gain,
			downstream.rospeed,
			downstream. `mjd-obs` ORDER BY ABS(downstream. `mjd-obs` - upstream. `mjd-obs`)) AS all_downstream_ranked
	FROM
		z_unique_product_sof_files downstream
		JOIN product_frames upstream 
			ON downstream."eso seq arm" = upstream."eso seq arm"
			AND downstream.recipe_order > upstream.recipe_order
			AND(downstream.rospeed = upstream.rospeed)
			AND((downstream.binning = upstream.binning and downstream.gain = upstream.gain) OR downstream. "eso seq arm" = 'NIR')
				WHERE
		upstream.complete = 1 and downstream_complete < 1  and (upstream.status = 'pass' or upstream.status is null)  AND upstream_recipe = 'mbias')
)
SELECT
	downstream_sof as sof,
	upstream_file as file,
	upstream_tag AS tag,
	upstream_filepath as filepath,
	downstream_complete as complete,
	*
FROM
	closest_upstream_match
WHERE (matched_upstream_ranked = 1)
ORDER BY
	downstream_recipe_order,
	'eso seq arm'
/* cal_mbias(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mbias(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mbias(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mbias(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mbias(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mbias(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */;

CREATE VIEW "legacy_cal_mdark" AS
WITH closest_upstream_match AS (
	-- THIS PARTITION GETS THE CLOSEST MATCHING upstream CALIBRATION FOR EACH downstream PRODUCT FRAME
	SELECT
		*,
		ROW_NUMBER() OVER (PARTITION BY downstream_sof,
			upstream_recipe ORDER BY all_downstream_ranked) AS matched_upstream_ranked
	FROM (
	SELECT
        upstream.file AS upstream_file,
		upstream. 'eso pro catg' AS upstream_tag,
		upstream.filepath AS upstream_filepath,
		upstream.recipe AS upstream_recipe,
		upstream.status as upstream_status,
		upstream.complete as upstream_complete,
		upstream.recipe_order AS upstream_recipe_order,
		downstream.sof AS downstream_sof,
		downstream.complete as downstream_complete,
		downstream.recipe AS downstream_recipe,
		downstream.recipe_order AS downstream_recipe_order,
		ABS(downstream. `mjd-obs` - upstream. `mjd-obs`) AS mjd_diff,
		ABS(downstream. `exptime` - upstream. `exptime`) AS exp_diff,
		ROW_NUMBER() OVER (PARTITION BY downstream.binning,
			downstream.slit,
			downstream.recipe,
			downstream.rospeed,
			downstream. `mjd-obs` ORDER BY ABS(downstream. `exptime` - upstream. `exptime`) , ABS(downstream. `mjd-obs` - upstream. `mjd-obs`)) AS all_downstream_ranked
	FROM
		z_unique_product_sof_files downstream
		JOIN product_frames upstream 
			ON downstream."eso seq arm" = upstream."eso seq arm"
			AND downstream.recipe_order > upstream.recipe_order
			AND(downstream.rospeed = upstream.rospeed)
				WHERE
		upstream.complete = 1 and  downstream_complete < 1 and (upstream.status = 'pass' or upstream.status is null)  AND upstream_recipe = 'mdark' and upstream."eso seq arm" = 'NIR')
)
SELECT
	downstream_sof as sof,
	upstream_file as file,
	upstream_tag AS tag,
	upstream_filepath as filepath,
	downstream_complete as complete,
	*
FROM
	closest_upstream_match
WHERE (matched_upstream_ranked = 1 )
ORDER BY
	downstream_recipe_order,
	'eso seq arm'
/* cal_mdark(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,exp_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mdark(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,exp_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mdark(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,exp_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mdark(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,exp_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mdark(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,exp_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mdark(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,exp_diff,all_downstream_ranked,matched_upstream_ranked) */;

CREATE VIEW "legacy_cal_disp_solution" as 
WITH closest_upstream_match AS (
	-- THIS PARTITION GETS THE CLOSEST MATCHING upstream CALIBRATION FOR EACH downstream PRODUCT FRAME
	SELECT
		*,
		ROW_NUMBER() OVER (PARTITION BY downstream_sof,
			upstream_recipe ORDER BY all_downstream_ranked) AS matched_upstream_ranked
	FROM (
	SELECT
		upstream.file AS upstream_file,
		upstream. 'eso pro catg' AS upstream_tag,
		upstream.filepath AS upstream_filepath,
		upstream.recipe AS upstream_recipe,
		upstream.status as upstream_status,
		upstream.complete as upstream_complete,
		upstream.recipe_order AS upstream_recipe_order,
		downstream.sof AS downstream_sof,
		downstream.complete as downstream_complete,
		downstream.recipe AS downstream_recipe,
		downstream.recipe_order AS downstream_recipe_order,		
		ABS(downstream. `mjd-obs` - upstream. `mjd-obs`) AS mjd_diff,
		ROW_NUMBER() OVER (PARTITION BY downstream.binning,
			downstream.slit,
			downstream.recipe,
			downstream.rospeed,
			downstream. `mjd-obs` ORDER BY ABS(downstream. `mjd-obs` - upstream. `mjd-obs`)) AS all_downstream_ranked
	FROM
		z_unique_product_sof_files downstream
		JOIN product_frames upstream ON downstream. "eso seq arm" = upstream. "eso seq arm"
			AND downstream.recipe_order > upstream.recipe_order 
	WHERE
		upstream.complete = 1 and  downstream_complete < 1 and (upstream.status = 'pass' or upstream.status is null)  and upstream_recipe = 'disp_solution' and downstream_recipe in ('order_centres', 'spat_solution'))
)
SELECT
	downstream_sof as sof,
	upstream_file as file,
	upstream_tag AS tag,
	upstream_filepath as filepath,
	downstream_complete as complete,
	*
FROM
	closest_upstream_match
WHERE (matched_upstream_ranked = 1)
ORDER BY
	upstream_recipe_order,
	'eso seq arm'
/* cal_disp_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_disp_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_disp_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_disp_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_disp_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_disp_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */;

CREATE VIEW "legacy_cal_order_centres" AS
WITH closest_upstream_match AS (
	-- THIS PARTITION GETS THE CLOSEST MATCHING upstream CALIBRATION FOR EACH downstream PRODUCT FRAME
	SELECT
		*,
		ROW_NUMBER() OVER (PARTITION BY downstream_sof,
			'eso pro catg',
			'eso pro type',
			'eso pro tech',
			upstream_recipe ORDER BY all_downstream_ranked) AS matched_upstream_ranked
	FROM (
	SELECT
		upstream.file AS upstream_file,
		upstream. 'eso pro catg' AS upstream_tag,
		upstream.filepath AS upstream_filepath,
		upstream.recipe AS upstream_recipe,
		upstream.status as upstream_status,
		upstream.complete as upstream_complete,
		upstream.recipe_order AS upstream_recipe_order,
		downstream.sof AS downstream_sof,
		downstream.complete as downstream_complete,
		downstream.recipe AS downstream_recipe,
		downstream.recipe_order AS downstream_recipe_order,
		ABS(downstream. `mjd-obs` - upstream. `mjd-obs`) AS mjd_diff,
		ROW_NUMBER() OVER (PARTITION BY downstream.binning,
			downstream.slit,
			downstream.recipe,
			downstream.rospeed,
			downstream. `mjd-obs` ORDER BY ABS(downstream. `mjd-obs` - upstream. `mjd-obs`)) AS all_downstream_ranked
	FROM
		z_unique_product_sof_files downstream
		JOIN product_frames upstream ON downstream. "eso seq arm" = upstream. "eso seq arm"
			AND downstream.recipe_order > upstream.recipe_order
	WHERE
		upstream.complete = 1  and downstream_complete < 1  and (upstream.status = 'pass' or upstream.status is null) and upstream_recipe = 'order_centres' and downstream_recipe = 'mflat')
)
SELECT
	downstream_sof as sof,
	upstream_file as file,
	upstream_tag AS tag,
	upstream_filepath as filepath,
	downstream_complete as complete,
	*
FROM
	closest_upstream_match
WHERE (matched_upstream_ranked = 1)
ORDER BY
	upstream_recipe_order,
	'eso seq arm'
/* cal_order_centres(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_order_centres(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_order_centres(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_order_centres(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_order_centres(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_order_centres(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */;

CREATE VIEW "legacy_cal_mflat" AS WITH closest_upstream_match AS (
	-- THIS PARTITION GETS THE CLOSEST MATCHING upstream CALIBRATION FOR EACH downstream PRODUCT FRAME
	SELECT
		*,
		ROW_NUMBER() OVER (PARTITION BY downstream_sof,
			upstream_recipe ORDER BY all_downstream_ranked) AS matched_upstream_ranked
	FROM (
	SELECT
		upstream.file AS upstream_file,
		upstream. 'eso pro catg' AS upstream_tag,
		upstream.filepath AS upstream_filepath,
		upstream.recipe AS upstream_recipe,
		upstream.status as upstream_status,
		upstream.complete as upstream_complete,
		upstream.recipe_order AS upstream_recipe_order,
		downstream.sof AS downstream_sof,
		downstream.complete as downstream_complete,
		downstream.recipe AS downstream_recipe,
		downstream.recipe_order AS downstream_recipe_order,
		ABS(downstream. `mjd-obs` - upstream. `mjd-obs`) AS mjd_diff,
		ROW_NUMBER() OVER (PARTITION BY downstream.binning,
			downstream.slit,
			downstream.recipe,
			downstream.gain,
			downstream.rospeed,
			downstream. `mjd-obs` ORDER BY ABS(downstream. `mjd-obs` - upstream. `mjd-obs`)) AS all_downstream_ranked
	FROM
		z_unique_product_sof_files downstream
		JOIN product_frames upstream ON downstream. "eso seq arm" = upstream. "eso seq arm"
			AND(downstream.rospeed = upstream.rospeed
			OR downstream_recipe in('spat_solution'))
			AND(downstream.slit = upstream.slit
			OR downstream_recipe in('spat_solution')
		OR(downstream. "eso seq arm" = 'NIR'
		AND upstream.slit = 'SLIT1.5'
		AND downstream.slit = 'SLIT5.0'))
			AND((downstream.binning = upstream.binning and downstream.gain = upstream.gain)
			OR downstream. "eso seq arm" = 'NIR')
			AND downstream.recipe_order > upstream.recipe_order
	WHERE
		upstream.complete = 1 and downstream_complete < 1 and (upstream.status = 'pass' or upstream.status is null) 
		AND upstream_recipe = 'mflat')
)
SELECT
	downstream_sof AS sof, upstream_file AS file, upstream_tag AS tag, upstream_filepath AS filepath, downstream_complete AS complete, upstream_recipe, downstream_recipe, *
FROM
	closest_upstream_match
WHERE (matched_upstream_ranked < 3)
and(downstream_recipe != 'spat_solution'
	OR tag LIKE 'ORDER%')
ORDER BY
	upstream_recipe_order,
	'eso seq arm'
/* cal_mflat(sof,file,tag,filepath,complete,upstream_recipe,downstream_recipe,upstream_file,upstream_tag,upstream_filepath,"upstream_recipe:1",upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,"downstream_recipe:1",downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mflat(sof,file,tag,filepath,complete,upstream_recipe,downstream_recipe,upstream_file,upstream_tag,upstream_filepath,"upstream_recipe:1",upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,"downstream_recipe:1",downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mflat(sof,file,tag,filepath,complete,upstream_recipe,downstream_recipe,upstream_file,upstream_tag,upstream_filepath,"upstream_recipe:1",upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,"downstream_recipe:1",downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mflat(sof,file,tag,filepath,complete,upstream_recipe,downstream_recipe,upstream_file,upstream_tag,upstream_filepath,"upstream_recipe:1",upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,"downstream_recipe:1",downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mflat(sof,file,tag,filepath,complete,upstream_recipe,downstream_recipe,upstream_file,upstream_tag,upstream_filepath,"upstream_recipe:1",upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,"downstream_recipe:1",downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_mflat(sof,file,tag,filepath,complete,upstream_recipe,downstream_recipe,upstream_file,upstream_tag,upstream_filepath,"upstream_recipe:1",upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,"downstream_recipe:1",downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */;

CREATE VIEW "legacy_cal_spat_solution" AS
WITH closest_upstream_match AS (
	-- THIS PARTITION GETS THE CLOSEST MATCHING upstream CALIBRATION FOR EACH downstream PRODUCT FRAME
	SELECT
		*,
		ROW_NUMBER() OVER (PARTITION BY downstream_sof,
			upstream_recipe ORDER BY all_downstream_ranked) AS matched_upstream_ranked
	FROM (
	SELECT
		upstream.file AS upstream_file,
		upstream. 'eso pro catg' AS upstream_tag,
		upstream.filepath AS upstream_filepath,
		upstream.recipe AS upstream_recipe,
		upstream.status as upstream_status,
		upstream.complete as upstream_complete,
		upstream.recipe_order AS upstream_recipe_order,
		downstream.sof AS downstream_sof,
		downstream.complete as downstream_complete,
		downstream.recipe AS downstream_recipe,
		downstream.recipe_order AS downstream_recipe_order,
		ABS(downstream. `mjd-obs` - upstream. `mjd-obs`) AS mjd_diff,
		ROW_NUMBER() OVER (PARTITION BY downstream.binning,
			downstream.slit,
			downstream.recipe,
			downstream.rospeed,
			downstream. `mjd-obs` ORDER BY ABS(downstream. `mjd-obs` - upstream. `mjd-obs`)) AS all_downstream_ranked
	FROM
		z_unique_product_sof_files downstream
		JOIN product_frames upstream ON downstream. "eso seq arm" = upstream. "eso seq arm"
			AND downstream.recipe_order > upstream.recipe_order
	WHERE
		upstream.complete = 1 and downstream_complete < 1 and (upstream.status = 'pass' or upstream.status is null)  and upstream_recipe = 'spat_solution')
)
SELECT
	downstream_sof as sof,
	upstream_file as file,
	upstream_tag AS tag,
	upstream_filepath as filepath,
	downstream_complete as complete,
	*
FROM
	closest_upstream_match
WHERE (matched_upstream_ranked<3)
ORDER BY
	upstream_recipe_order,
	'eso seq arm'
/* cal_spat_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_spat_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_spat_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_spat_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_spat_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */
/* cal_spat_solution(sof,file,tag,filepath,complete,upstream_file,upstream_tag,upstream_filepath,upstream_recipe,upstream_status,upstream_complete,upstream_recipe_order,downstream_sof,downstream_complete,downstream_recipe,downstream_recipe_order,mjd_diff,all_downstream_ranked,matched_upstream_ranked) */;

CREATE VIEW "legacy_cal_std_flux" AS

WITH closest_upstream_match AS (
	-- THIS PARTITION GETS THE CLOSEST MATCHING upstream CALIBRATION FOR EACH downstream PRODUCT FRAME
	SELECT
		*,
		ROW_NUMBER() OVER (PARTITION BY downstream_sof ORDER BY all_downstream_ranked) AS matched_upstream_ranked
	FROM (
	SELECT
		upstream.file AS upstream_file,
		upstream. 'eso pro catg' AS upstream_tag,
		upstream.filepath AS upstream_filepath,
		upstream.recipe AS upstream_recipe,
		upstream.status as upstream_status,
		upstream.complete as upstream_complete,
		upstream.recipe_order AS upstream_recipe_order,
		downstream.sof AS downstream_sof,
		downstream.complete as downstream_complete,
		downstream.recipe AS downstream_recipe,
		downstream.recipe_order AS downstream_recipe_order,
		ABS(downstream. `mjd-obs` - upstream. `mjd-obs`) AS mjd_diff,
		ROW_NUMBER() OVER (PARTITION BY downstream.binning,
			downstream.slit,
			downstream.recipe,
			downstream.rospeed,
			downstream. `mjd-obs` ORDER BY ABS(downstream. `mjd-obs` - upstream. `mjd-obs`)) AS all_downstream_ranked
	FROM
		z_unique_product_sof_files downstream
		JOIN product_frames upstream ON downstream. "eso seq arm" = upstream. "eso seq arm"
			AND (upstream.slit = 'SLIT5.0')
			AND upstream.'eso pro catg' like 'RESP_TAB%'
			AND downstream.recipe_order > upstream.recipe_order
	WHERE
		upstream.complete = 1 and downstream_complete < 1 and (upstream.status = 'pass' or upstream.status is null)  and upstream_recipe in ('stare_std','nod_std','offset_std') and downstream_recipe not in ('stare_std','nod_std','offset_std')))
SELECT
	downstream_sof as sof,
	upstream_file as file,
	upstream_tag AS tag,
	upstream_filepath as filepath,
	downstream_complete as complete,
	upstream_recipe,
	downstream_recipe,
	*
FROM
	closest_upstream_match
WHERE (matched_upstream_ranked < 2)
ORDER BY
	upstream_recipe_order,
	'eso seq arm';
