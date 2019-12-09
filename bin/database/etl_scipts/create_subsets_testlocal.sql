-- SQL script for creating subsets stuff --

-- cleaning part
TRUNCATE TABLE sequencing_subset_machine CASCADE;
TRUNCATE TABLE sequencing_subset_ebi_project_data CASCADE;
TRUNCATE TABLE sequencing_subset CASCADE;
TRUNCATE TABLE subset_group CASCADE;
DELETE FROM sequencing_subset purge;
DELETE FROM subset_group purge;

COMMIT;

INSERT INTO subset_group (
    base_name,
    description
) VALUES (
    'T0',
    'Test dataset for debugging the massive concurent access errror'
);

select * from subset_group;


COMMIT;

SELECT
    *
FROM
    subset_group;
--params: subset size, lets say its 3
-- currently only fetchin 10
-- exluding the HMP data
-- filtering to stool samples
-- act_group_id: 

DROP TABLE tmp_sequencing_subsets;

CREATE TABLE tmp_sequencing_subsets
    AS
        SELECT
            1 subset_group_id,
            'T0_'
            || t00.dataset_name
            || '_sscount'
            || floor(ROWNUM / 10) AS subset_name,
            t00.ebi_project_data_id
        FROM
            (
                SELECT
                    t2.id ebi_project_data_id,
                    t2.dataset_name
                FROM
                    metadata t0
                    INNER JOIN metadata_ebi_project_data t1 ON t1.metadata_id = t0.id
                    INNER JOIN ebi_project_data t2 ON t2.id = t1.ebi_project_data_id
                WHERE
                    t0.body_site = 'stool'
                    AND   t0.dataset_id != 6
                ORDER BY
                    t2.dataset_id,
                    t2.id
            ) t00
--        FETCH FIRST 10 ROWS
--only;
;
--association table
SELECT * FROM tmp_sequencing_subsets;



insert into SEQUENCING_SUBSET
(subset_group_id, subset_name)
SELECT distinct subset_group_id, subset_name
from tmp_sequencing_subsets;

INSERT INTO SEQUENCING_SUBSET_EBI_PROJECT_DATA
(SEQUENCING_SUBSET_ID, SUBSET_NAME,EBI_PROJECT_DATA_ID)
SELECT 
t0.id,
t0.SUBSET_NAME,
t1.EBI_PROJECT_DATA_ID
FROM SEQUENCING_SUBSET t0
INNER JOIN tmp_sequencing_subsets t1
ON t1.subset_group_id = t0.subset_group_id 
and t0.subset_name = t1.subset_name
;
DROP TABLE tmp_sequencing_subsets;

commit;

-- check results
SELECT * FROM SEQUENCING_SUBSET ORDER BY ID desc;
SELECT * FROM SEQUENCING_SUBSET_EBI_PROJECT_DATA;
SELECT * FROM SEQUENCING_SUBSET_machine;

-- truncate mapping
--TRUNCATE TABLE SEQUENCING_SUBSET_machine;
--TRUNCATE TABLE LOG_MACHINE_PROCESSING;

SELECT * FROM SEQUENCING_SUBSET_machine ORDER BY ASSIGNING_TIME DESC;
SELECT * FROM LOG_MACHINE_PROCESSING order by id desc;


SELECT * FROM LOG_MACHINE_PROCESSING
WHERE machine_name = 'bioinfo275'
ORDER BY id
;

SELECT 
id,
machine_name,
SEQUENCING_SUBSET_id,
log_entry,
log_date,
TO_CHAR(log_date, 'YYYY-MON-DD HH24:MI:SS') log_date2
FROM LOG_MACHINE_PROCESSING 
WHERE 
--machine_name = 'bioinfo_515'
log_entry LIKE '%Error%'
ORDER BY log_date2 desc;
;
SELECT * FROM subset_group;


-- debugging: bioinfo_987

SELECT *
FROM LOG_MACHINE_PROCESSING
--WHERE machine_name = 'bioinfo_987'
ORDER BY ID desc;

SELECT *
FROM SEQUENCING_SUBSET_machine
WHERE SEQUENCING_SUBSET_id = 504



-- 33867	bio-node-2_163	504	Error, [Errno 17] File exists: '/gfs/computations/phage/results/RPH/T0_BritoIL_2016_sscount254_results/metaphlan'	26-MAR-18	2018-MAR-26 18:19:58

SELECT *
FROM LOG_MACHINE_PROCESSING
WHERE 
machine_name = 'bio-node-2_163'
and 
SEQUENCING_SUBSET_id = 504
ORDER BY ID desc;

RENAME COLUMN LOG_MACHINE_PROCESSING.SUBSET_id TO SEQUENCING_SUBSET_id
;
alter table LOG_MACHINE_PROCESSING 
rename column SUBSET_id to SEQUENCING_SUBSET_id;


