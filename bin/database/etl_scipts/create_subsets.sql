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
    'TMP1',
    'Test dataset'
);

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
            'TMP1_'
            || t00.dataset_name
            || '_sscount'
            || floor(ROWNUM / 3) AS subset_name,
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
        FETCH FIRST 10 ROWS

only;

--association table

SELECT
    *
FROM
    tmp_sequencing_subsets;

INSERT INTO sequencing_subset (
    subset_group_id,
    subset_name
)
    SELECT DISTINCT
        subset_group_id,
        subset_name
    FROM
        tmp_sequencing_subsets;

INSERT INTO sequencing_subset_ebi_project_data (
    sequencing_subset_id,
    subset_name,
    ebi_project_data_id
)
    SELECT
        t0.id,
        t0.subset_name,
        t1.ebi_project_data_id
    FROM
        sequencing_subset t0
        INNER JOIN tmp_sequencing_subsets t1 ON t1.subset_group_id = t0.subset_group_id
                                                AND t0.subset_name = t1.subset_name;

DROP TABLE tmp_sequencing_subsets;

COMMIT;

-- check results

SELECT
    *
FROM
    sequencing_subset;

SELECT
    *
FROM
    sequencing_subset_ebi_project_data;

SELECT
    *
FROM
    sequencing_subset_machine;

-- truncate mapping

TRUNCATE TABLE sequencing_subset_machine;

TRUNCATE TABLE log_machine_processing;

SELECT
    *
FROM
    subset_group;

SELECT
    *
FROM
    sequencing_subset_machine;

SELECT
    *
FROM
    log_machine_processing
ORDER BY
    id DESC;

SELECT
    status,
    COUNT(*)
FROM
    sequencing_subset_machine
GROUP BY
    status;

SELECT
    machine_namess,
    status,
    COUNT(*)
FROM
    (
        SELECT
            sequencing_subset_id,
            id,
            machine_name,
            substr(machine_name,0,0) machine_namess,
            status
        FROM
            sequencing_subset_machine
    )
GROUP BY
    machine_namess,
    status;

SELECT
    *
FROM
    log_machine_processing
WHERE
    machine_name LIKE 'bio-a%'
ORDER BY
    id DESC;

SELECT
    *
FROM
    sequencing_subset_machine t0
    INNER JOIN log_machine_processing t1 ON t1.sequencing_subset_id = t0.sequencing_subset_id
WHERE
    status = 'FAILED'
ORDER BY
    t0.machine_name,
    t1.id DESC;
    
    
    
SELECT
    status, count(*)
FROM
    sequencing_subset_machine t0
GROUP BY status
    